#!/bin/bash

#SBATCH -t 1-00:00:00
#SBATCH --job-name=fastsimcoal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=10GB
#SBATCH --partition=normal
#SBATCH --mail-type=end
#SBATCH --mail-user=daniel.schoenberger@uky.edu
#SBATCH --account=coa_jdu282_uksr
#SBATCH --error=SLURM_JOB_%j.err
#SBATCH --output=SLURM_JOB_%j.out

set -euo pipefail
shopt -s nullglob


################################################################################
#### Settings
################################################################################
WORKDIR="/pscratch/jdu282_uksr/Daniel_Hemi/fastsimcoal2"
IMG="/share/singularity/images/ccs/fastsimcoal2/28/mcc-fastsimcoal2-28-rocky9.sinf"
FSC="singularity exec ${IMG} fsc28"

# Biological / simulation settings
N_TDIV=301
MAX_TDIV=250000
CONTACT_FRACTION=0.5
MIGS=("0" "1e-6" "4e-6" "7e-6")

# Paired seed settings
FSC_SEED_BASE=1
DOWNSAMPLE_SEED_BASE=2

# Population and sampling settings
NPOP_HAPLOID=500000
SAMPLE_SIZE=72

# Sequence simulation settings
N_CHROMS=200000
N_LINKAGE_BLOCKS=1
N_LOCI_PER_BLOCK=1
RECOMB_RATE=0
MUT_RATE=1e-8

# VCF downsampling settings
TARGET_SNPS=3000
REPS=1

# Directory
OUTDIR="${WORKDIR}/vcf_outfiles"


################################################################################
#### Helper functions
################################################################################
write_tpl_file() {
  local tpl="${WORKDIR}/simulation.tpl"

  cat > "${tpl}" <<EOF2
//Number of population samples (demes)
2
//Population effective sizes (haploid number of genes)
${NPOP_HAPLOID}
${NPOP_HAPLOID}
//Sample sizes
${SAMPLE_SIZE}
${SAMPLE_SIZE}
//Growth rates
0
0
//Number of migration matrices : 0 means no migration
2
//Migration matrix 0: active migration
0 MIG
MIG 0
//Migration matrix 1: no migration
0 0
0 0
//Historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix
2 historical event
TCONTACT 0 0 0 1 0 1
TDIV 0 1 1 1 0 1
//Number of independent chromosomes and chromosome structure flag
${N_CHROMS} 0
//Number of contiguous linkage blocks
${N_LINKAGE_BLOCKS}
//Per block: data type, number of loci, rec rate, mut rate
DNA ${N_LOCI_PER_BLOCK} ${RECOMB_RATE} ${MUT_RATE}
EOF2
}

make_tdiv_values() {
  local tdiv_values_file="$1"

  awk -v n="${N_TDIV}" \
      -v max="${MAX_TDIV}" '
  BEGIN {
    if (n < 1) {
      print "ERROR: N_TDIV must be at least 1" > "/dev/stderr"
      exit 1
    }

    if (n == 1) {
      print 0
      exit 0
    }

    for (i = 0; i < n; i++) {
      t = int((i * max / (n - 1)) + 0.5)
      print t
    }
  }' | awk '!seen[$0]++' > "${tdiv_values_file}"

  local n_lines
  n_lines=$(wc -l < "${tdiv_values_file}")

  if [[ "${n_lines}" -ne "${N_TDIV}" ]]; then
    echo "ERROR: expected ${N_TDIV} TDIV values but found ${n_lines}" >&2
    exit 1
  fi
}

make_def() {
  local def="${WORKDIR}/simulation.def"
  local tdiv_values_file="${WORKDIR}/TDIV_values.txt"

  make_tdiv_values "${tdiv_values_file}"

  {
    echo "TDIV TCONTACT MIG"
    for mig in "${MIGS[@]}"; do
      while read -r tdiv; do
        tcontact="$(awk -v t="${tdiv}" -v contact_fraction="${CONTACT_FRACTION}" 'BEGIN { printf "%.0f", t * contact_fraction }')"
        printf "%s %s %s\n" "${tdiv}" "${tcontact}" "${mig}"
      done < "${tdiv_values_file}"
    done
  } > "${def}"
}

make_tmp_tpl() {
  local in_tpl="$1"
  local out_tpl="$2"
  local nchr="$3"

  awk -v nchr="${nchr}" '
  BEGIN { replace_next = 0 }
  {
    if ($0 ~ /\/\/Number of independent chromosomes and chromosome structure flag/) {
      print
      replace_next = 1
      next
    }

    if (replace_next == 1 && $0 !~ /^\/\//) {
      print nchr " 0"
      replace_next = 0
      next
    }

    print
  }' "${in_tpl}" > "${out_tpl}"
}

sanitize_tdiv() {
  local x="$1"
  x="${x%.000000}"
  x="${x//-/}"
  echo "${x}"
}

sanitize_mig() {
  local x="$1"
  x="${x%.000000}"
  echo "${x}"
}

gen_to_vcf() {
  local gen_file="$1"
  local vcf_file="$2"

  awk '
  BEGIN { FS="[[:space:]]+"; OFS="\t" }

  function gtcode(x, y) {
    y = x
    gsub(/\r/, "", y)
    if (y == "0") return "0/0"
    if (y == "1") return "0/1"
    if (y == "2") return "1/1"
    return "./."
  }

  NR == 1 {
    expected_nf = NF
    print "##fileformat=VCFv4.2"
    print "##source=fastsimcoal2_gen_to_vcf"
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"

    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    for (i = 5; i <= expected_nf; i++) {
      gsub(/\r/, "", $i)
      header = header OFS $i
    }
    print header
    next
  }

  NF < 4 { next }
  $1 !~ /^[0-9]+$/ { next }
  $2 !~ /^[0-9]+$/ { next }

  {
    printf "%s\t%s\t.\t%s\t%s\t.\tPASS\t.\tGT", $1, $2, $3, $4
    for (i = 5; i <= expected_nf; i++) {
      val = (i <= NF ? $i : "")
      printf OFS gtcode(val)
    }
    printf "\n"
  }' "${gen_file}" > "${vcf_file}"
}

downsample_vcf_to_target() {
  local in_vcf="$1"
  local out_vcf="$2"
  local target="$3"
  local seed="$4"

  local nvar
  nvar=$(grep -vc '^#' "${in_vcf}" || true)

  if [[ "${nvar}" -lt "${target}" ]]; then
    echo "ERROR: ${in_vcf} has only ${nvar} SNPs, fewer than target ${target}" >&2
    return 1
  fi

  {
    grep '^##' "${in_vcf}"
    grep '^#CHROM' "${in_vcf}"
    grep -v '^#' "${in_vcf}" | \
      awk -v seed="${seed}" 'BEGIN { srand(seed) } { print rand() "\t" $0 }' | \
      sort -k1,1n | \
      awk -v target="${target}" 'NR <= target { sub(/^[^\t]*\t/, ""); print }' | \
      sort -k1,1n -k2,2n
  } > "${out_vcf}"
}


################################################################################
#### Clean previous results
################################################################################
cd "${WORKDIR}"

rm -rf "${WORKDIR}/vcf_outfiles" \
       "${WORKDIR}"/tmp_sim[0-9]*_tdiv*_mig* \
       "${WORKDIR}"/sim[0-9]*_tdiv*_mig*

rm -f "${WORKDIR}/simulation.def" \
      "${WORKDIR}/simulation.tpl" \
      "${WORKDIR}/TDIV_values.txt" \
      "${WORKDIR}/TDIV_seeds.txt"

mkdir -p "${OUTDIR}"


################################################################################
#### Rebuild template, evenly spaced TDIV values, and parameter table
################################################################################
write_tpl_file
make_def


################################################################################
#### Check settings
################################################################################
if [[ "${REPS}" -ne 1 ]]; then
  echo "ERROR: this paired-seed script is written for REPS=1" >&2
  exit 1
fi

awk -v fsc_seed_base="${FSC_SEED_BASE}" \
    -v downsample_seed_base="${DOWNSAMPLE_SEED_BASE}" '
BEGIN {
  OFS = "\t"
  print "tdiv_index", "tdiv", "fsc_seed", "downsample_seed"
}
{
  print NR, $1, fsc_seed_base + $1, downsample_seed_base + $1
}' "${WORKDIR}/TDIV_values.txt" > "${WORKDIR}/TDIV_seeds.txt"


################################################################################
#### Run fastsimcoal2, convert to VCF, and downsample to target SNP number
################################################################################
echo "===== Running secondary gene flow simulations ====="

TPL="${WORKDIR}/simulation.tpl"

while read -r tdiv; do
  tdiv_tag="$(sanitize_tdiv "${tdiv}")"

  fsc_seed=$((FSC_SEED_BASE + tdiv))
  downsample_seed=$((DOWNSAMPLE_SEED_BASE + tdiv))

  for mig_index in "${!MIGS[@]}"; do
    mig="${MIGS[$mig_index]}"
    sim_tag=$((mig_index + 1))
    mig_tag="$(sanitize_mig "${mig}")"

    run_name="sim${sim_tag}_tdiv${tdiv_tag}_mig${mig_tag}"
    TMPDIR_RUN="${WORKDIR}/tmp_${run_name}_${SLURM_JOB_ID:-$$}"
    TMP_TPL="${TMPDIR_RUN}/${run_name}.tpl"
    TMP_DEF="${TMPDIR_RUN}/${run_name}.def"

    mkdir -p "${TMPDIR_RUN}"

    make_tmp_tpl "${TPL}" "${TMP_TPL}" "${N_CHROMS}"

    tcontact="$(awk -v t="${tdiv}" -v contact_fraction="${CONTACT_FRACTION}" 'BEGIN { printf "%.0f", t * contact_fraction }')"

    {
      echo "TDIV TCONTACT MIG"
      printf "%s %s %s\n" "${tdiv}" "${tcontact}" "${mig}"
    } > "${TMP_DEF}"

    echo "Running ${run_name} with fastsimcoal2 seed ${fsc_seed}"

    eval "${FSC} -t ${TMP_TPL} -f ${TMP_DEF} -n ${REPS} -G -g --seed ${fsc_seed}"

    GEN_FILES=( "${WORKDIR}/${run_name}"/*.gen "${TMPDIR_RUN}/${run_name}"/*.gen )

    if [[ ${#GEN_FILES[@]} -eq 0 ]]; then
      echo "ERROR: no .gen files found for ${run_name}" >&2
      exit 1
    fi

    if [[ ${#GEN_FILES[@]} -ne 1 ]]; then
      echo "ERROR: expected 1 .gen file for ${run_name}, found ${#GEN_FILES[@]}" >&2
      exit 1
    fi

    gen="${GEN_FILES[0]}"
    tmp_vcf="${TMPDIR_RUN}/${run_name}.tmp.vcf"
    final_vcf="${OUTDIR}/${run_name}.vcf"

    echo "Converting ${gen} -> temporary VCF"
    gen_to_vcf "${gen}" "${tmp_vcf}"

    echo "Downsampling ${tmp_vcf} -> ${final_vcf} with seed ${downsample_seed}"
    downsample_vcf_to_target "${tmp_vcf}" "${final_vcf}" "${TARGET_SNPS}" "${downsample_seed}"

    rm -rf "${WORKDIR}/${run_name}" "${TMPDIR_RUN}"
  done
done < "${WORKDIR}/TDIV_values.txt"

rm -rf "${WORKDIR}"/tmp_sim[0-9]*_tdiv*_mig*
rm -rf "${WORKDIR}"/sim[0-9]*_tdiv*_mig*

echo "Done."
echo "Final VCFs are in: ${OUTDIR}"
echo "TDIV values are in: ${WORKDIR}/TDIV_values.txt"
echo "TDIV seed table is in: ${WORKDIR}/TDIV_seeds.txt"
