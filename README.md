# delim-SOM
Using multi-layer Kohonen Self-Organizing Maps ("SuperSOMs") to delimit species and produce integrative taxonomies using Unsupervised Machine Learning (UML).

This repository expands the use of Kohonen maps as described in Pyron et al. (2023). It uses multi-layer Self-Organizing Maps ("SuperSOMs") in the R package 'kohonen' to delimit species based on allelic, spatial, climatic, and phenotypic data.

The contribution of each layer to the final model output is recorded, along with the clustering assignment of each individual over multiple learning replicates. The results therefore mirror a 'STRUCTURE'-type analysis including individual ancestry coefficients, but represent an unified delimitation model that incorporates various dimensions of ecological and evolutionary divergence for integative taxonomy. 

# Example - _Desmognathus monticola_, the Seal Salamander

We provide a sample dataset and analysis for Seal Salamanders (_Desmognathus monticola_), which represent two species in the Appalachian Mountains and Piedmont of the eastern United States, based on four datasets comprising a SNP matrix from Genotype-By-Sequencing (GBS) analysis, lat/long/elevation (xyz), environment (climate, hydrology, and ecoregion), and phenotype (17 linear morphometric measurements).

The genetic, spatial, and environmental data come from 71 individuals from 71 sites, while the phenotypic data are expanded to include 156 specimens from those sites, with the mean of each measurement taken after size correction. The allele frequencies come from a GBS matrix of 7,809 SNPs, trimmed to 80% completeness.

The climate variables are Level IV Ecoregion (https://www.epa.gov/eco-research/level-iii-and-iv-ecoregions-continental-united-states), HUC4 watershed (https://www.usgs.gov/national-hydrography/watershed-boundary-dataset), ENVIREM - monthCountByTemp10 (https://envirem.github.io/), and WorldClim - BIO15 (https://www.worldclim.org/data/bioclim.html).

The phenotype variables are 17 linear morphometric measurements to 0.01mm precision: SVL (snout-vent length), TL (tail length), AG (axilla-groin length), CW (chest width), FL (femur length [rather than hindlimb length]), HL (humerus length [rather than forelimb length]), SG (snout-gular length), TW(tail width at rear of vent), TO (length of third toe), FI (length of third finger), HW (head width), ED (eye diameter), IN (internarial distance), ES (eye-snout distance), ON (orbito-narial distance), IO (inter-orbital distance), and IC (inter-canthal distance). Here, I size-correct these by SLV using pooled groups ("population2") in 'GroupStruct' (Chan and Grismer 2021, 2022): https://github.com/chankinonn/GroupStruct, then take the mean by site.



# References

Chan, K.O. & Grismer, L. L. (2021). A standardized and statistically defensible framework for quantitative morphological analyses in taxonomic studies. Zootaxa, 5023: 293-300.

Chan, K.O. and Grismer, L.L., 2022. GroupStruct: an R package for allometric size correction. Zootaxa, 5124(4), pp.471-482.

Pyron, R.A., O’Connell, K.A., Duncan, S.C., Burbrink, F.T. and Beamer, D.A., 2023. Speciation hypotheses from phylogeographic delimitation yield an integrative taxonomy for Seal Salamanders (Desmognathus monticola). Systematic Biology, 72(1), pp.179-197.
