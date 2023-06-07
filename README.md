# delim-SOM
Using multi-layer Kohonen Self-Organizing Maps ("SuperSOMs") to delimit species and produce integrative taxonomies using Unsupervised Machine Learning (UML).

This repository expands the use of Kohonen maps as described in Pyron et al. (2023). It uses multi-layer Self-Organizing Maps ("SuperSOMs") in the R package 'kohonen' to delimit species based on allelic, spatial, climatic, and phenotypic data.

We provide a sample dataset and analysis for Seal Salamanders (Desmognathus monticola), which represent two species in the Appalachian Mountains and Piedmont of the eastern United States, based on four datasets comprising a SNP matrix from Genotype-By-Sequencing (GBS) analysis, lat/long/elevation (xyz), environment (climate, hydrology, and ecoregion), and phenotype (17 linear morphometric measurements).

The contribution of each layer to the final model output is recorded, along with the clustering assignment of each individual over multiple learning replicates. The results therefore mirror a 'STRUCTURE'-type analysis including individual ancestry coefficients, but represent an unified delimitation model that incorporates various dimensions of ecological and evolutionary divergence for integative taxonomy. 

# Example - Desmognathus monticola, the Seal Salamander

# References

Pyron, R.A., O’Connell, K.A., Duncan, S.C., Burbrink, F.T. and Beamer, D.A., 2023. Speciation hypotheses from phylogeographic delimitation yield an integrative taxonomy for Seal Salamanders (Desmognathus monticola). Systematic Biology, 72(1), pp.179-197.
