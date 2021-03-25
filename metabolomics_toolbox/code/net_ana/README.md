## Uncovering in vivo metabolic associations from time series: an empirical network approach

This project provides solutions for analysis of time-series metabolic dynamics. Time dynamics were visualized by dimensionality reduction through FDA-PCA (functional data analysis principle component analysis). Functional states were found by network construction and clustering by CausalKinetiX and community clustering.

./src contains functions, ./scripts contains working scripts, ./tests contains test for functions.


To run the MATLAB workflow, you need to git the and add to MATLAB path for following repositories:

1. [Metabolomics toolbox](https://github.com/artedison/Edison_Lab_Shared_Metabolomics_UGA) is the metabolomics toolbox of [Edison Lab](http://edison.ccrc.uga.edu) and contains many useful MATLAB functions for metabolic data preprocessing and statistic analysis
2. [fda_learn](https://github.com/mikeaalv/fda_learn) is a wrapper for functions in functional data analysis library [fdaM](http://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/). fdaM library also need to be installed for the program to run properly.

To run the R workflow:

1. please install related packages: stringr, magrittr, R.matlab, RCy3, httr, jsonlite, reticulate, testthat, deSolve, ggplot2, dplyr, ComplexHeatmap, CausalKinetiX, foreach, doMC, igraph, reshape2
2. Cytoscape need to be installed for some R script.
3. Conda environment need to be constructed for network_clust_annotation.R
4. HPC environment is needed for simulation scripts.
