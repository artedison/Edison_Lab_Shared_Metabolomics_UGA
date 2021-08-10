## Uncovering in vivo metabolic associations from time series: an empirical network approach

This project provides solutions for analysis of time-series metabolic dynamics. Time dynamics were visualized by dimensionality reduction through FDA-PCA (functional data analysis principle component analysis). Functional states were found by network construction and clustering by CausalKinetiX and community clustering.

./src contains functions, ./scripts contains working scripts, ./tests contains test for functions.


To run the MATLAB workflow, you need to git the and add to MATLAB path for following repositories:

1. [Metabolomics toolbox](https://github.com/artedison/Edison_Lab_Shared_Metabolomics_UGA) is the metabolomics toolbox of [Edison Lab](http://edison.ccrc.uga.edu) and contains many useful MATLAB functions for metabolic data preprocessing and statistic analysis
2. [fda_learn](https://github.com/mikeaalv/fda_learn) is a wrapper for functions in functional data analysis library [fdaM](http://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/). fdaM library also need to be installed for the program to run properly.
3. [subplot_tight](https://www.mathworks.com/matlabcentral/fileexchange/30884-controllable-tight-subplot) is needed for some process. 

To run the R workflow:

1. please install related packages: stringr(1.4.0), magrittr(1.5), R.matlab(3.6.2), RCy3(2.2.9), httr(1.4.1), jsonlite(1.6), reticulate(1.18-9000), testthat(2.2.1), deSolve(1.28), ggplot2(3.2.1), dplyr(0.8.5), ComplexHeatmap(1.18.1), CausalKinetiX(0.2.3), foreach(1.4.7), doMC(1.3.6), igraph(1.2.5), reshape2(1.4.3)
2. Cytoscape (3.8.0) need to be installed for some R script. Clustermaker (1.3.1) is needed for some pipeline.
3. Conda (4.7.11) environment need to be constructed for network_clust_annotation.R
4. HPC environment is needed for simulation scripts.



