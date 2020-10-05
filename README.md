# ChIPanalyser - Analysis Code

This repository contains all code used for the analysis presented in the following manuscript.


## Analysis 

This folder contains source code for all analysis presented in the paper. This includes processing of overlapping regions, combining of ChIP data sets, extraction of relevant data for plot production.

## Batch Run 

This folder contains the main script used for this analysis. The data input table contains all the parameters that will be parse to R and ChIPanalyser. 

It should be noted that the analysis was performed on a computer cluster and each "run" (different set of combinations) were run in parallel. It is for this reason that we only provide a subset of the analysis as it is just a question of changing parameters in this table when necessary.


NOTE: The path towards the data input file should be changed in both the `submit.sh` file as well as the `ChIPanalPerformAnalysis.sh` file.  


## Plotting 

This folder contains code specific to plots produced in the manuscript. In order to build manuscript worthy plots, we plotted data separately. 