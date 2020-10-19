# ChIPanalyser - CSBJ

The following repository contains R code related to the ChIPanalyser manuscript. It should be noted that given the wide range of parameters tested for the manuscript, we only present the working skeleton. Parameters should be changed accordingly. For batch submission - see section below. When batch submission was not necessary, R code was sourced directly through the command line. It should be noted that most of the code presented here only serves to handle large amount of data.

We highly recommend using the [ChIPanalyser Vignette](https://www.bioconductor.org/packages/release/bioc/html/ChIPanalyser.html) for a simplified ChIPanalyser pipeline. 

## Data Table \& Job Submission 
The bulk of the analysis was carried out of the University of Essex HPC (Ceres). In order to facilitate the analysis of many data sets simultaneously, we batch-submitted individual jobs using the following pipeline. 

1. First, we created a data input  table (space separated and no quotes) containing all arguments that will be parsed to our ChIPanalyser wrapper script (performAnalysis.R - see below). We provide a dummy Data input table (DataInputDHS.txt). Each column represents a separate argument and can be described as:
	* Cell line 
	* Data set name
	* Path to PFM/PWM file 
	* Format of PFM/PWM file 
	* Path to ChIP enrichment signal file (.bed, .wig,...)
	* Path to DNA Sequence - DNASequenceSet object extracted from the BSgenome package and saved as Rda file (optional)
	* Path to DNA accessibility data saved as Granges Object (.Rda) - if no accessibility provided, set column value to "NULL".
	* Path towards loci of interest saved as Granges Object (.Rda) - if no Loci of interest provided, set column value to "NULL".
	* Number of regions to select for the analysis (optional) - In current analysis, regions are extracted later - thus validation is carried out on all regions - set argument to "NULL".
	* Path towards peak location file (.bed, .gff,...) - if no Peaks provided, set column value to "NULL".

	* Size of genome bins ( default is set to 20kbp)
	* Number of cores to use 
	* Depreciated File naming argument 
	* Noise filter method to be applied to data
	* Path towards output folder 
	* Goodness of fit method to be used 

2. Once the data table has been saved, the path towards this file should be changed in both the `submit.sh` file and the `ChIPanalPerformAnalysis.sh` file. 

3. From the command line: 
``` 
./submit.sh

```

The `submit.sh` script will call the  `ChIPanalPerformAnalysis.sh` script for each row in the Data input table and submit each job individually.  `ChIPanalPerformAnalysis.sh` makes a call to  `ChIPanalPerformAnalysis.R` which is the R script used to perform our analysis. For details on the analysis  - see below.
	

## Analysis and Extraction

The Analysis and Extraction folder contains all scripts included in our analysis. 

### Main Functions 

* `DataHand.R`: R script containing all functions used to handle and parse data after analysis. These Functions mainly serve the purpose of handling many data sets simultaneously.
* `ChIPanalPerformAnalysis.R` : R script parsing argument from the submission files. 
* `performAnalysis.R` : wrapper functions to run full analysis. This script saves the output of the analysis in the form of `.Rda` files. These files contain all information (N - lambda - estimated profiles etc) required for further analysis and plotting. 

### Main analysis Runs when no batch submission was required 

* `HoxAnalysis.R`: R script for the Analysis of Hox TFs .
* `pwmMotif.R`: R script used to generate PWM Motifs plot.
* `rnaNew.R`: R script used in the RNA rescaling analysis. 
* `human.R`: R script used to run ChIPanalyser on human data sets (tool comparison was run on human data as not all tools provide support for _Drosophila_).
* `toopComp_train10_va20.R`: R script used to assess the performance of other tools compared to ChIPanalyser (NOTE: number of regions used for training/validation should be changed accordingly).
* `noiseFilter.R`: R script used to analyse the effect of our noise filtering methods. 
* `chromosomeWitholding.R`: R script used to validate the model between chromosomes.
* `toolProfileComp.R`: R script used to extract and produce ChIP like profiles for each tool.

### Data handling 

* `allData.R`: R script used to import saved .Rda files (output of `performAnalysis.R`) and extract relevant  validation data (saved as .Rda files).
* `allDataParsing.R` R script used to import and handle extracted validation data. The functions provided here only clean and reformat data for faster and automated plotting. 
* `dataConsist.R`: R script use to overlay optimal parameters from multiple data sets and plot the results in the form of a heat map. 

## Preprocessing 

* `QDA.R`: R script used to generate QDA accessibility data used in the Hox TF analysis.
* `peakOverlap.R`: R script used to generate loci of interest.
* `combi.R`: R script used to combine S2 CTCF ChIP data sets.



## Plotting 

The plotting folder contains the scripts used for plotting.

* `accessPlots.R`: R script used to produce difference in accessibility between cell lines. 
* `HoxPlotting.R`: R script used to produce Hox profile plots after running `HoxAnalysis.R`.
* `main_plots.R`: R script used to produce main manuscript plots (must run `allData.R`,`allDataParsing.R`,`dataConsist.R` beforehand).
* `methodPlots.R`: R script used to produce noise filter method plots and metric difference plots .
* `chromosomeWitholdingPlots.R`: R script used to plot validation of model between chromosomes ( must run `chromosomeWitholding.R` beforehand).
* `DataConsist.R`: R script used to create overlay matrices. 
* `toolCompPlots.R`: R script used to create plots comparing the performance of various tools and frameworks. 
* `rna.R`: R script used to produce the RNA rescaling plots (must run `rnaNew.R` beforehand).
* `justifPlots.R`: R script used to produce preliminary analysis and justification on Pearson drop metric drop.
* `pwmMotif.R`: R script used to produce PWM motif plots. 
* `toolProfilePlot.R`: R script used to produce profiles for each tool. 
