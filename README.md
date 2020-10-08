# ChIPanalyser - CSBJ

The following repository contain R code related to the ChIPanalyser manuscript. It should be noted than given the wide range of parameters tested for the manuscript. We only present the working skeleton. Parameters should be changed accordingly. For batch submission, see section below. When batch submission was not necessary, R code was sourced directly through the command line. It should be noted that the code presented here extends the analysis capabilities offered by the ChIPanalyser package. The code used here serves the purpose of automation of analysis for many data sets. 


## Data Table \& Job Submission 
The bulk of the analysis was carried out of the University of Essex cluster (Ceres). In order to facilitate the analysis of many data sets simultaneously, we batch submitted individual jobs using the following pipeline. 

1. First, we created a data input  table (space separated and no quotes) containing all arguments that will be parsed to our ChIPanalyser wrapper script (performAnalysis.R - see below). Each column, represents a separate argument and can be described as:
	* Cell line 
	* Data set name
	* Path to PFM/PWM file 
	* Format of PFM/PWM file 
	* Path to ChIP enrichment signal file (.bed, .wig,...)
	* Path to DNA Sequence - DNASequenceSet object extracted from the BSgenome package and saved as Rda file (optional)
	* Path to DNA accessibility data saved as Granges Object (.Rda)
	* Path towards loci of interest saved as Granges Object (.Rda)
	* Number of regions to select for the analysis (optional)
	* Path towards peak location file (.bed, .gff,...)
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

The `submit.sh` script will call the  `ChIPanalPerformAnalysis.sh` script for each row in the Data input table and submit each job individually.  `ChIPanalPerformAnalysis.sh` makes a call to  `ChIPanalPerformAnalysis.R` which is the R script used to perform our analysis. For detail on the analysis see below.
	

## Analysis 

The analysis folder contains all scripts include in our analysis. 

* `ChIPanalPerformAnalysis.R` : R script parsing argument from the submission files. 
* `performAnalysis.R` : wrapper functions to run full analysis. This script saves the output of the analysis in the form of Rda files. These files contain all information (N - lambda - estimated profiles etc) required for further analysis and plotting. 
* `HoxAnalysis.R`: R script for the Analysis of Hox TFs 
* `pwmMotif.R`: R script used to generate PWM Motifs plot
* `QDA.R`: R script used to generate QDA accessibility data used in the Hox TF analysis.
* `peakOverlap.R`: R script used to generate loci of interest.
* `noiseFilter.R`: R script used to analyse the effect of our noise filtering methods. 
* `extractParam.R`: R script containing functions as well as code used to extract relevant optimal parameters and export in a Latex format style. 
* `allData.R`: R code used to import saved .Rda files (output of `performAnalysis.R`) and extract relevant  validation data (saved as .Rda files)
* `allDataParsing.R` R code used to import and handle extracted validation data. The functions provided here only clean and reformat data for faster and automated plotting. 
* `combi.R`: R code used to combine all S2 data sets.
* `dataConsist.R`: R code use to overlay optimal parameters from multiple data sets and plot the results in the form of a heat map. 
* `rnascaling.R`: R code used in the RNA rescaling analysis. This section of the code also contains the plotting.
## Plotting 

The plotting folder contains the scripts used in order to produce plots.

* 