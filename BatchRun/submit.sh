#!/bin/bash
# File: subit.sh
## ChIPanalyser Perform Analysis ##

##### All files listed here are the core files containing path to files used for analysis ######

## Step 1: assessing number of files that will analysed

files=$(wc -l DataInputNULLgeometric.txt | awk '{print $1}')


## Creatin looping vector (skipping colnames)

row=$(seq 1 1 $files)


## Looping over files
## NOTE: argument MUST BE in the CORRECT ORDER! See: performAnalysis.R

for i in $row
    do
    qsub ./ChIPanalPerformAnalysis.sh $i
    echo $1
done
