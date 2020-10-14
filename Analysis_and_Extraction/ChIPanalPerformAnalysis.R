#!/usr/bin/Rscript3.5.0

### Loading Libraries and Scripts

direc<-getwd()

library(BSgenome.Dmelanogaster.UCSC.dm6)
library(BSgenome)
library(RcppRoll)
library(GenomicRanges)
library(ROCR)
library(ChIPanalyser)

source("DataHand.R")

args <- commandArgs(TRUE)
print(args)
# Parsing NULL's because R is doenst want to

if(grepl("NULL",args[7])){args7<-NULL}else{args7<-args[7]}

if(grepl("NULL",args[8])){args8<-NULL}else{args8<-get(load(args[8]))}

if(grepl("NULL",args[9])){args9<-NULL}else{args9<-as.numeric(args[9])}


if(grepl("NULL",args[10])){args10<-NULL}else{args10<-args[10]}

if(grepl("NULL",args[13])){args13<-NULL}else{args13<-args[13]}

args11<-as.numeric(args[11])
args12<-as.numeric(args[12])

args14<-as.character(args[14])
args15<-as.character(args[15])
## Adding Other methods 
args16<-c(as.character(args[16]),"AUC","recall","spearman")

## Building TF list for performAnalysis
TFList <- list(as.character(args[3]), as.character(args[4]))

## Loading Dna Accessibility RDA file 
if(!is.null(args7)){
    splits<-strsplit(args7, "/")[[1]]
    access<-splits[length(splits)]
    access<-strsplit(access,".Rda")[[1]]
} else {
    access<-"NULL"
}


## Setting up file names for outpur
if(is.null(args9)){
    name <- paste0(args[1],"_",args[2],"_step100_",access)
}else {
    name <- paste0(args[1],"_",args[2],"_reduce_top10opti",args9,"_",access)
}
print(name)
rm(loci)

# Change directory to output folder
setwd(paste0(direc,"/",args15))



#setwd(paste0(direc,""))
performAnalysis(TFList=TFList,
                ChIP=args[5],
                DNASequenceSet=args[6],
                Access=args7,
                setSequence=args8,
                reduce=args9,
                peaks=args10,
                tileSize=args11,
                cores=args12,
                filename=name,
                OP=args13,
                noiseFilter=args14,
                method=args16)
