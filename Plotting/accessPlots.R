################################################################################
######################## Access plots ##########################################
################################################################################
direc <- getwd()

library(BSgenome.Dmelanogaster.UCSC.dm6)
library(BSgenome)
library(RcppRoll)
library(parallel)
library(GenomicRanges)
library(ROCR)
## data load
kc<-get(load("/home/pm16057/DNAaccess/cellAccess/Kc_DHS_005.Rda"))
s2<-get(load("/home/pm16057/DNAaccess/cellAccess/S2_DHS_005.Rda"))
bg3<-get(load("/home/pm16057/DNAaccess/cellAccess/BG3_DHS_005.Rda"))

embryo<-get(load("/home/pm16057/ChIPanalyser/ChIPanalyserFinal/ChIPanalyserTesting/RaduData/AccessOld.Rdata"))

DNA<-get(load("/home/pm16057/ChIPanalyser/ChIPanalyserFinal/performAnalysis/DNASequenceSet.Rda"))
DNA<-sum(width(DNA))
### Intersecting with themselves

kc<-intersect(kc,kc)
s2<-intersect(s2,s2)
bg3<-intersect(bg3,bg3)

embryo<-intersect(embryo,embryo)

#
AccessProp<-c("Kc167"=(sum(width(kc))/DNA)*100,
                                 "Embryo"=(sum(width(embryo))/DNA)*100,
                                 "S2"=(sum(width(s2))/DNA)*100,
                                 "BG3"=(sum(width(bg3))/DNA)*100)

pdf("AccessProp.pdf",height=7,width=7)
par(family="sans")
par(mar=c(8,8,4,4))
cols2<-c("#233142","#facf5a","#4f9da6","#ff5959")
barplot(AccessProp,col=cols2,cex.axis=1.5,cex.lab=1.5,las=2,cex=1.5)
title(ylab="% of Accessible DNA",cex.lab =1.8, line =6)

dev.off()
