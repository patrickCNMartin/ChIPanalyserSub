############################# Data Handling Functions     ######################
direc<-getwd()

library(BSgenome.Dmelanogaster.UCSC.dm6)
library(BSgenome)
library(RcppRoll)
library(GenomicRanges)
library(ROCR)
library(ChIPanalyser)

###### loading ChIP data script####


chipLoading<-function(ChIP){
    print("Loading ChIP data")
    
    if(class(ChIP)!="data.frame"){
        if(grepl(".Rda",ChIP)){
            ChIP<-get(load(ChIP))
        } else if(grepl(".bed",ChIP)){
            # Tab Format
            
            if(length(grep(x=readLines(ChIP,5),pattern="\t"))>0){
                if(length(grep(x=readLines(ChIP,5),pattern="\t"))!=5){
                    buffer<-5-length(grep(x=readLines(ChIP,5),pattern="\t"))
                    ChIP<-read.delim(ChIP, header=F,sep="\t", stringsAsFactors=F, skip=buffer)
                } else {
                    ChIP<-read.delim(ChIP, header=F,sep="\t", stringsAsFactors=F)
                }
                # if "chr" is missing
                if(length(grep(x=ChIP[,1], pattern="chr"))==0 ){
                    ChIP[,1]<-paste0("chr",ChIP[,1])
                }
                # bed format test
                if(ncol(ChIP)<4){
                    stop(paste0(deparse(substitute(ChIP)),
                                " is not in bed file Format. Please convert to bed/bedGraph file format"))}
                # Remove extra columns from format change
                if(ncol(ChIP)>4){
                    ChIP<-ChIP[,-4]
                }
            } else if(length(grep(x=readLines(ChIP,5),pattern=" "))>0){
                #spacing format
                if(length(grep(x=readLines(ChIP,5),pattern="\t"))!=5){
                    buffer<-5-length(grep(x=readLines(ChIP,5),pattern="\t"))
                    ChIP<-read.delim(ChIP, header=F,sep="\t", stringsAsFactors=F, skip=buffer)
                } else {
                    ChIP<-read.delim(ChIP, header=F,sep="\t", stringsAsFactors=F)
                }
                # if "chr" is missing
                if(length(grep(x=ChIP[,1], pattern="chr"))==0 ){
                    ChIP[,1]<-paste0("chr",ChIP[,1])
                }
                # bed format test
                if(ncol(ChIP)<4){
                    stop(paste0(deparse(substitute(CHIP)),
                                " is not in bed file Format. Please convert to bed/bedGraph file format"))}
                # Remove extra columns from format change
                if(ncol(ChIP)>4){
                    ChIP<-ChIP[,-4]
                }
            } else if(length(grep(x=readLines(ChIP,5),pattern=","))>0){
                # CSV format
                if(length(grep(x=readLines(ChIP,5),pattern="\t"))!=5){
                    buffer<-5-length(grep(x=readLines(ChIP,5),pattern="\t"))
                    ChIP<-read.delim(ChIP, header=F,sep="\t", stringsAsFactors=F, skip=buffer)
                } else {
                    ChIP<-read.delim(ChIP, header=F,sep="\t", stringsAsFactors=F)
                }
                # if "chr" is missing
                if(length(grep(x=ChIP[,1], pattern="chr"))==0 ){
                    ChIP[,1]<-paste0("chr",ChIP[,1])
                }
                # bed format test
                if(ncol(ChIP)<4){
                    stop(paste0(deparse(substitute(CHIP)),
                                " is not in bed file Format. Please convert to bed/bedGraph file format"))}
                # Remove extra columns from format change
                if(ncol(ChIP)>4){
                    ChIP<-ChIP[,-4]
                    
                }
            }
            
        } else if(grepl(".bedGraph",ChIP)){
            if(length(grep(x=readLines(ChIP,5),pattern="\t"))>0){
                if(length(grep(x=readLines(ChIP,5),pattern="\t"))!=5){
                    buffer<-5-length(grep(x=readLines(ChIP,5),pattern="\t"))
                    ChIP<-read.delim(ChIP, header=F,sep="\t", stringsAsFactors=F, skip=buffer)
                } else {
                    ChIP<-read.delim(ChIP, header=F,sep="\t", stringsAsFactors=F)
                }
                # if "chr" is missing
                if(length(grep(x=ChIP[,1], pattern="chr"))==0 ){
                    ChIP[,1]<-paste0("chr",ChIP[,1])
                }
                # bed format test
                if(ncol(ChIP)<4){
                    stop(paste0(deparse(substitute(ChIP)),
                                " is not in bed file Format. Please convert to bed/bedGraph file format"))}
                # Remove extra columns from format change
                if(ncol(ChIP)>4){
                    ChIP<-ChIP[,-4]
                }
            } else if(length(grep(x=readLines(ChIP,5),pattern=" "))>0){
                #spacing format
                if(length(grep(x=readLines(ChIP,5),pattern=" "))!=5){
                    buffer<-5-length(grep(x=readLines(ChIP,5),pattern=" "))
                    ChIP<-read.delim(ChIP, header=F,sep=" ", stringsAsFactors=F, skip=buffer)
                } else {
                    ChIP<-read.delim(ChIP, header=F,sep=" ", stringsAsFactors=F)
                }
                # if "chr" is missing
                if(length(grep(x=ChIP[,1], pattern="chr"))==0 ){
                    ChIP[,1]<-paste0("chr",ChIP[,1])
                }
                # bed format test
                if(ncol(ChIP)<4){
                    stop(paste0(deparse(substitute(CHIP)),
                                " is not in bed file Format. Please convert to bed/bedGraph file format"))}
                # Remove extra columns from format change
                if(ncol(ChIP)>4){
                    ChIP<-ChIP[,-4]
                }
            } else if(length(grep(x=readLines(ChIP,5),pattern=","))>0){
                # CSV format
                if(length(grep(x=readLines(ChIP,5),pattern=","))!=5){
                    buffer<-5-length(grep(x=readLines(ChIP,5),pattern=","))
                    ChIP<-read.csv(ChIP, header=F, stringsAsFactors=F, skip=buffer)
                } else {
                    ChIP<-read.csv(ChIP, header=F, stringsAsFactors=F)
                }
                # if "chr" is missing
                if(length(grep(x=ChIP[,1], pattern="chr"))==0 ){
                    ChIP[,1]<-paste0("chr",ChIP[,1])
                }
                # bed format test
                if(ncol(ChIP)<4){
                    stop(paste0(deparse(substitute(CHIP)),
                                " is not in bed file Format. Please convert to bed/bedGraph file format"))}
                # Remove extra columns from format change
                if(ncol(ChIP)>4){
                    ChIP<-ChIP[,-4]
                    
                }
            }
        } else if(grepl(".wig", ChIP)){
            if(length(grep(x=readLines(ChIP,5),pattern="\t"))>0){
                if(length(grep(x=readLines(ChIP,5),pattern="\t"))!=5){
                    buffer<-5-length(grep(x=readLines(ChIP,5),pattern="\t"))
                    ChIP<-read.delim(ChIP, header=F,sep="\t", stringsAsFactors=F, skip=buffer)
                } else {
                    ChIP<-read.delim(ChIP, header=F,sep="\t", stringsAsFactors=F)
                }
                # if "chr" is missing
                if(length(grep(x=ChIP[,1], pattern="chr"))==0 ){
                    ChIP[,1]<-paste0("chr",ChIP[,1])
                }
                # bed format test
                if(ncol(ChIP)<4){
                    stop(paste0(deparse(substitute(ChIP)),
                                " is not in bed file Format. Please convert to bed/bedGraph file format"))}
                # Remove extra columns from format change
                if(ncol(ChIP)>4){
                    ChIP<-ChIP[,-4]
                }
            } else if(length(grep(x=readLines(ChIP,5),pattern=" "))>0){
                #spacing format
                if(length(grep(x=readLines(ChIP,5),pattern=" "))!=5){
                    buffer<-5-length(grep(x=readLines(ChIP,5),pattern=" "))
                    ChIP<-read.delim(ChIP, header=F,sep=" ", stringsAsFactors=F, skip=buffer)
                } else {
                    ChIP<-read.delim(ChIP, header=F,sep=" ", stringsAsFactors=F)
                }
                # if "chr" is missing
                if(length(grep(x=ChIP[,1], pattern="chr"))==0 ){
                    ChIP[,1]<-paste0("chr",ChIP[,1])
                }
                # bed format test
                if(ncol(ChIP)<4){
                    stop(paste0(deparse(substitute(CHIP)),
                                " is not in bed file Format. Please convert to bed/bedGraph file format"))}
                # Remove extra columns from format change
                if(ncol(ChIP)>4){
                    ChIP<-ChIP[,-4]
                }
            } else if(length(grep(x=readLines(ChIP,5),pattern=","))>0){
                # CSV format
                if(length(grep(x=readLines(ChIP,5),pattern=","))!=5){
                    buffer<-5-length(grep(x=readLines(ChIP,5),pattern=","))
                    ChIP<-read.csv(ChIP, header=F, stringsAsFactors=F, skip=buffer)
                } else {
                    ChIP<-read.csv(ChIP, header=F, stringsAsFactors=F)
                }
                # if "chr" is missing
                if(length(grep(x=ChIP[,1], pattern="chr"))==0 ){
                    ChIP[,1]<-paste0("chr",ChIP[,1])
                }
                # bed format test
                if(ncol(ChIP)<4){
                    stop(paste0(deparse(substitute(CHIP)),
                                " is not in bed file Format. Please convert to bed/bedGraph file format"))}
                # Remove extra columns from format change
                if(ncol(ChIP)>4){
                    ChIP<-ChIP[,-4]
                    
                }
            }
        }
    }
    return(ChIP)
}


######################################
########### Annex functions ##########
######################################

peaksLoading<-function(x){
    # Laoding files based on files extension
    
    if(grepl(x=x,pattern=".bed")& !grepl(x=x,pattern=".gff")){
        
        x<-read.table(x, stringsAsFactors=F)
        if(length(grep(x=x[,1], pattern="chr"))==0){
            x[,1]<-paste0("chr",x[,1])
            x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,2],x[,3]))
        }else{
            x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,2],x[,3]))
        }
        
    }else if(grepl(x=x, pattern=".Rda")){
        x<-get(load(x))
    }else if(grepl(x=x, pattern=".gff3")){
        x<-read.table(x, skip=30,stringsAsFactors=F)
        if(length(grep(x=x[,1], pattern="chr"))==0){
            x[,1]<-paste0("chr",x[,1])
            x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,4],x[,5]))
        }else{
            x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,4],x[,5]))
        }
        
    } else if(grepl(x=x, pattern=".gff")){
        x<-read.table(x, skip=30,stringsAsFactors=F)
        if(length(grep(x=x[,1], pattern="chr"))==0){
            x[,1]<-paste0("chr",x[,1])
            x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,4],x[,5]))
        }else{
            x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,4],x[,5]))
        }
    } else{
        x<-read.table(x, header=T,stringsAsFactors=F)
        if(length(grep(x=x[,1], pattern="chr"))==0){
            x[,1]<-paste0("chr",x[,1])
            x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,2],x[,3]))
        }else{
            x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,2],x[,3]))
        }
    }
    return(x)
}




######################################
########## Perform analysis ##########
######################################

# Single Run analysis - looping should be done in bash so you can run thing fatser


performAnalysis<-function(TFList,ChIP,DNASequenceSet,
                          Access=NULL,setSequence=NULL,
                          reduce=NULL,peaks=NULL,tileSize=20000, cores=1,
                          filename=NULL,OP=NULL,noiseFilter="sigmoid",method="pearsonMean",normValues=NULL,withValidation=T){
    
    
    print("Loading DNA Access")
    if(!is.null(Access) & class(Access)!="GRanges"){
        Access<-get(load(Access))
        name <- as.character(seqnames(Access))
        Access <- GRanges(seqnames=name,
                          range= IRanges(start= start(Access), end=end(Access)))
    }
    ### extracting relavalnt data and loading it
    
    PO<-parameterOptions(noiseFilter=noiseFilter,stepSize=100)
    
    # DNASequenceSet
    print("Loading DNA Seq")
    if(class(DNASequenceSet)=="character" ){
        DNASequenceSet<-get(load(DNASequenceSet))
    }else if(class(DNASequenceSet)=="BSgenome"){
        DNASequenceSet<-getSeq(DNASequenceSet)
    }
    ### Loading TF PFM ###
    print("Loading PFM")
    if(class(TFList)!="list" & length(TFList)!=2){
        stop(paste0(deparse(substitute(TFList))," is not a List of 2 elements"))
    }
    if(grepl(".Rda",TFList[[1]])){
        PFM<-get(load(TFList[[1]]))
        GPP<-genomicProfiles(PFM=PFM,PFMFormat=TFList[[2]],
                             BPFrequency=DNASequenceSet,PWMThreshold=0.7)
    } else{
        GPP<-genomicProfiles(PFM=TFList[[1]],PFMFormat=TFList[[2]],
                             BPFrequency=DNASequenceSet,PWMThreshold=0.7)
    }
    
    
    ## precessing ChIP seq
    chipProfile <- processingChIP(profile=ChIP,loci=setSequence, reduce=reduce,
                                  peaks=NULL, chromatinState=NULL,
                                  parameterOptions=PO,cores=cores,normValues=normValues)
    if(withValidation){
        
        
        if(length(setSequence)<=10){
            TrainScore<-scores(chipProfile)[seq_along(setSequence)]
            ValidationScore<-scores(chipProfile)
            
            Trainloci<-loci(chipProfile)[seq_along(setSequence)]
            Validationloci<-loci(chipProfile)
        } else{
            TrainScore<-scores(chipProfile)[1:10]
            ValidationScore<-scores(chipProfile)
            
            Trainloci<-loci(chipProfile)[1:10]
            Validationloci<-loci(chipProfile)
        }
        
        
        
        
        .scores(chipProfile)<-TrainScore
        .loci(chipProfile)<-Trainloci
        
        
        # compute optimal
        optimal <- computeOptimal(GPP,DNASequenceSet,chipProfile,Access,PO,cores=cores)
        
        
        ### Computing Optimal Parameters ###
        
        if(!is.null(filename)){
            save(optimal,file=paste0(filename,noiseFilter,"_OptimalOutputTraining.Rda"))
            save(chipProfile, file=paste0(filename,noiseFilter,"_ChIPTraining.Rda"))
        } else {
            save(optimal,file=paste0(noiseFilter,"optimalOutputTraining.Rda"))
            save(chipProfile, file=paste0(noiseFilter,"ChIPTraining.Rda"))
        }
        
        ## optimal param for validation
        ## let assume you have more than one
        
        for(meth in method){
            
            param<-optimal[[1]][[1]][[meth]]
            
            lambda<-param[1]
            bm<-param[2]
            
            GPP<-genomicProfiles(PFM=TFList[[1]],PFMFormat=TFList[[2]],
                                 BPFrequency=DNASequenceSet,PWMThreshold=0.7,lambdaPWM=lambda,boundMolecules=bm,stepSize=100)
            gw<-computeGenomeWideScores(GPP,DNASequenceSet,Access,cores=cores)
            .scores(chipProfile)<-ValidationScore
            .loci(chipProfile)<-Validationloci
            pwm<-computePWMScore(gw,DNASequenceSet,chipProfile,Access,cores=cores)
            occup<-computeOccupancy(pwm)
            chip<-computeChIPProfile(occup,chipProfile,cores=cores)
            gof<-profileAccuracyEstimate(chip,chipProfile,cores=cores)
            optimalList<-list("Occupancy"=occup,"ChIPProfile"=chip,"Gof"=gof)
            
            if(!is.null(filename)){
                save(optimalList,file=paste0(filename,noiseFilter,"_",meth,"OptimalOutputValidation.Rda"))
                save(chipProfile, file=paste0(filename,noiseFilter,"_",meth,"ChIPValidation.Rda"))
            } else {
                save(optimalList,file=paste0(noiseFilter,"optimalOutputValidation.Rda"))
                save(chipProfile, file=paste0(noiseFilter,"ChIPValidation.Rda"))
            }
            
        }
        
        
    }else{
        chipProfile <- processingChIP(profile=ChIP,loci=setSequence, reduce=reduce,
                                      peaks=NULL, chromatinState=NULL,
                                      parameterOptions=PO,cores=cores,normValues=normValues)
        
        
        optimal <- computeOptimal(GPP,DNASequenceSet,chipProfile,Access,PO,cores=cores)
        
        method=method[1]
        ### Computing Optimal Parameters ###
        if(!is.null(filename)){
            save(optimal,file=paste0(filename,noiseFilter,"_",method,"OptimalOutput.Rda"))
            save(chipProfile, file=paste0(filename,noiseFilter,"_",method,"ChIP.Rda"))
        } else {
            save(optimal,file=paste0(noiseFilter,"optimalOutput.Rda"))
            save(chipProfile, file=paste0(noiseFilter,"ChIP.Rda"))
        }
        
        
    }
}



extractParam<-function(data,
                       bm=c(1,10,20,50,100,200,500, 1000, 2000,5000,10000,20000, 50000, 100000,200000,500000, 1000000),
                       lambda=seq(0.25,5, by=0.25),
                       method="MSE"){
    #buildTemplate
    
    temp <- vector("list", ncol(data))
    names(temp)<-colnames(data)
    for(i in seq_along(temp)){
        buffer<-as.vector(data[,i])
        dim(buffer)<-c(length(lambda), length(bm))
        colnames(buffer)<-bm
        rownames(buffer)<-lambda
        #temp[[i]]<-buffer
        if(method %in% c("MSE","ks","geometric")){
            best <- which(buffer==min(buffer), arr.ind=TRUE)
            if(nrow(best)>1)best<-matrix(best[1,],ncol=2)
            temp[[i]]<-c(buffer[best],bm[best[,2]], lambda[best[,1]])
            names(temp[[i]])<-c("Score","bm","lambda")
        } else{
            best <- which(buffer==max(buffer), arr.ind=TRUE)
            if(nrow(best)>1)best<-matrix(best[1,],ncol=2)
            temp[[i]]<-c(buffer[best],bm[best[,2]], lambda[best[,1]])
            names(temp[[i]])<-c("Score","bm","lambda")
        }
    }
    return(temp)
    
}

latexFormat <- function(param, filename){
    ## we will expect to collapse all type of access
    filename<-paste0(filename,".txt")
    files<-gsub("_"," ",names(param[[1]]))
    #files<-gsub("modEncode","",files)
    cat("",file=filename)
    for(i in seq_along(param[[1]])){
        both<-paste0(files[i],"&",param[[1]][[i]]["bm"],"&",param[[1]][[i]]["lambda"],"&",
                     round(param[[1]][[i]]["Score"],digits=3),"&",
                     param[[2]][[i]]["bm"],"&",param[[2]][[i]]["lambda"],"&",
                     round(param[[2]][[i]]["Score"],digits=3),"&",
                     param[[3]][[i]]["bm"],"&",param[[3]][[i]]["lambda"],"&",
                     round(param[[3]][[i]]["Score"],digits=3),"\\","\\","\n","\\hline","\n")
        
        cat(both,file=filename,append=TRUE)
    }
}



.splitDNARanges <- function(DNASequenceSet,cores){
    SplitSeq <- floor(seq(1,length(DNASequenceSet),length.out=cores+1))
    splitDNASequenceSet<-vector("list",(cores))
    start <- SplitSeq[1:(length(SplitSeq)-1)]
    end <- c(SplitSeq[2:(length(SplitSeq)-1)]-1,SplitSeq[length(SplitSeq)])
    
    for(i in seq_len(cores)){
        splitDNASequenceSet[[i]]<-DNASequenceSet[start[i]:end[i]]
    }
    return(splitDNASequenceSet)
}

dataLoadOld<-function(directory,method=c("AUCMean","MSEMean")){
    sets<-vector("list",length(directory))
    for(i in seq_along(directory)){
        sets[[i]]<-vector("list",length(method))
        buffer<-get(load(directory[i]))
        for(j in seq_along(method)){
            sets[[i]][[j]]<-as.vector(buffer[[1]][[2]][[method[j]]])
        }
        names(sets[[i]])<-method
        
    }
    names(sets)<-sapply(strsplit(directory,"sigmoid"),"[[",1)
    return(sets)
}


## function associated with the whole thingy
nameMatch <- function(dhs,null){
    
    dhs<-dhs[,colnames(dhs) %in% colnames(null)]
    null<-null[,colnames(null) %in% colnames(dhs)]
    dhs<-dhs[,!duplicated(colnames(dhs))]
    null<-null[,!duplicated(colnames(null))]
    return(list("dhs"=dhs,"null"=null))
}

bytf<-function(mats,tforder=c("CTCF","BEAF","Hw","Pc","other"),tfdrop=c("Pho2","Pho3","SuH","BCD","abdb","CAD","GT","KR","HB")){
    ## we are going to assume that the input is the output of the previous done
    dhs<-mats$dhs
    null<-mats$null
    
    ## tf dropping
    for(i in tfdrop){
        dhs<-dhs[,-grep(i,colnames(dhs))]
        null<-null[,-grep(i,colnames(null))]
    }
    
    dhssplit<-strsplit(colnames(dhs),"reduce")
    nullsplit<-strsplit(colnames(null),"reduce")
    
    dhsUni<-unique(sapply(dhssplit,"[[",1))
    nullUni<-unique(sapply(nullsplit,"[[",1))
    
    # put in a check
    if(!all.equal(dhsUni,nullUni)){
        stop("Non matching matrices - doublecheck again")
    }
    ordered<-c()
    for(i in tforder){
        if(i != "other"){
            ordered<-c(ordered,grep(i,colnames(dhs)))
            
        } else {
            ordered <-c(ordered,which(
                !grepl(tforder[1],colnames(dhs)) &
                    !grepl(tforder[2],colnames(dhs)) &
                    !grepl(tforder[3],colnames(dhs))&
                    !grepl(tforder[4],colnames(dhs))))
            
        }
    }
    dhs<-dhs[,ordered]
    null<-null[,ordered]
    
    ## new set of namesplit
    dhssplit<-strsplit(colnames(dhs),"reduce")
    
    dhsUni<-unique(sapply(dhssplit,"[[",1))
    
    
    reduce<-c()
    
    for(i in seq_along(dhsUni)){
        
        idx<-which(sapply(strsplit(colnames(dhs),"reduce"),"[[",1)==dhsUni[i])
        
        reduce<-c(reduce,idx[order(as.numeric(sapply(dhssplit[idx],"[[",2)))])
        
    }
    dhs<-dhs[,reduce]
    null<-null[,reduce]
    return(list("dhs"=dhs,"null"=null))
}


metricCompute<-function(mats,region="all"){
    dhs<-mats$dhs
    null<-mats$null
    
    
    if(region=="all"){
        tfs<-unique(sapply(strsplit(colnames(dhs),"reduce"),"[[",1))
        med<-vector("list",length(tfs))
        sds<-vector("list",length(tfs))
        for(i in seq_along(tfs)){
            med[[i]]<-vector("list",2)
            sds[[i]]<-vector("list",2)
            idx<-which(sapply(strsplit(colnames(dhs),"reduce"),"[[",1)==tfs[i])
            med[[i]][[1]]<-median(apply(dhs[,idx],2,max))
            med[[i]][[2]]<-median(apply(null[,idx],2,max))
            sds[[i]][[1]]<-median(apply(dhs[,idx],2,sd))
            sds[[i]][[2]]<-median(apply(null[,idx],2,sd))
            
        }
        names(med)<-gsub("_"," ",tfs)
        names(sds)<-gsub("_"," ",tfs)
    } else{
        buffer<-grep(as.charcter(region),colnames(dhs))
        tfs<-unique(sapply(strsplit(colnames(dhs[buffer]),"reduce"),"[[",1))
        med[[i]]<-vector("list",2)
        sds[[i]]<-vector("list",2)
        for(i in seq_along(tfs)){
            idx<-which(sapply(strsplit(colnames(dhs),"reduce"),"[[",1)==tfs[i])
            med[[i]][[1]]<-max(dhs[,idx])
            med[[i]][[2]]<-max(null[,idx])
            sds[[i]][[1]]<-sd(dhs[,idx])
            sds[[i]][[2]]<-sd(null[,idx])
        }
    }
    return(list("med"=med,"sd"=sds))
}

metricExtract<-function(mats){
    med<-mats$med
    sds<-mats$sd
    ## extracting based on access
    medNULL<-sapply(med,function(x){x[[2]]})
    medDHS<-sapply(med,function(x){x[[1]]})
    sdNULL<-sapply(sds,function(x){x[[2]]})
    sdDHS<-sapply(sds,function(x){x[[1]]})
    return(list("medNULL"=medNULL,"medDHS"=medDHS,"sdNULL"=sdNULL,"sdDHS"=sdDHS))
}

catBuilder<-function(metrics,n=5){
    #creating splits
    catmed<-seq(min(min(metrics$medNULL),min(metrics$medDHS)),
                max(max(metrics$medNULL),max(metrics$medDHS)),length.out=n+1)
    catvar<-seq(min(min(metrics$sdNULL),min(metrics$sdDHS)),
                max(max(metrics$sdNULL),max(metrics$sdDHS)),length.out=n+1)
    # replacemnbt list
    buffer<-metrics
    for(i in seq_len(n)){
        if(i!=n){
            buffer$medNULL[metrics$medNULL >= catmed[i] & metrics$medNULL<catmed[i+1]]<-i
            buffer$medDHS[metrics$medDHS >= catmed[i] & metrics$medDHS<catmed[i+1]]<-i
            buffer$sdNULL[metrics$sdNULL >= catvar[i] & metrics$sdNULL<catvar[i+1]]<-i
            buffer$sdDHS[metrics$sdDHS >= catvar[i] & metrics$sdDHS<catvar[i+1]]<-i
        } else{
            buffer$medNULL[metrics$medNULL > catmed[i] & metrics$medNULL<=catmed[i+1]]<-i
            buffer$medDHS[metrics$medDHS > catmed[i] & metrics$medDHS<=catmed[i+1]]<-i
            buffer$sdNULL[metrics$sdNULL > catvar[i] & metrics$sdNULL<=catvar[i+1]]<-i
            buffer$sdDHS[metrics$sdDHS > catvar[i] & metrics$sdDHS<=catvar[i+1]]<-i
        }
    }
    
    return(buffer)
    
}



################################################
## need to source the script extractParam.R#####
################################################

latexFormat<-function(param,filename){
    
    
    for(i in seq_along(param)){
        both<-paste0(files[i],"&",param[[i]][[1]]$bm,"&",param[[i]][[1]]$lambda,"&",
                     round(param[[i]][[1]]$score,digits=3),"&",param[[i]][[1]]$region,"&",
                     param[[i]][[2]]$bm,"&",param[[i]][[2]]$lambda,"&",
                     round(param[[i]][[2]]$score,digits=3),"&",param[[i]][[2]]$region,
                     "\\","\\","\n","\\hline","\n")
        
        cat(both,file=filename,append=TRUE)
        
    }
}


## quick function
steps<-function(cats,tfs=c("CTCF","Hw","BEAF","other")){
    name<-names(cats$medDHS)
    start<-c()
    end<-c()
    for(i in tfs){
        if(i!="other"){
            start<-c(start,min(grep(i,name)))
            end<-c(end,max(grep(i,name)))
        } else {
            start<-c(start,min(which(!grepl("CTCF",name)
                                     & !grepl("BEAF",name)
                                     & !grepl("Hw",name))))
            end<-c(end,max(which(!grepl("CTCF",name)
                                 & !grepl("BEAF",name)
                                 & !grepl("Hw",name))))
        }
        
    }
    return(list(start,end))
}


regionExtract <- function(mats){
    dhs<-mats$dhs
    null<-mats$null
    
    tags<-unique(sapply(strsplit(colnames(dhs),"reduce"),"[[",1))
    decay<-vector("list", length(tags))
    
    for(i in seq_along(tags)){
        idx<-which(sapply(strsplit(colnames(dhs),"reduce"),"[[",1)==tags[i])
        decay[[i]]$medDHS<-apply(dhs[,idx],2,max)
        decay[[i]]$medNULL<-apply(null[,idx],2,max)
        decay[[i]]$sdDHS<-apply(dhs[,idx],2,sd)
        decay[[i]]$sdNULL<-apply(null[,idx],2,sd)
    }
    names(decay)<-tags
    return(decay)
}

NANFilter<-function(decay){
    nans<-sapply(decay, function(x){
        which(is.na(x))
    })
    if(length(nans)>0){
        nanDrop<-mapply(function(x,y){
            if(length(y)==0)return(x)
            if(length(y)>0)return(x[-y])
        },x=decay,y=nans,SIMPLIFY=F)
        return(nanDrop)
    } else{
        return(decay)
    }
    
}



deltaMetric <- function(decay, tfs=c("CTCF","BEAF","Hw")){
    subs<-vector("list", length(tfs))
    for(i in seq_along(subs)){
        subs[[i]]<-decay[grep(tfs[i],names(decay))]
    }
    names(subs)<-tfs
    
    for(i in seq_along(subs)){
        buffer<-rep(0,7)
        for(j in seq_along(subs[[i]])){
            compVec<-c(20,50,100,150,200,300,500)
            regions<-as.numeric(sapply(strsplit(names(subs[[i]][[j]][[1]]),"reduce"),"[[",2))
            matching<-match(compVec,regions)
            buffer<-cbind(buffer,matching)
            buffer[!is.na(matching),j+1]<-subs[[i]][[j]][[1]]-subs[[i]][[j]][[2]]
        }
        buffer<-buffer[,-1]
        colnames(buffer)<-names(subs[[i]])
        subs[[i]]<-t(buffer)
    }
    return(subs)
    
}

deltaPval<-function(decay){
    res<-matrix(0,ncol=ncol(decay),nrow=ncol(decay))
    for(i in seq_len(ncol(decay))){
        for(j in seq_len(ncol(decay))){
            buffer<-t.test(decay[,i],decay[,j])$p.value
            if(buffer<0.05){
                res[i,j]<-0
            } else{
                res[i,j]<-1
            }
            
            
            
        }
    }
    return(res)
}



overlayMat<-function(mat, method="geometricMean"){
    over<-matrix(0,nrow=14,ncol=17)
    mats<-vector("list",ncol(mat))
    for(i in seq_along(mats)){
        buffer<-mat[,i]
        dim(buffer)<-c(14,17)
        mats[[i]]<-buffer
    }
    
    for(k in seq_along(mats)){
        buffer<-mats[[k]]
        idx<-seq_along(as.vector(mats[[k]]))
        
        
        # define top hits
        top<-head(idx,floor(length(idx)/5))
        
        # Top hits
        if(method %in% c("geometricMean","MSEMean","ksMean","recallMean")){
            ord<-order(buffer,decreasing=F)
        } else{
            ord<-order(buffer,decreasing=T)
        }
        ##creating a ordered value matrix
        optimalMatrix<-matrix(match(idx,ord),ncol=17, nrow=14)
        
        ## Location of top hits in ordered value matrix
        topHits<-which(optimalMatrix<max(top),arr.ind=T)
        
        # Creating empty matrix
        mat<- matrix(0,ncol=17,nrow=14)
        
        ## Keep in mind that there are more Similarity methods
        ## This will pull it towards those optimal Parameters
        ## For better results balence it
        ## or only take two Parameters
        
        ## Overlaying top hits
        mat[topHits]<-1
        
        
        ## Adding top hits to matrix
        over<-over+mat
    }
    over<-over/max(over)
    colnames(over)<-c(1, 10, 20, 50, 100,
                      200, 500,1000,2000, 5000,10000,20000,50000, 100000,
                      200000, 500000, 1000000)
    rownames(over)<-c(0.25, 0.5, 0.75, 1, 1.25,
                      1.5, 1.75, 2, 2.5, 3, 3.5 ,4 ,4.5, 5)
    return(over)
    
}

bytfmat<-function(mat,tf=c("CTCF","BEAF","Hw"),cell=c("BG3","Kc167","S2")){
    res<-vector("list",length(tf)*length(cell))
    count<-1
    tag<-c()
    for(i in seq_along(tf)){
        for(j in seq_along(cell)){
            res[[count]]<-mat[,which(grepl(tf[i],colnames(mat)) & grepl(cell[j],colnames(mat)))]
            tag<-c(tag,paste(tf[i],cell[j]))
            count<-count+1
        }
    }
    names(res)<-tag
    return(res)
}







AccessExtract<-function(subject,query){
    setLocal<-vector("list",length(subject))
    
    for(i in seq_along(subject)){
        localIntersect<-setdiff(subject[i], query)
        setLocal[[i]]<-data.frame("chr"=as.character(seqnames(localIntersect)),"start"=start(localIntersect), "end"=end(localIntersect))
    }
    names(setLocal)<-names(subject)
    return(setLocal)
}



