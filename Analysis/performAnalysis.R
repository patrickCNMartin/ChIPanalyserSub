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
