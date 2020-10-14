###############################################################################
############################# all data parisng ################################
###############################################################################


library(BSgenome.Dmelanogaster.UCSC.dm6)
library(BSgenome)
library(RcppRoll)
library(parallel)
library(GenomicRanges)
library(ROCR)

library(ChIPanalyser)

source("DataHand.R")

AUCdhs<-get(load("AUCDHSMatrix.Rda"))
#Auc500<-get(load("AUCDHSMatrix500bp.Rda"))
AUCnull<-get(load("AUCNULLMatrix.Rda"))
MSEdhs<-get(load("mseDHSMatrix.Rda"))
#MSE500<-get(load("mseDHSMatrix500bp.Rda"))
MSEnull<-get(load("mseNULLMatrix.Rda"))
GEOdhs<-get(load("geoDHSMatrix.Rda"))
GEOnull<-get(load("geoNULLMatrix.Rda"))


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

## Optimal PAramter extraction
AUCOpti<- bytf(nameMatch(AUCdhs,AUCnull))
MSEOpti<-bytf(nameMatch(MSEdhs,MSEnull))


## need to source the script extractParam.R
AUCparam<-extractParam(AUCOpti,method="AUC")[1:33]
MSEparam<-extractParam(MSEOpti,method="MSE")[1:33]

latexFormat(AUCparam,"AUCOptimal")
latexFormat(MSEparam,"MSEOptimal")


## build your plotting values
aucMats<-metricExtract(metricCompute(bytf(nameMatch(AUCdhs,AUCnull))))
mseMats<-metricExtract(metricCompute(bytf(nameMatch(MSEdhs,MSEnull))))
mins<-min(sapply(mseMats,min))
mseMats<-lapply(mseMats,function(x,y){buf<-x+y;return(log(buf))},1)
# set number of cats
n<-5

# build cats
catsAUC<-catBuilder(aucMats,n=n)
catsMSE<-catBuilder(mseMats,n=n)
#extracting actual values
metrics<-aucMats
catmed<-seq(min(min(metrics$medNULL),min(metrics$medDHS)),
max(max(metrics$medNULL),max(metrics$medDHS)),length.out=n+1)

metrics<-mseMats
catvar<-seq(min(min(metrics$sdNULL),min(metrics$sdDHS)),
max(max(metrics$sdNULL),max(metrics$sdDHS)),length.out=n+1)

## new build with auc and MSEdhs

cats<-c(catsAUC[1:2],catsMSE[3:4])

# prilim plot param
ptCol<-seq_len(n)
ptSize<-seq(4,12,length.out=n)

colfunc<-colorRampPalette(c("#facf5a","#ff5959","#233142"))
cols<-colfunc(5)


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

step<-steps(cats)

start<-step[[1]]
end<-step[[2]]

tags<-names(cats[[1]])
tags<-gsub("modEncode ","",tags,ignore.case=T)
tags<-gsub("Su","su",tags)


## plotting
pdf("AccessMapsAllRegionsMeanAUCreduce.pdf",width=18,height=13.5)
layout(matrix(c(rep(1,6),rep(2,6),5,rep(3,5),6,4,4,4,4,4),byrow=F,ncol=4),
width=c(1,1,1,1))
par(family="mono")
par(xpd=NA)
j=1
par(mar=c(22,18.5,8,4))
plot(0,type="n",ylim=c(start[1],end[1]+1),xlim=c(0.4,2),axes=F,xlab="",ylab="")
axis(2,at=seq_along(cats[[1]][start[j]:end[j]]),labels=tags[start[j]:end[j]],las=2,cex.axis=1.6)
axis(1,at=c(1,2),labels=F)
text(x=c(1,2),y =(0), srt = 45, adj = 1,labels = c("No Accessibility","DHS Accessibility"),
xpd = TRUE,cex=2.8)
text(x=0,y=length(cats[[1]][start[j]:end[j]])+1,labels="A",cex=4)

buf1<-cats[[1]][start[j]:end[j]]
buf2<-cats[[2]][start[j]:end[j]]
buf3<-cats[[3]][start[j]:end[j]]
buf4<-cats[[4]][start[j]:end[j]]
for(i in seq_along(cats[[1]][start[j]:end[j]])){

    points(x=1,y=i,cex=ptSize[buf3[i]],col=cols[buf1[i]],pch=19)
    points(x=2,y=i,cex=ptSize[buf4[i]],col=cols[buf2[i]],pch=19)
}

j=2
par(mar=c(22,18.5,8,4))
plot(0,type="n",ylim=c(start[1],end[1]+1),xlim=c(0.4,2),axes=F,xlab="",ylab="")
axis(2,at=seq_along(cats[[1]][start[j]:end[j]]),labels=tags[start[j]:end[j]],las=2,cex.axis=1.6)
axis(1,at=c(1,2),labels=F)
text(x=c(1,2),y =(0), srt = 45, adj = 1,labels = c("No Accessibility","DHS Accessibility"),
xpd = TRUE,cex=2.8)
text(x=0,y=length(cats[[1]][start[j]:end[j]])+1,labels="B",cex=4)
buf1<-cats[[1]][start[j]:end[j]]
buf2<-cats[[2]][start[j]:end[j]]
buf3<-cats[[3]][start[j]:end[j]]
buf4<-cats[[4]][start[j]:end[j]]
for(i in seq_along(cats[[1]][start[j]:end[j]])){

    points(x=1,y=i,cex=ptSize[buf3[i]],col=cols[buf1[i]],pch=19)
    points(x=2,y=i,cex=ptSize[buf4[i]],col=cols[buf2[i]],pch=19)
}

j=3
par(mar=c(22,18.5,0,4))
plot(0,type="n",ylim=c(start[1],end[1]+1),xlim=c(0.4,2),axes=F,xlab="",ylab="")
axis(2,at=seq_along(cats[[1]][start[j]:end[j]]),labels=tags[start[j]:end[j]],las=2,cex.axis=1.6)
axis(1,at=c(1,2),labels=F)
text(x=c(1,2),y =(0), srt = 45, adj = 1,labels = c("No Accessibility","DHS Accessibility"),
xpd = TRUE,cex=2.8)
text(x=0,y=length(cats[[1]][start[j]:end[j]])+1,labels="C",cex=4)
buf1<-cats[[1]][start[j]:end[j]]
buf2<-cats[[2]][start[j]:end[j]]
buf3<-cats[[3]][start[j]:end[j]]
buf4<-cats[[4]][start[j]:end[j]]
for(i in seq_along(cats[[1]][start[j]:end[j]])){

    points(x=1,y=i,cex=ptSize[buf3[i]],col=cols[buf1[i]],pch=19)
    points(x=2,y=i,cex=ptSize[buf4[i]],col=cols[buf2[i]],pch=19)
}
#j=4
#par(mar=c(22,18.5,0,4))
#plot(0,type="n",ylim=c(start[1],end[1]+1),xlim=c(0.4,2),axes=F,xlab="",ylab="")
#axis(2,at=seq_along(cats[[1]][start[j]:end[j]]),labels=tags[start[j]:end[j]],las=2,cex.axis=1.6)
#axis(1,at=c(1,2),labels=F)
#text(x=c(1,2),y =(0), srt = 45, adj = 1,labels = c("No Accessibility","DHS Accessibility"),
#xpd = TRUE,cex=2.8)
#text(x=0,y=length(cats[[1]][start[j]:end[j]])+1,labels="D",cex=4)
#buf1<-cats[[1]][start[j]:end[j]]
#buf2<-cats[[2]][start[j]:end[j]]
#buf3<-cats[[3]][start[j]:end[j]]
#buf4<-cats[[4]][start[j]:end[j]]
#for(i in seq_along(cats[[1]][start[j]:end[j]])){

#    points(x=1,y=i,cex=ptSize[buf3[i]],col=cols[buf1[i]],pch=19)
#    points(x=2,y=i,cex=ptSize[buf4[i]],col=cols[buf2[i]],pch=19)
#}




par(mar=c(57,8.3,10,18.8))
par(xpd=NA)
plot(0,type="n",ylim=c(0,3),xlim=c(0,1),axes=F,xlab="",ylab="")
axis(4,at=seq(-6,-1),labels=round(catmed,digits=2),las=2,cex.axis=2)
rect(0,-6,1,-5,col=cols[1])
rect(0,-5,1,-4,col=cols[2])
rect(0,-4,1,-3,col=cols[3])
rect(0,-3,1,-2,col=cols[4])
rect(0,-2,1,-1,col=cols[5])
text(x=0.5,y=0,"Median AUC",cex=2)

par(mar=c(4,4,4,4))
plot(0,type="n",axes=F,xlab="",ylab="")

par(mar=c(9.5,6,2,13))
par(xpd=NA)
plot(0,type="n",ylim=c(0,1),xlim=c(0,1.2),axes=F,xlab="",ylab="")
axis(4,at=seq(-5,-1,by=1),labels=paste(round(catvar[1:5],digits=2),"-",round(catvar[2:6],digits=2)),las=2,cex.axis=2,tick=F,line=(-3))
points(x=0.5,y=-5,cex=ptSize[1])
points(x=0.5,y=-4,cex=ptSize[2])
points(x=0.5,y=-3,cex=ptSize[3])
points(x=0.5,y=-2,cex=ptSize[4])
points(x=0.5,y=-1,cex=ptSize[5])
text(y=0.1,x=0.5,"log(SD MSE)",cex=2)
dev.off()


## decay plots

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

### creatin difference between curves

decayMse<-bytf(nameMatch(MSEdhs,MSEnull))
decayAUC<-bytf(nameMatch(AUCdhs,AUCnull))
decaygeo<-bytf(nameMatch(GEOdhs,GEOnull))

decay<-regionExtract(decayAUC)

for(i in seq_along(decay)){
    decay[[i]]<-NANFilter(decay[[i]])
}

decay<-decay[1:33]



pdf("decayAUC.pdf",width=28, height=34)
cols<-c("#233142","#ff5959")
par(family="mono")
par(xpd=NA)
mat<-matrix(c(1:14,0,15:22,0,0,23:33,0,0,0,0),ncol=5,nrow=8, byrow=T)
mat[which(mat==0,arr.ind=T)]<-34:40
layout(mat)
count<-1
labels<-LETTERS[1:3]
for(i in seq_along(decay)){
   par(mar=c(8.5,8,9,11))
   ## for geo and mse
   #ylims<-c(min(min(log(decay[[i]]$medDHS)),min(log(decay[[i]]$medNULL))),
          #max(max(log(decay[[i]]$medDHS)),max(log(decay[[i]]$medNULL))))
  # for AUC
  #ylims<-c(min(min(decay[[i]]$medDHS),min(decay[[i]]$medNULL)),
         #max(max(decay[[i]]$medDHS),max(decay[[i]]$medNULL)))
         ylims<-c(0.5,1)
   xlims<-c(1, min(length(decay[[i]]$medDHS),length(decay[[i]]$medNULL)))

   plot(0,type="n",axes=F,xlab="",ylab="",xlim=xlims,
   ylim=ylims)
   if(i %in% c(1,15,23)){text(x=-0.6,y=ylims[2]+0.21*ylims[2],labels[count],cex=7);count<-count+1}
   tag1<-sapply(strsplit(names(decay[[i]]$medDHS),"reduce"),"[[",2)
   tag2<-sapply(strsplit(names(decay[[i]]$medNULL),"reduce"),"[[",2)
   tag<-intersect(tag1,tag2)
   main<-sapply(strsplit(names(decay[[i]]$medDHS),"reduce"),"[[",1)[1]
   main<-gsub("_"," ",main)
   main<-gsub("modEncode ","",main)
   main<-gsub("Su","su",main)
   axis(1,at=seq_len(xlims[2]),labels=tag,cex.axis=1.6)
   axis(2,at=round(seq(ylims[1],ylims[2],length.out=5),digits=2),labels=round(seq(ylims[1],ylims[2],length.out=5),digits=2),
   cex.axis=1.6,las=2)
   title(xlab="# of Loci",cex.lab=1.6)
   title(ylab="min log(GEO) score",cex.lab=1.6,line=6)
   title(main=main,cex.main=2)
   ## for geo and mse
   #lines(seq_along(decay[[i]]$medDHS),log(decay[[i]]$medDHS),col=cols[1],type="b",lwd=2)
   #lines(seq_along(decay[[i]]$medNULL),log(decay[[i]]$medNULL),col=cols[2],type="b",lwd=2)
   ##for auc
   lines(seq_along(decay[[i]]$medDHS),decay[[i]]$medDHS,col=cols[1],type="b",lwd=2)
   lines(seq_along(decay[[i]]$medNULL),decay[[i]]$medNULL,col=cols[2],type="b",lwd=2)

    legY<-ylims[1]+((ylims[2]-ylims[1])/2)
    print(legY)
   legend(x=xlims[2],y=legY,legend=c("DHS","No Access"),fill=c(cols[1],cols[2]),bty="n",cex=1.6)
}

par(xpd=NA)
plot(x=c(-50,50),y=c(0.8,0),type="n",xlim=c(0,1),ylim=c(1,2),lty=2,lwd=2,axes=F,xlab="",ylab="")
plot(x=c(-50,50),y=c(0.8,0),type="n",xlim=c(0,1),ylim=c(1,2),lty=2,lwd=2,axes=F,xlab="",ylab="")
plot(x=c(-50,50),y=c(0.8,0),type="n",xlim=c(0,1),ylim=c(1,2),lty=2,lwd=2,axes=F,xlab="",ylab="")
plot(x=c(-50,50),y=c(0.8,0),type="n",xlim=c(0,1),ylim=c(1,2),lty=2,lwd=2,axes=F,xlab="",ylab="")
plot(x=c(-50,50),y=c(0.8,0),xlim=c(0,1),ylim=c(1,2),type="l",lty=2,lwd=2,axes=F,xlab="",ylab="")
plot(x=c(-50,50),y=c(0.8,0),xlim=c(0,1),ylim=c(1,2),type="l",lty=2,lwd=2,axes=F,xlab="",ylab="")
dev.off()




decayMse<-regionExtract(bytf(nameMatch(MSEdhs,MSEnull)))
decayAUC<-regionExtract(bytf(nameMatch(AUCdhs,AUCnull)))
#decaygeo<-regionExtract(bytf(nameMatch(GEOdhs,GEOnull)))
decaylist<-list(decayAUC,decayMse)
decay<-vector("list", 3)

for(j in seq_along(decay)){
   decay[[j]]<-vector("list", length(decaylist[[j]]))

   for(i in seq_along(decay[[j]])){
      decay[[j]][[i]]<-NANFilter(decaylist[[j]][[i]])
   }
   names(decay[[j]])<-names(decaylist[[j]])
   decay[[j]]<-decay[[j]][grepl("282",names(decay[[j]])) | grepl("3665",names(decay[[j]])) | grepl("3716",names(decay[[j]]))]
}
names(decay)<-c("AUC","MSE")

#### decay for top three using delta

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

decayMse<-bytf(nameMatch(MSEdhs,MSEnull))
decayAUC<-bytf(nameMatch(AUCdhs,AUCnull))
decaygeo<-bytf(nameMatch(GEOdhs,GEOnull))

decay<-regionExtract(decayMse)

for(i in seq_along(decay)){
    decay[[i]]<-NANFilter(decay[[i]])
}

decayAUC<-deltaMetric(decay)
decayMse<-deltaMetric(decay)


decayAUCpval<-lapply(decayAUC,deltaPval)




### plotting baby

pdf("boxplotDecayAUCpval.pdf",width=18,height=14)
par(family="mono")
par(xpd=NA)
layout(matrix(1:6,ncol=2,byrow=F),width=c(4,1.6))
par(mar=c(8,8,8,4))
mains<-c("CTCF","BEAF-32","su(Hw)")
cols<-c("#facf5a","#ff5959","#4f9da6")
tagsauc<-LETTERS[1:3]
for(i in seq_along(decayAUC)){
  colnames(decayAUC[[i]])<-c(20,50,100,150,200,300,500)
  boxplot(decayAUC[[i]],las=2,col=cols[i],frame=F,cex.axis=1.4,cex.names=1.4,ylim=c(-0.12,0.4))
  lines(y=c(0,0),x=c(0.3,7.4),lty=2,lwd=1.8)
  title(xlab="Number of Regions Selected", line=6,cex.lab=1.6)
  title(ylab="Delta max AUC Score (DHS - NULL)", line=6,cex.lab=1.6)
  title(main=mains[i],cex.main=1.6)

  text(x=0,y=0.6,tagsauc[i],cex=4)

}
tagsauc<-LETTERS[4:6]
colfunc<-colorRampPalette(c("#4f9da6","#edf5f6"))
cols<-colfunc(30)
par(mar=c(8,2,4,4))
for(i in seq_along(decayAUCpval)){

    image(x=1:7,y=1:7,t(decayAUCpval[[i]]),zlim=c(0,1),col=cols,axes=F,xlab="",ylab="")
    axis(1,at=1:7,labels=c(20,50,100,150,200,300,500),cex.axis=1.6)
    axis(2,at=1:7,labels=c(20,50,100,150,200,300,500),cex.axis=1.6)
    text(x=-0.8,y=8.3,tagsauc[i],cex=4)
    title(main="Differences between # regions - (p<0.05) ",cex.main=1.6)

}

dev.off()



# removing raster images ( if you add them again check layout)
for(i in seq_along(decayAUCpval)){
  legend_image<-as.raster(matrix(rev(cols),ncol=1))
  par(mar=c(7,0.5,3.1,0.5))
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
  text(x=1.6, y =seq(0,1,l=5) , labels = seq(0,1,l=5),cex=1.6)
  rasterImage(legend_image, 0, 0, 1,1)
}





###




#decaysub<-decay[grepl("282",names(decay)) | grepl("3665",names(decay)) | grepl("3716",names(decay))]

pdf("decaytop3.pdf",width=16, height=7.5)
cols<-c("#233142","#ff5959")
par(family="mono")
par(xpd=NA)
par(mfrow=c(2,3))
par(oma=c(0,5,0,0))
lets<-LETTERS[1:9]
count<-1
for(j in seq_along(decay)){

for(i in seq_along(decay[[j]])){

   par(mar=c(5,8,5,11))
   ## for geo and mse
   if(names(decay)[j] %in% c("MSE","GEO")){
   #ylims<-c(min(min(log(decay[[j]][[i]]$medDHS)),min(log(decay[[j]][[i]]$medNULL))),
          #max(max(log(decay[[j]][[i]]$medDHS)),max(log(decay[[j]][[i]]$medNULL))))
   ylims<-c(1,6)

  } else{
  #ylims<-c(min(min(decay[[j]][[i]]$medDHS),min(decay[[j]][[i]]$medNULL)),
         #max(max(decay[[j]][[i]]$medDHS),max(decay[[j]][[i]]$medNULL)))
         ylims<-c(0.5,1)
  }
   xlims<-c(1, min(length(decay[[j]][[i]]$medDHS),length(decay[[j]][[i]]$medNULL)))

   plot(0,type="n",axes=F,xlab="",ylab="",xlim=xlims,
   ylim=ylims)




   tag1<-sapply(strsplit(names(decay[[j]][[i]]$medDHS),"reduce"),"[[",2)
   tag2<-sapply(strsplit(names(decay[[j]][[i]]$medNULL),"reduce"),"[[",2)
   tag<-intersect(tag1,tag2)
   main<-sapply(strsplit(names(decay[[j]][[i]]$medDHS),"reduce"),"[[",1)[1]
   main<-gsub("_"," ",main)
   main<-gsub("modEncode ","",main)
   axis(1,at=seq_len(xlims[2]),labels=tag,cex.axis=1.6)
   axis(2,at=round(seq(ylims[1],ylims[2],length.out=5),digits=2),labels=round(seq(ylims[1],ylims[2],length.out=5),digits=2),
   cex.axis=1.6,las=2)
   title(xlab="# of Loci",cex.lab=1.6)
   text(x=-1,y=ylims[2]+((ylims[2]-ylims[1])/5),labels=lets[count],cex=4)
   title(main=main,cex.main=2)


   if(names(decay)[j] %in% c("MSE","GEO")){
   lines(seq_along(decay[[j]][[i]]$medDHS),log(decay[[j]][[i]]$medDHS),col=cols[1],type="b",lwd=2)
   lines(seq_along(decay[[j]][[i]]$medNULL),log(decay[[j]][[i]]$medNULL),col=cols[2],type="b",lwd=2)
 }else{
   lines(seq_along(decay[[j]][[i]]$medDHS),decay[[j]][[i]]$medDHS,col=cols[1],type="b",lwd=2)
   lines(seq_along(decay[[j]][[i]]$medNULL),decay[[j]][[i]]$medNULL,col=cols[2],type="b",lwd=2)
  }
    legY<-ylims[1]+((ylims[2]-ylims[1])/2)
    print(legY)
   legend(x=xlims[2],y=legY,legend=c("DHS","No Access"),fill=c(cols[1],cols[2]),bty="n",cex=1.6)
   if(names(decay)[j]=="MSE"){
   #title(ylab="min log(MSE) score",cex.lab=1.6,line=6)
   if(i==1){

     text(x=(-2),y=legY,"log(MSE)",cex=2.5,srt=90)
   }
 } else if(names(decay)[j]=="GEO"){
   #title(ylab="min log(GEO) score",cex.lab=1.6,line=6)
   if(i==1){

     text(x=(-2),y=legY,"log(Geo)",cex=2.5,srt=90)
   }
 } else{
    #title(ylab="max AUC score",cex.lab=1.6,line=6)
    if(i==1){

      text(x=(-2),y=legY,"AUC",cex=2.5,srt=90)
    }
 }

count<-count+1
}
}
dev.off()
