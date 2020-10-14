################################################################################
####################### Lets do some fucking plotting ##########################
################################################################################

#### just in case let's load shit
direc <- getwd()

library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome)
library(RcppRoll)
library(parallel)
library(GenomicRanges)
library(ROCR)
library(rtracklayer)
library(MotifDb)


## sourcing scripts for analysis
setwd("/home/pm16057/ChIPanalyser/ChIPanalyserFinal/ChIPdev")
files <- dir()
for (i in files) source(i)
setwd(direc)

## load data

catshit <- dir()[grepl("Validation", dir()) & grepl("Catchitt", dir()) & grepl("train10", dir())]
centi <- dir()[grepl("Validation", dir()) & grepl("msCENTIPEDE", dir()) & grepl("train10", dir())]
PIQ <- dir()[grepl("Validation", dir()) & grepl("PIQ", dir()) & grepl("train10", dir())]
ChIPanal <- paste0("ChIPanal/",dir("ChIPanal/")[grepl("Validation", dir("ChIPanal/"))  & grepl("train10", dir("ChIPanal/")) & !grepl("ChIPProfile", dir("ChIPanal/"))])


reord <- c(1,5,2,4,6,3,7)

catshit <- catshit[reord]
catshit <- lapply(catshit,function(x)return(x))

centi <- centi[reord]
centi <- lapply(centi,function(x)return(x))

PIQ <- PIQ[reord]
PIQ <- lapply(PIQ,function(x)return(x))

ChIPanal <- ChIPanal[reord]
ChIPanal <- lapply(ChIPanal,function(x)return(x))


for(i in seq_along(reord)){
	catshit[[i]] <- get(load(catshit[[i]]))
	centi[[i]] <- get(load(centi[[i]]))

	PIQ[[i]] <- get(load(PIQ[[i]]))
	ChIPanal[[i]] <- get(load(ChIPanal[[i]]))

}


#catshit<-get(load("/home/pm16057/methodComp/Catchitt_scores_Validation_chr11.Rda"))
#centi<-get(load("/home/pm16057/methodComp/msCENTIPEDE_scores_Validation_chr11.Rda"))
#PIQ<-get(load("/home/pm16057/methodComp/PIQ_scores_Validation_chr11.Rda"))
#ChIPanal <-get(load("/home/pm16057/methodComp/ChIPanal/Validation_chr11.Rda"))





data <-list("Catchitt"=catshit,"msCENTIPEDE"=centi,"PIQ"=PIQ,"ChIPanalyser"=ChIPanal)

valReg <- c("20","50","100","200","500","1000","6755")
methods <- c("AUC","recall","MSE")
tools <- c("Catchitt","msCENTIPEDE","PIQ","ChIPanalyser")

meth<-vector("list",length(methods))
names(meth)<-methods

meth <- lapply(meth, function(x){
			   ret <- vector("list", length(data))
			   names(ret) <- names(data)
			   return(ret)
})

meth <- lapply(meth, function(x){
			   ret <- lapply(x, function(y){
			   				int <- vector("list", length(reord))
			   				names(int) <- valReg
			   				return(int)
			   })
			   return(ret)
})



   
      
for(i in seq_along(data)){
	 for(j in seq_along(methods)){
     		for(k in seq_along(reord)){
			
             if(names(data)[i]!="ChIPanalyser"){
                 if(names(meth)[j]=="AUC"){
                    meth[[j]][[i]][[k]]<-data[[i]][[k]]$auc

                 } else if(names(meth)[j]=="recall"){
                     meth[[j]][[i]][[k]]<-data[[i]][[k]]$rec

                 } else if(names(meth)[j]=="MSE"){
                   meth[[j]][[i]][[k]]<-data[[i]][[k]]$mse

                 }
             }else{

                 loc<-profiles(data[[i]][[k]]$gof)[[1]]
                 sco <- sapply(loc,"[[",methods[j])
                 meth[[j]][[i]][[k]]<-sco
             }

			
		}
	}
}



cats <- get(load( "Catchitt_scores_Validation_chr11Full.Rda"))
ms <- get(load("msCENTIPEDE_scores_Validation_chr11Full.Rda"))
piq <- get(load("PIQ_scores_Validation_chr11Full.Rda"))
chips <- get(load("ChIPanal/Validation_chr11Full.Rda"))

data <-list("Catchitt"=cats,"msCENTIPEDE"=ms,"PIQ"=piq,"ChIPanalyser"=chips)

valReg <- c("20","50","100","200","500","1000","6755")
methods <- c("AUC","recall","MSE")
tools <- c("Catchitt","msCENTIPEDE","PIQ","ChIPanalyser")

methbox<-vector("list",length(methods))
names(methbox)<-methods

methbox <- lapply(methbox, function(x){
			   ret <- vector("list", length(data))
			   names(ret) <- names(data)
			   return(ret)
})




   
      
for(i in seq_along(data)){
	 for(j in seq_along(methods)){
     					
             if(names(data)[i]!="ChIPanalyser"){
                 if(names(methbox)[j]=="AUC"){
                    methbox[[j]][[i]]<-data[[i]]$auc

                 } else if(names(methbox)[j]=="recall"){
                     methbox[[j]][[i]]<-data[[i]]$rec

                 } else if(names(methbox)[j]=="MSE"){
                   methbox[[j]][[i]]<-data[[i]]$mse

                 }
             }else{

                 loc<-profiles(data[[i]]$gof)[[1]]
                 sco <- sapply(loc,"[[",methods[j])
                 methbox[[j]][[i]]<-sco
             }

			
		
	}
}

lboxlim <- list(c(0.48,0.52),c(0.48,0.52),c(0,0.006))
for(k in seq_along(methbox)){
	for(m in seq_along(methbox[[k]])){
		#print(summary(methbox[[k]][[m]]))
		print(lboxlim[[k]])
		methbox[[k]][[m]][methbox[[k]][[m]] <lboxlim[[k]][1] | methbox[[k]][[m]] > lboxlim[[k]][2]] <- NA
	}
}

cols<-c("#4f9da6","#ff5959","#facf5a","#7dd477")
pdf("tool_comp.pdf",width=25,height=13)
#par(mfrow=c(2,2),family="sans",xpd=NA,mar=c(8,4,5,4))
par(mfrow=c(2,3),family="sans",xpd=NA,mar=c(7.5,8.5,6,16.5))
tits<-c("AUC","Recall","MSE")
ylims <-list(c(0.3,1),c(0.3,1),c(0,0.1))
ytext <- c(1.1,1.1,0.115)
mid <- c(0.7,0.7,0.058)
for(i in seq_along(meth)){
	plot(0, type="n", axes = F, xlab ="",ylab="",
		xlim = c (1,7), ylim= ylims[[i]])
	title(xlab = "# of regions in validation", cex.lab = 2, line= 5)
	title(ylab = tits[i], cex.lab=2, line=6.5)
	title(main = tits[i], cex.main =2)
	axis(1, at = seq(1,7),labels = valReg, cex.axis=2)
	axis(2, at =  seq(ylims[[i]][1],ylims[[i]][2], length.out=5),cex.axis=2, las=2)
	text(x = 0, y = ytext[i], labels = LETTERS[i], cex =4)
	for(j in seq_along(meth[[i]])){
		yp <- sapply(meth[[i]][[j]], function(x){
					mean(x, na.rm=T)
		})
		xp <- seq_along(reord)
		points(xp,yp,type="b",col = cols[j], lwd =2.5)
	}
	legend(x= 7,y = mid[i],legend = names(meth[[i]]),col = cols, bty ="n", fill = cols , cex =2)
}

lims <-list(c(0.45,0.55),c(0.45,0.55),c(0,0.006))
labs<-LETTERS[4:6]
yl <- c(0.565,0.565,0.007)
par(mar =c(13.5,8.5,5,8))
for(i in seq_along(methbox)){
  boxplot(methbox[[i]],main=tits[i],col=cols,frame=F,cex.axis=2,cex.main=2,ylim=lims[[i]],las=2)
  #text(x=0,y=max(lims[[i]])+0.15*max(lims[[i]]),labels=labs[i],cex=4)
  text(x=0,y=yl[i],labels=labs[i],cex=4)

}

dev.off()



#lims <-list(c(0,1),c(0,1),c(0,0.17))
#labs<-LETTERS[1:3]
#for(i in seq_along(meth)){
 # boxplot(meth[[i]],main=tits[i],col=cols,frame=F,cex.axis=1,cex.main=1.2,ylim=lims[[i]],las=2)
  #text(x=0,y=max(lims[[i]])+0.15*max(lims[[i]]),labels=labs[i],cex=3)
#}
#bar <-c(155,278,23,80)
#names(bar)<-names(data)
#barplot(bar,col=cols,cex.axis=1,las=2,cex.main=1.2,main="Total Run Time in minutes",cex.main=1)
#text(x=0,y=max(bar)+0.13*max(bar),"D",cex=3)

