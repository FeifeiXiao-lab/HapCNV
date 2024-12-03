library("DNAcopy")
source("/FLCNA_compensate_functions.R")
CBS.function<-function(data.mat,sample.name=NULL,threshold=1,alpha=0.01){
  N<-nrow(data.mat)
  cnv.list<-data.frame() #CNVRuler format
  if (is.null(sample.name)){
    sample.name<-paste0("sample",seq(1:N))
  }
  i<-1
  for(i in 1:N){
    print(i)
    tryCatch(
      # if error occurs, reduce tmp.optimal.nCNV by 1
      {
        CNA.object1 <- CNA(as.vector(data.mat[i,]),rep(1,length(as.vector(data.mat[i,]))),seq(1,length(as.vector(data.mat[i,]))),data.type="logratio")
        smoothed.CNA.object1 <- smooth.CNA(CNA.object1)
        segment1 <- segment(smoothed.CNA.object1, verbose=1, alpha = alpha)
        segment1$output$CN<-0
        segment1$output$CN[which(abs(segment1$output$seg.mean)>threshold)]<-1
        
        tmp.dat<-segment1$output[which(segment1$output$CN==1),]
        tmp.dat$CNV<-rep("Loss",nrow(tmp.dat))
        tmp.dat$CNV[which(tmp.dat$seg.mean>0)]<-"Duplication"
        tmp.cnv.list<-data.frame(Chr=rep(1,nrow(tmp.dat)),
                                 Start=tmp.dat$loc.start,
                                 End=tmp.dat$loc.end,
                                 Event=tmp.dat$CNV,
                                 Sample_ID=sample.name[i])
        cnv.list<-rbind(cnv.list,tmp.cnv.list)
      },
      error=function(e){
        # if there is still an error, skip to next iteration
        # next
      }
    )
    
    
  }
  return(cnv.list)
}

# ------------------------------------------------------------- #

library("modSaRa")
modSaRaC<-function(data.mat,sample.name=NULL,alpha=0.05){
  N<-nrow(data.mat)
  cnv.list<-data.frame() #CNVRuler format
  if (is.null(sample.name)){
    sample.name<-paste0("sample",seq(1:N))
  }
  # i<-1
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = ncol(data.mat), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  tmp.sara.res.dat<-data.frame()
  i<-1
  for(i in 1:N){
    setTxtProgressBar(pb, i)
    tryCatch(
      # if error occurs, reduce tmp.optimal.nCNV by 1
      {
        tmp.sara.res<-modifiedSaRa(data.mat[i,],alpha=alpha)
        if (any(is.na(tmp.sara.res$cnv.end))){
          tmp.sara.res$cnv.end[length(tmp.sara.res$cnv.end)]<-length(data.mat[i,])
        }
        
        tmp.sara.res.dat<-rbind(tmp.sara.res.dat,
                                data.frame(Sample_ID=sample.name[i],
                                           Start=tmp.sara.res$cnv.start,
                                           End=tmp.sara.res$cnv.end,
                                           CNV.state=tmp.sara.res$cnv.state))
        
      },
      error=function(e){
        # if there is still an error, skip to next iteration
        # next
      }
    )
    close(pb)
    
  }
  return(tmp.sara.res.dat)
}

# -------------------------------------------------------------------- #
# -------------------------------------------------------------------- #
# -------------------------------------------------------------------- #
# -------------------------------------------------------------------- #
# HapCNV
library("FLCNA")
library("Seurat")
# also use the upadated function in FLCNA
library("GenomicRanges")
setwd("/20240628_downsampling")
norl.res<-get(load(file="normalized.2.cell.Rdata"))


HapCNV<-function(count.mat,
                 scenario.name,
                 # FLCNA funciton params
                 nclust,
                 lambda  =10,
                 N       =50,
                 iter.LQA=20,
                 # CNA.out params
                 cutoff  =0.35,
                 L       =200
){


log2R<-apply(count.mat,2,function(x){log2((x+0.01)/median(x))})
log2R<-scale(log2R,center = T,scale = FALSE)

# rows represents gene, cols represent cells
output <- FLCNA(K = c(nclust), lambda = lambda, Y=log2R,N=N,iter.LQA = iter.LQA)
setwd("/processed_data")
saveRDS(output,paste0(scenario.name,".FLCNA.clust.RDS"))
# output<-readRDS(paste0("FLCNA.3clust.RDS"))

start=seq(1:nrow(log2R))
end  =seq(1:nrow(log2R))
chr=rep("chr1",nrow(log2R))
tmp.ref.dat<-data.frame(chr=chr,start=start,end=end)
ref<-makeGRangesFromDataFrame(tmp.ref.dat)

# CNA clustering
CNVdata <- CNA.out(mean.matrix = output$mu.hat.best, ref=ref, cutoff=cutoff, L=L,max.cp=50000)
# saveRDS(CNVdata,paste0(scenario.name,"FLCNA.mean3.RDS"))
# CNVdata<-readRDS(paste0("FLCNA.mean.RDS"))
# -------------------------------------------------------------------------- #
# normalization using the FLCNA output
# cluster level CNV indicator matrix
# 1 for no CNV, 2 for there is a CNV
CNV.clust.mat<-matrix(1,nrow = nrow(log2R),ncol = nclust)
# i<-1
for (i in 1:length(CNVdata)){
  tmp.CNVdata<-CNVdata[[i]]
  for (j in 1:nrow(tmp.CNVdata)) {
    CNV.clust.mat[tmp.CNVdata$start[j]:tmp.CNVdata$end[j],i]<-2
  }
}
# CVN indicators for each sample, according to the FLCNA output
# the columns in "alpha.hat.best" corresponds to the rows in "mu.hat.best"
# "s.hat.best" is the results in "alpha.hat.best"
# the order in CNVdata correspond to the row order in "mu.hat.best"
# thus, the columns in "alpha.hat.best" = the rows in "mu.hat.best" = the order in CNVdata
CNV.mat<-CNV.clust.mat[,output[["s.hat.best"]]]
dim(CNV.mat)
# --------------------------------------------------------------------- #
# use median to build ref matrix, 
# each cell c and each bin j, pick bin j showing normal state in cells without c as reference bins
ref_seq.mat.median<-matrix(1,nrow = nrow(log2R),ncol = ncol(log2R))
# dim(ref_seq.mat.median)
colnames(ref_seq.mat.median)<-colnames(log2R)
rownames(ref_seq.mat.median)<-rownames(log2R)

# normalization
count.mat.lib<-Seurat::NormalizeData(count.mat,normalization.method = "RC",scale.factor = 1000000)
count.mat.lib<-as.matrix(count.mat.lib)
# Initializes the progress bar
print("Construct reference for each cell...")
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = ncol(ref_seq.mat.median), # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar
# for each cell, construct reference sequence
for (i in 1:ncol(ref_seq.mat.median)) {
  setTxtProgressBar(pb, i)
  tmp.CNV.mat<-CNV.mat[,-i]
  tmp.seq<-apply(tmp.CNV.mat, 1, function(x){any(x==1)})
  # seq1:there is at least one cell has normal CNV state in these sequence
  seq1   <-which(tmp.seq)
  # seq2:no cell has normal CNV state in these sequence
  seq2   <-which(!tmp.seq)
  
  median.calcul<-function(tmp.RC.mat,
                          tmp.CNV.mat2){
    median.list<-c()
    for (k in 1:nrow(tmp.RC.mat)) {
      median.list<-c(median.list,median(tmp.RC.mat[k,which(tmp.CNV.mat2[k,]==1)]))
    }
    return(median.list)
  }
  # one/not all of cells in this bin has/ have CNV, the referene bin is the median of normal bins
  ref_seq.mat.median[seq1,i]<-median.calcul(tmp.RC.mat=count.mat.lib[seq1,-i],
                                            tmp.CNV.mat2=tmp.CNV.mat[seq1,])
  # all of cells in this bin have CNV, the reference is the median of normal bins in this cell
  seq3  <- which(CNV.mat[,i]==1) # normal sequence in cell i
  ref_seq.mat.median[seq2,i]<-median(count.mat.lib[seq3,i])
  
}
close(pb)

# --------------------------------------------------------------------- #
# normalization
count.mat.lib<-Seurat::NormalizeData(count.mat.lib,normalization.method = "RC",scale.factor = 1000000)
ref.mat.lib<-Seurat::NormalizeData(ref_seq.mat.median,normalization.method = "RC",scale.factor = 1000000)
# log2Rdata<-log2((0.01+count.mat.lib)/(0.001+ref.mat.lib))
log2Rdata<-log2((1+count.mat.lib)/(1+ref.mat.lib))

return(list(ref.mat=ref.mat.lib,
            log2R.norm.mat=log2Rdata))
}

HapCNVres<-HapCNV(count.mat=as.matrix(norl.res$RC_norm[,-c(1,2,3)]),
                  # dim(count.mat)
                  colnames(count.mat),
                  scenario.name="RealData",
                  nclust=c(2,3),
                  # nclust=c(3),
                  lambda  =c(2,5,10),
                  N       =50,
                  iter.LQA=20,
                  # CNA.out params
                  cutoff  =0.35,
                  L       =200)

log2Rdata<-HapCNVres$log2R.norm.mat
log2Rdata<-scale(log2Rdata,center = TRUE,scale = FALSE)
log2Rdata<-t(log2Rdata)
# =============================================================================================== #
# smooth and segmentation
# rows represent sample
modSaRa.smooth.log2Rdata<-modSaRa::smooth(log2Rdata, R = 5, t = 2)
rownames(modSaRa.smooth.log2Rdata)<-rownames(log2Rdata)
# CBS
CBS.res    <-CBS.function(modSaRa.smooth.log2Rdata,
                          sample.name=rownames(modSaRa.smooth.log2Rdata),
                          alpha=0.001,
                          chr=norl.res$RC_norm[,1],
                          maploc=norl.res$RC_norm[,2])
saveRDS(CBS.res,"CBS.res.RDS")


# ==============================================================================================#
# CNVRuler was used to construct CNV regions.
# reformat CBS.res into the data that can be used as the input of CNVRuler
res<-readRDS("CBS.res.RDS")

res$Event[which(res$Event=="Duplication")]<-"Gain"
write.table(res,"res.CNVRuler.txt",row.names = FALSE,quote = FALSE,sep = "\t")

clinical.data<-res
clinical.data$Phenotype<-1
clinical.data$Phenotype[which(substr(clinical.data$Sample_ID,1,1) %in% c("u","U"))]<-0
clinical.data<-clinical.data[,which(colnames(clinical.data) %in% c("Sample_ID","Phenotype"))]
write.table(clinical.data,"clinic.CNVRuler.txt",row.names = FALSE,quote = FALSE,sep = "\t")

