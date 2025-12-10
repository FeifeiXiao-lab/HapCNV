library("DNAcopy")
library("modSaRa")
library("mclust")
library("FLCNA")
# data.mat=modSaRa.smooth.log2Rdata
# sample.name=rownames(modSaRa.smooth.log2Rdata)
# alpha=0.001
# chr=norl.res$RC_norm[,1]
# maploc=norl.res$RC_norm[,2]

CBS.function<-function(data.mat,sample.name=NULL,alpha=0.01,chr,maploc,priors){
  N<-nrow(data.mat)
  cnv.list<-data.frame() #CNVRuler format
  if (is.null(sample.name)){
    sample.name<-paste0("sample",seq(1:N))
  }
  # chr    = rep(1,ncol(data.mat))
  # maploc = seq(1,ncol(data.mat))
  loc.info<-paste0(chr,"_",maploc)
  Chr_index <- as.numeric(as.factor(chr))
  Chr_cp_index <- 1+which(abs(diff(Chr_index))>=1)
  # i<-1
  # change points
  cp.list<-c()
  # vector of intensities
  Y<-c()
  # i=1
  for(i in 1:N){
    print(i)
    tryCatch(
      # if error occurs, reduce tmp.optimal.nCNV by 1
      {
        CNA.object1 <- CNA(as.vector(data.mat[i,]),chrom=chr,maploc=maploc,data.type="logratio")
        smoothed.CNA.object1 <- smooth.CNA(CNA.object1)
        segment1 <- segment(smoothed.CNA.object1, verbose=1, alpha = alpha)
        
        tmp.dat<-segment1$output
        tmp.cnv.list<-data.frame(Chr   = tmp.dat$chrom,
                                 Start = tmp.dat$loc.start,
                                 End   = tmp.dat$loc.end,
                                 seg.mean  = tmp.dat$seg.mean,
                                 Sample_ID = sample.name[i],
                                 index=i)
        
        loc.info<-paste0(chr,"_",maploc)
        cp=c(match(paste0(tmp.cnv.list$Chr,"_",tmp.cnv.list$Start),loc.info),
             match(paste0(tmp.cnv.list$Chr,"_",tmp.cnv.list$end)[nrow(tmp.cnv.list)],loc.info))
        cp2<-sort(unique(c(cp,Chr_cp_index)))
        # priors <- GMM.init(seg_means=tmp.cnv.list$seg.mean,G=3)
        
        CNAcluster.res<-CNAcluster(Y=as.vector(data.mat[i,]), cp=cp2, L=100, st=3, priors)
        
        
        tmp.cnv.list2<-data.frame(Chr=chr[CNAcluster.res$CNA.start],
                                  Start=maploc[CNAcluster.res$CNA.start],
                                  End=maploc[CNAcluster.res$CNA.end],
                                  # seg.mean = tmp.dat$seg.mean[match(CNAcluster.res$CNA.start,cp)],
                                  CNV.state= CNAcluster.res$CNA.state,
                                  Sample_ID=sample.name[i],
                                  index=i)
        cnv.list<-rbind(cnv.list,tmp.cnv.list2)
        
      },
      error=function(e){
        # if there is still an error, skip to next iteration
        # next
      }
    )
    
  }
  
  cnv.list$width=cnv.list$End-cnv.list$Start
  # CNV.res<-cnv.list[,c("Sample_ID", "Chr","Start","End","width","seg.mean","CNV.state")]
  CNV.res<-cnv.list[,c("Sample_ID", "Chr","Start","End","width","CNV.state")]
  
  
  return(CNV.res)
}


# mean.matrix = output$mu.hat.best
# ref=ref
# cutoff=cutoff
# L=L
# max.cp=50000
CNA.out<-function (mean.matrix, ref, cutoff = 0.80, L = 100,max.cp=50000) {
  # specify the data sepcific prior using a clustering method
  priors <- mean.init(mean.matrix=mean.matrix, ref=ref, cutoff = cutoff)
  
  # add chromosome start/end postions as the change points
  Chr_index <- as.numeric(rep(ref@seqnames@values, ref@seqnames@lengths))
  # Chr_index <- as.numeric(gsub("^.{0,3}", "", rep(ref@seqnames@values, ref@seqnames@lengths)))
  Chr_cp_index <- 1+which(abs(diff(Chr_index))>=1)
  
  
  CNAdata <- vector(mode = "list", length = nrow(mean.matrix))
  # g<-1
  for (g in 1:nrow(mean.matrix)) {
    cp.index1 <- 1 + which(abs(diff(mean.matrix[g, ], 1)) >= 
                             cutoff)
    cp.index <- unique(c(cp.index1, Chr_cp_index))
    cp.index <- sort(cp.index)
    if ((length(cp.index) < 1) | (length(cp.index) > max.cp)) {
      next
    }
    x.inv <- try(res <- CNAcluster(Y = mean.matrix[g, ], 
                                   cp = cp.index, L,priors=priors), silent = TRUE)
    if ("try-error" %in% class(x.inv)) 
      next
    if (length(x.inv$CNA.end) == 0) {
      next
    }
    CNAdata[[g]] <- data.frame(state = res$CNA.state, start = ref@ranges@start[res$CNA.start], 
                               end = (ref@ranges@start + ref@ranges@width)[res$CNA.end], 
                               chr = rep(ref@seqnames@values, ref@seqnames@lengths)[res$CNA.start], 
                               width_bins = (res$CNA.end - res$CNA.start + 1))
  }
  return(CNAdata)
}

# ----------------------------------------------------------------------- #
# update the initiation of mu using GMM given number of CNV states 
# Y = mean.matrix[1, ]
# cp = cp.index
# L
# dim(mean.matrix)
CNAcluster<-function (Y, cp, L, priors,st = st) 
{
  
  EM = gausianMixture(x = Y, cp, priors = priors, L, st = st)
  
  newcp = EM$cp.final
  h = EM$index.final
  CNA.state <- getState(EM = EM)
  return(list(newcp = newcp, h = h, CNA.state = CNA.state$CNA.state, 
              CNA.start = CNA.state$CNA.start, CNA.end = CNA.state$CNA.end))
}

# ----------------------------------------------------------------------- #
getState<- function (EM = EM) 
{
  state = EM$state.new
  cp.f = EM$cp.final
  
  start.index = which(state != 2)
  CNA.start = cp.f[start.index]
  CNA.start = CNA.start[!is.na(CNA.start)]
  CNA.end = cp.f[start.index + 1]
  CNA.end = CNA.end[!is.na(CNA.end)]
  CNA.state = state[start.index]
  # CNA.state[which(CNA.state == 0)] = "del"
  CNA.state[which(CNA.state == 1)] = "del"
  CNA.state[which(CNA.state == 3)] = "dup"
  # CNA.state[which(CNA.state == 4)] = "dup"
  return(list(CNA.state = CNA.state, CNA.start = CNA.start, 
              CNA.end = CNA.end))
}


# ----------------------------------------------------------------------- #
# update the initiation of mu using clustering method given number of CNV states 
# mean.matrix = mean.matrix
# ref = ref
# cutoff=0.8
# mean.matrix is from the FLCNVA output

mean.init<-function (mean.matrix, ref, cutoff = 0.8){
  library("mclust")
  Chr_index <- as.numeric(rep(ref@seqnames@values, ref@seqnames@lengths))
  # Chr_index <- as.numeric(gsub("^.{0,3}", "", rep(ref@seqnames@values, ref@seqnames@lengths)))
  Chr_cp_index <- 1+which(abs(diff(Chr_index))>=1)
  
  
  CNAdata <- vector(mode = "list", length = nrow(mean.matrix))
  # g<-1
  seg_means<-c()
  for (g in 1:nrow(mean.matrix)) {
    cp.index1 <- 1 + which(abs(diff(mean.matrix[g, ], 1)) >= 
                             cutoff)
    cp.index <- unique(c(cp.index1, Chr_cp_index))
    cp.index <- sort(cp.index)
    for (s in 1:(length(cp.index)-1)){
      seg_means<-c(seg_means,mean(mean.matrix[g,cp.index[s]:cp.index[s+1]]))
    }
  }
  
  st = 3
  yMclust <- Mclust(seg_means,G=st,verbose = FALSE)
  init.clust<-yMclust$classification
  mus<-c()
  sds<-c()
  ps <-c()
  for (s1 in 1:st){
    mus<-c(mus,mean(seg_means[init.clust==s1]))
    sds<-c(sds,sd(seg_means[init.clust==s1]))
    ps<-c(ps,sum(init.clust==s1)/length(init.clust))
  }
  reorder<-order(mus)
  mu=mus[reorder]
  sd=sds[reorder]
  p=ps[reorder]
  priors <- list(p = p, mu = mu, sigma = sd)
  return(priors=priors)
  # return(list(priors=priors,init.clust=init.clust))
}
# # -------------------------------------------------------------------- #
# # -------------------------------------------------------------------- #
# # -------------------------------------------------------------------- #
# # -------------------------------------------------------------------- #
# GapCNV
source("FLCNA_compensate_functions.R")

library("Seurat")
# also use the upadated function in FLCNA
library("GenomicRanges")



scenario.name="RealData"
nclust<-c(2,3)
lambda  =c(10,20,50)
N       =50
iter.LQA=20
cutoff  =0.8
L       =200



# find cluster information according to previous FLCNA output
setwd("/RealDataAnlys/Data")
norl.res<-get(load(file="norl.res.DC.map.Rdata"))
count.mat<-as.matrix(norl.res$RC_norm[,-c(1,2,3)])

sample.filter<-which(colnames(count.mat) %in% c("t-g2_S70", #10 cell samples#
                                                # "t-e10", not in this data
                                                "t-f4_S96", # low coverage
                                                "u-h3_S76","u-h6_S11","u-h7_S45" # low coverage
))
count.mat<-count.mat[,-sample.filter]
# -------------------------------------------------------------------------- #
# generate cnv satus from output
log2R<-apply(count.mat,2,function(x){log2((x+0.01)/median(x))})
log2R<-scale(log2R,center = T,scale = FALSE)
output <- FLCNA(K = c(nclust), lambda = lambda, Y=log2R,N=N,iter.LQA = iter.LQA)


chrs<-as.numeric(substr(norl.res$RC_norm[,1],7,8))
start=seq(1:ncol(output$mu.hat.best))
end  =seq(1:ncol(output$mu.hat.best))
tmp.ref.dat<-data.frame(chr=chrs,start=start,end=end)
ref<-makeGRangesFromDataFrame(tmp.ref.dat)
# CNVdata <- FLCNA::CNA.out(mean.matrix = output$mu.hat.best, ref=ref, cutoff=cutoff, L=L)
CNVdata <- CNA.out(mean.matrix = output$mu.hat.best, ref=ref, cutoff=cutoff, L=L,max.cp=50000)
priors <- mean.init(mean.matrix=output$mu.hat.best, ref=ref, cutoff = cutoff)

# -------------------------------------------------------------------------- #
# normalization using the FLCNA output
# cluster level CNV indicator matrix
# 1 for no CNV, 2 for there is a CNV
nclust<-output$K.best
s.hat.best<-output$s.hat.best
CNV.clust.mat<-matrix(1,nrow = nrow(count.mat),ncol = nclust)
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
CNV.mat<-CNV.clust.mat[,s.hat.best]

# kept samples
setwd("/RealDataAnlys/Data/downsampling")
norl.res.down<-get(load(file="normalized.2.cell.Rdata"))
keep.samples.ind<-match(colnames(norl.res.down$RC_norm[,-c(1,2,3)]),colnames(count.mat))
count.mat<- count.mat[,keep.samples.ind]
CNV.mat  <- CNV.mat[,keep.samples.ind]

# --------------------------------------------------------------------- #
# use median to build ref matrix,
# each cell c and each bin j, pick bin j showing normal state in count.mat without c as reference bins
ref_seq.mat.median<-matrix(1,nrow = nrow(count.mat),ncol = ncol(count.mat))
colnames(ref_seq.mat.median)<-colnames(count.mat)
rownames(ref_seq.mat.median)<-rownames(count.mat)

# normalization
count.mat.lib<-Seurat::NormalizeData(count.mat,normalization.method = "RC",scale.factor = 1000000)
count.mat.lib<-as.matrix(count.mat.lib)
print("Construct reference for each cell...")
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = ncol(ref_seq.mat.median), # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar
# for each cell, construct reference sequence
i<-1
for (i in 1:ncol(ref_seq.mat.median)) {
  setTxtProgressBar(pb, i)
  tmp.CNV.mat<-CNV.mat[,-i]
  tmp.seq<-apply(tmp.CNV.mat, 1, function(x){any(x==1)})
  # seq1:there is at least one cell has normal CNV state in these sequence
  seq1   <-which(tmp.seq)
  # seq2:no cell has normal CNV state in these sequence, all cluster have same CNV state, either del(0) or dup(1)
  seq2   <- which( apply(tmp.CNV.mat, 1, function(x){all(x==0) | all(x==2)}) )
  
  # seq3: no cell has normal CNV state in these sequence, all clusters have different CNV state, del(0) or dup(1)
  seq3       <- which( apply(tmp.CNV.mat, 1, function(x){ !any(x==1) & (any(x==0) & any(x==2))}) )
  
  # tmp.RC.mat=count.mat.lib[seq1,-i]
  # tmp.CNV.mat2=tmp.CNV.mat[seq1,]
  median.calcul<-function(tmp.RC.mat,
                          tmp.CNV.mat2){
    median.list<-c()
    # k=1
    for (k in 1:nrow(tmp.RC.mat)) {
      median.list<-c(median.list,median(tmp.RC.mat[k,which(tmp.CNV.mat2[k,]==1)]))
    }
    return(median.list)
  }
  # one/not all of cells in this bin has/ have CNV, the referene bin is the median of normal bins
  if (length(seq1)>0){
    ref_seq.mat.median[seq1,i]<-median.calcul(tmp.RC.mat=count.mat.lib[seq1,-i],
                                              tmp.CNV.mat2=tmp.CNV.mat[seq1,])
  }
  
  # seq2:no cell has normal CNV state in these sequence, all cluster have same CNV state, either del(0) or dup(1)
  if (length(seq2)>0){
    seq4  <- which(CNV.mat[,i]==1) # normal sequence in cell i
    ref_seq.mat.median[seq2,i]<-median(count.mat.lib[seq4,i])
  }
  
  # seq3: no cell has normal CNV state in these sequence, all clusters have different CNV state, del(0) or dup(1)
  # the reference sequence is the median of bins with different CNV states
  median.calcul2<-function(cellind,
                           tmp.RC.mat,
                           tmp.CNV.mat2){
    median.list<-c()
    for (k in 1:nrow(tmp.RC.mat)) {
      if (cellind[k]==0){
        median.list<-c(median.list,median(tmp.RC.mat[k,which(tmp.CNV.mat2[k,]==2)]))
      }else  if (cellind[k]==2){
        median.list<-c(median.list,median(tmp.RC.mat[k,which(tmp.CNV.mat2[k,]==0)]))
      }
    }
    return(median.list)
  }
  
  if (length(seq3)>0){ # avoid error when no element is in seq3
    ref_seq.mat.median[seq3,i]<-median.calcul2(cellind=CNV.mat[seq3,i],
                                               tmp.RC.mat=count.mat.lib[seq3,-i],
                                               tmp.CNV.mat2=tmp.CNV.mat[seq3,])
  }
  
  
}
close(pb)
# --------------------------------------------------------------------- #
# normalization
dim(count.mat.lib)
dim(ref_seq.mat.median)
count.mat.lib<-Seurat::NormalizeData(count.mat.lib,normalization.method = "RC",scale.factor = 1000000)
ref.mat.lib<-Seurat::NormalizeData(ref_seq.mat.median,normalization.method = "RC",scale.factor = 1000000)
# log2Rdata<-log2((0.01+count.mat.lib)/(0.001+ref.mat.lib))
log2Rdata<-log2((1+count.mat.lib)/(1+ref.mat.lib))


log2Rdata<-as.matrix(log2Rdata)
dim(log2Rdata)
log2Rdata<-log2Rdata[,!(substr(colnames(log2Rdata),1,1) %in% c("f","F"))]
log2Rdata<-scale(log2Rdata,center = TRUE,scale=FALSE)
log2Rdata<-t(log2Rdata)
priors<-readRDS(paste0("priors.RDS"))
library("modSaRa")
# smooth and segmentation
# rows represent sample
modSaRa.smooth.log2Rdata<-modSaRa::smooth(log2Rdata, R = 5, t = 2)
rownames(modSaRa.smooth.log2Rdata)<-rownames(log2Rdata)
# CBS
CBS.res<-CBS.function(data.mat=modSaRa.smooth.log2Rdata,
                      sample.name=rownames(modSaRa.smooth.log2Rdata),
                      alpha=0.001,
                      chr=norl.res$RC_norm[,1],
                      maploc=norl.res$RC_norm[,2],
                      priors=priors)

# ========================================================================================== #
# evaluation based on known CNVs

find.seq<-function(res,chr,CNVR.start,CNVR.end){
  samples<-names(table(res$Sample_ID))
  # length(samples)
  n.treated<-length(which(substr(samples,1,1) %in% c("t","T")))
  n.untreated<-length(which(substr(samples,1,1) %in% c("u","U")))
  overlap.seq<-which(toupper(res$Chr)==toupper(chr) & (!(res$Start > CNVR.end | res$End < CNVR.start)))
  tmp.res<-res[overlap.seq,]
  return(tmp.res)
}

F1est<-function(data,
                true.name,
                est.name){
  TP<-which(data[,true.name]==1 & data[,est.name]==1)
  FP<-which(data[,true.name]==0 & data[,est.name]==1)
  FN<-which(data[,true.name]==1 & data[,est.name]==0)
  precison.rate<-length(TP)/(length(TP)+length(FP))
  recall.rate  <-length(TP)/(length(TP)+length(FN))
  F1<-2*(precison.rate*recall.rate)/(precison.rate+recall.rate)
  
  return(c(F1=F1,
           precison.rate=precison.rate,
           recall.rate=recall.rate,
           TP=length(TP),
           FP=length(FP),
           FN=length(FN)))
}

setwd("/RealDataAnlys/Data")
norl.res<-get(load(file="norl.res.DC.map.Rdata"))
pos.data<-norl.res$RC_norm[,c(1,2,3)]

samples<-unique(CBS.res$Sample_ID)
TRUE.CNV.info<-data.frame(cell=samples,GCH1=0, Pf11_1=0, Pf332=1, MDR1=1)
TRUE.CNV.info$est.GCH1<-0
TRUE.CNV.info$est.Pf11_1<-0
TRUE.CNV.info$est.Pf332<-0
TRUE.CNV.info$est.MDR1<-0

res<-CBS.res

# Pf332
tmp.res<-find.seq(res,chr="Pf3D7_11_v3",
                  CNVR.start=1953000,
                  CNVR.end=1966000)
TRUE.CNV.info$est.Pf332[match(unique(tmp.res$Sample_ID),TRUE.CNV.info$cell)] <-1


# MDR1
tmp.res<-find.seq(res,chr="Pf3D7_05_v3",
                  CNVR.start=888000,
                  CNVR.end=971000)
TRUE.CNV.info$est.MDR1[match(unique(tmp.res$Sample_ID),TRUE.CNV.info$cell)] <-1

# Pf11_1
tmp.res<-find.seq(res,chr="Pf3D7_10_v3",
                  CNVR.start=1521000,
                  CNVR.end=1549000)
TRUE.CNV.info$est.Pf11_1[match(unique(tmp.res$Sample_ID),TRUE.CNV.info$cell)] <-1

F1est(data=TRUE.CNV.info[which(substr(TRUE.CNV.info$cell,1,1)=="u"),], true.name="Pf11_1", est.name="est.Pf11_1")
F1est(data=TRUE.CNV.info[which(substr(TRUE.CNV.info$cell,1,1)=="u"),], true.name="Pf332", est.name="est.Pf332")
F1est(data=TRUE.CNV.info[which(substr(TRUE.CNV.info$cell,1,1)=="u"),], true.name="MDR1", est.name="est.MDR1")


# ==============================================================================================#
# CNVRuler input
res<-CBS.res
res$CNV.state[which(res$CNV.state=="dup")]<-"Gain"
res$CNV.state[which(res$CNV.state=="del")]<-"Loss"
colnames(res)<-c("Sample_ID","Chr","Start","End","Width","Event")
res<-res[,c(2,3,4,6,1)]

write.table(res,"res.CNVRuler.GMM.txt",row.names = FALSE,quote = FALSE,sep = "\t")

clinical.data<-res
clinical.data$Phenotype<-1
clinical.data$Phenotype[which(substr(clinical.data$Sample_ID,1,1) %in% c("u","U"))]<-0
clinical.data<-clinical.data[,which(colnames(clinical.data) %in% c("Sample_ID","Phenotype"))]
write.table(clinical.data,"clinic.CNVRuler.GMM.txt",row.names = FALSE,quote = FALSE,sep = "\t")


