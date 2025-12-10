library("modSaRa2")
library("DNAcopy")
library("mclust")
library("FLCNA")
source("/FLCNA_compensate_functions.R")

CBS.function<-function(data.mat,sample.name=NULL,alpha=0.01,chr,maploc){
  N<-nrow(data.mat)
  cnv.list<-data.frame() #CNVRuler format
  if (is.null(sample.name)){
    sample.name<-paste0("sample",seq(1:N))
  }
  pos<-paste0(chr,"_",maploc)
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
        # segment1$output$CN<-0
        # segment1$output$CN[which(abs(segment1$output$seg.mean)>threshold)]<-1
        # which(segment1$output$loc.start>segment1$output$loc.end)
        tmp.dat<-segment1$output
        tmp.cnv.list<-data.frame(Chr=tmp.dat$chrom,
                                 Start=tmp.dat$loc.start,
                                 End=tmp.dat$loc.end,
                                 seg.mean= tmp.dat$seg.mean,
                                 Sample_ID=sample.name[i],
                                 index=i)
        
        # calculate sum of sqaures for each segments
        tmp.pos.start<-paste0(tmp.cnv.list$Chr,"_",tmp.cnv.list$Start)
        tmp.pos.end  <-paste0(tmp.cnv.list$Chr,"_",tmp.cnv.list$End)
        
        sum.x.sq<-rep(0,nrow(tmp.cnv.list))
        num.marker <- rep(0,nrow(tmp.cnv.list))
        for (j in 1:nrow(tmp.cnv.list)) {
          tmp.start<-match(tmp.pos.start[j],pos)
          tmp.end  <-match(tmp.pos.end[j],pos)
          ##Find segment sum.x.sq
          sum.x.sq[j] <- sum(as.vector(data.mat[i,tmp.start:tmp.end])^2)
          num.marker[j]      <- tmp.end-tmp.start
        }
        
        tmp.cnv.list$sum.x.sq=sum.x.sq
        tmp.cnv.list$num.marker = num.marker
        cnv.list<-rbind(cnv.list,tmp.cnv.list)
        
        
      },
      error=function(e){
        # if there is still an error, skip to next iteration
        # next
      }
    )
  }
  
  priors = mean.init_GapCNV(mean.vector=cnv.list$seg.mean, G=3)
  
  GMM.res<-gausianMixture_GapCNV(CBS.seg.res=cnv.list, priors=priors, L=200, criteria=0.001, G=3)
  
  CNV.res<-GMM.res$CNV.res[,c(c("Sample_ID", "Chr","Start","End","width","seg.mean","CNV.state"))]
  # colnames(GMM.res$CNV.res)
  CNV.res<-CNV.res[which(CNV.res$CNV.state!="Norm"),]
  
  return(CNV.res)
}


# ---------------------------------------------------------------------------------------#
# this function is for GapCNV CNV calling, based on the segments derived from CBS
mean.init_GapCNV<-function (mean.vector, G=3){
  library("mclust")
  
  yMclust <- Mclust(mean.vector,G=G,verbose = FALSE)
  init.clust<-yMclust$classification
  mus<-c()
  sds<-c()
  ps <-c()
  for (s1 in 1:G){
    mus<-c(mus,mean(mean.vector[init.clust==s1]))
    sds<-c(sds,sd(mean.vector[init.clust==s1]))
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



gausianMixture_GapCNV <- function(CBS.seg.res, priors, L, criteria=0.01, G) {
  
  means=CBS.seg.res$seg.mean
  sum.x.sq=CBS.seg.res$sum.x.sq
  N=nrow(CBS.seg.res)
  len=CBS.seg.res$num.marker
  st=G
  state    <- vector()
  para.new <- updateEM_GapCNV(p.in=priors$p, mu.in=priors$mu, sigma.in=priors$sigma, means, sum.x.sq, N, len, st)
  tmp.diff <- sum((priors$mu-para.new$mu.new)**2)
  iter = 1
  while(iter<L & tmp.diff>criteria){
    # print(iter)
    old.mu   <- para.new$mu.new
    para.new <- updateEM_GapCNV(p.in=para.new$p.new, mu.in=para.new$mu.new, sigma.in=para.new$sigma.new, means, sum.x.sq, N, len, st)
    tmp.diff <- sum(c(old.mu-para.new$mu.new)**2)
    iter = iter + 1
  }
  
  for (pt in 1:N) {
    state[pt] <- which.max(para.new$p[pt,])
  }
  
  state.new    <- state
  # assign CNV state for each segments
  CNV.state = state.new
  CNV.state[which(state.new == 1)] = "Del"
  CNV.state[which(state.new == 2)] = "Norm"
  CNV.state[which(state.new == 3)] = "Dup"
  CBS.seg.res$CNV.state=CNV.state
  CBS.seg.res$width = CBS.seg.res$End - CBS.seg.res$Start
  
  return(list(p.final = para.new$p.new, mu.final = para.new$mu.new, sigma.final = para.new$sigma.new, CNV.res=CBS.seg.res))
}

#' @title Update parameters using EM algorithm
#' 
#' @description In the Gaussian Mixture Model, parameters will be updated based on EM algorithm.
#'
#' @param p.in Initial probability for each CNA cluster.
#' @param mu.in Initial mean value for each CNA cluster.
#' @param sigma.in Initial variance for each CNA cluster.
#' @param means Mean value vector for each segment.
#' @param sum.x.sq Sum of squared mean values for each segment.
#' @param N Number of candiate CNAs.
#' @param len Width of candiate CNAs.
#' @param st Number of assumed states in the EM algorithm.
#' @return The return is the updated parameters using EM algorithm
#' \item{p.new}{Updated probability for each CNA cluster.}
#' \item{mu.new}{Updated mean value for each CNA cluster.}
#' \item{sigma.new}{Updated variance for each CNA cluster.}
#' 

updateEM_GapCNV <- function (p.in, mu.in, sigma.in, means, sum.x.sq, N, len, st) {
  ##Calcute the prob of each segment belong to each state
  p <- dens <- matrix(NA, N, st)
  for (i in 1:N) {
    a <- rep(NA, st)
    for (j in 1:st) {
      dens[i,j] <- dnorm(means[i], mu.in[j], sqrt(sigma.in[j]/len[i]), log=TRUE)
      a[j]      <- log(p.in[j]) + dens[i,j]
    }
    max.a     <- max(a)
    for (k in 1:st) {
      p[i,k] <- (exp(a[k]-max.a))/sum(exp(a-max.a))
    }
  }
  ##Update p, mu, sigma
  p <- na.omit(p)
  p.new <- mu.new <- sigma.new <- rep(NA, st)
  for (k in 1:st) {
    p.new[k]  <- sum(p[,k])/N
    mu.new[k] <- (sum(means*len*p[,k]))/(sum(len*p[,k]))
    sigma.new[k] <- sum(p[,k]*(sum.x.sq-2*len*mu.new[k]*means+len*mu.new[k]^2))/(sum(p[,k]*len))
  }
  return(list(p.new = p.new, mu.new = mu.new, sigma.new = sigma.new, p = p))
}





# -------------------------------------------------- #
# !!!functions mean.init,getState,CNAcluster,CNA.out are adjusted from FLCNA, to make FLCNA suitable for haploid cells

CNA.out<-function (mean.matrix, ref, cutoff = 0.35, L = 100,max.cp=50000) {
  # specify the data sepcific prior using a clustering method
  priors <- mean.init(mean.matrix=mean.matrix, ref=ref, cutoff = cutoff)
  
  # add chromosome start/end postions as the change points
  Chr_index <- as.numeric(rep(ref@seqnames@values, ref@seqnames@lengths))
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
getState <- function (EM = EM) 
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

mean.init<-function (mean.matrix, ref, cutoff = 0.8){
  library("mclust")
  Chr_index <- as.numeric(rep(ref@seqnames@values, ref@seqnames@lengths))
  Chr_cp_index <- 1+which(abs(diff(Chr_index))>=1)
  
  
  CNAdata <- vector(mode = "list", length = nrow(mean.matrix))
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
}

# # -------------------------------------------------------------------- #
# # -------------------------------------------------------------------- #
# # -------------------------------------------------------------------- #
# # -------------------------------------------------------------------- #
# GapCNV

library("Seurat")
# also use the upadated function in FLCNA
library("GenomicRanges")



scenario.name="RealData"
nclust<-c(2,3)
lambda  =c(2,5,10,20,50)
N       =50
iter.LQA=20
cutoff  =0.35
L       =200



# rows represents gene, cols represent cells
setwd("/RealDataAnlys/Data")
norl.res<-get(load(file="norl.res.DC.map.Rdata"))
count.mat<-as.matrix(norl.res$RC_norm[,-c(1,2,3)])

sample.filter<-which(colnames(count.mat) %in% c("t-g2_S70", #10 cell samples#
                                                # "t-e10", not in this data
                                                "t-f4_S96", # low coverage
                                                "u-h3_S76","u-h6_S11","u-h7_S45" # low coverage
))
count.mat<-count.mat[,-sample.filter]
log2R<-apply(count.mat,2,function(x){log2((x+0.01)/median(x))})
log2R<-scale(log2R,center = T,scale = FALSE)
output <- FLCNA(K = c(nclust), lambda = lambda, Y=log2R,N=N,iter.LQA = iter.LQA)

# -------------------------------------------------------------------------- #
# generate cnv satus from output
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
# dim(ref_seq.mat.median)
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
  
  # which(!is.na(match(seq1,seq2)))
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
log2Rdata<-log2((1+count.mat.lib)/(1+ref.mat.lib))

log2Rdata<-as.matrix(log2Rdata)
log2Rdata<-log2Rdata[,!(substr(colnames(log2Rdata),1,1) %in% c("f","F"))]
log2Rdata<-scale(log2Rdata,center = TRUE,scale=FALSE)
log2Rdata<-t(log2Rdata)

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
CBS.res<-get(load(file="CBS.res.GMM.Rdata"))



# ==============================================================================================#
# CNVRuler was used to construct CNV regions.
# reformat CBS.res into the data that can be used as the input of CNVRuler

res <- CBS.res
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


