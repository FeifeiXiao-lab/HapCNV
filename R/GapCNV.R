#' this function construct pseudo-reference for each cell and perform normalization
#' @param count.mat the count matrix after QC and normalization by GC content and mappbility
#' @param log2R the matrix of log2 Ratio after QC and normalization by GC content and mappbility
#' @param ref the reference information for bins in GRange format
#' @param nclust vector of opitonal number of clusters
#' @param lambda the vector of optional penalty terms
#' @param cutoff the cutoff value to determine whether a bin is a change point
#' @return ref.mat the matrix of pseufo-reference
#' @return log2R.norm.mat the normalized log2 Ratio matrix by pseufo-reference
#'
#'@export
#'
#'
GapCNV<-function(count.mat,
                 log2R,
                 ref,
                 nclust  = 2,
                 lambda  = 10,
                 cutoff  = 0.35
){
  log2R<-scale(log2R,center = T,scale = FALSE)

  # dim(log2R)
  # rows represents gene, cols represent cells
  output <- FLCNA_mod(K = c(nclust), lambda = lambda, y=log2R)

  start=seq(1:nrow(log2R))
  end  =seq(1:nrow(log2R))
  chr=ref@seqnames
  tmp.ref.dat<-data.frame(chr=chr,start=start,end=end)
  ref2<-makeGRangesFromDataFrame(tmp.ref.dat)

  # CNA clustering
  CNVdata <- CNA.out(mean.matrix = output$mu.hat.best, ref=ref2, cutoff=cutoff, max.cp=50000)
  # -------------------------------------------------------------------------- #
  # normalization using the FLCNA output
  # cluster level CNV indicator matrix
  # 1 for no CNV, 2 for there is a CNV
  nclust<-output$K.best
  s.hat.best<-output$s.hat.best
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
  CNV.mat<-CNV.clust.mat[,s.hat.best]
  # dim(CNV.mat)
  # dim(count.mat)
  # --------------------------------------------------------------------- #
  # use median to build ref matrix,
  # each cell c and each bin j, pick bin j showing normal state in count.mat without c as reference bins
  ref_seq.mat.median<-matrix(1,nrow = nrow(count.mat),ncol = ncol(count.mat))
  # dim(ref_seq.mat.median)
  if (!is.null(colnames(count.mat))){
    colnames(ref_seq.mat.median)<-colnames(count.mat)
    rownames(ref_seq.mat.median)<-rownames(count.mat)
  }

  # normalization
  count.mat.lib<-Seurat::NormalizeData(count.mat,normalization.method = "RC",scale.factor = 1000000)
  count.mat.lib<-as.matrix(count.mat.lib)
  print("Construct pseudo-reference for each cell...")
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
  count.mat.lib<-Seurat::NormalizeData(count.mat.lib,normalization.method = "RC",scale.factor = 1000000)
  ref.mat.lib<-Seurat::NormalizeData(ref_seq.mat.median,normalization.method = "RC",scale.factor = 1000000)
  # log2Rdata<-log2((0.01+count.mat.lib)/(0.001+ref.mat.lib))
  log2Rdata<-log2((1+count.mat.lib)/(1+ref.mat.lib))
  log2Rdata<-scale(log2Rdata,center = TRUE,scale=FALSE)
  modSaRa.smooth.log2Rdata<-modSaRa::smooth(as.matrix(log2Rdata), R = 5, t = 2)
  if(!is.null(colnames(log2Rdata))){
    colnames(modSaRa.smooth.log2Rdata)<-colnames(log2Rdata)
  }

  return(list(ref.mat=ref.mat.lib,
              log2R.norm.mat=modSaRa.smooth.log2Rdata))

}

