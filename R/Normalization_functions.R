#' quality control by removing bins of low calling rate and then normalize read counts to remove bias introduced by variabilties in GC content and mappability
#' this function performed quality control where bins with calling rate smaller than the thredhold, GC content out of the range, and mappbility smaller than the threshold were removed.
#' @param RC the read count matrix
#' @param ref_qc the reference genomoe data frame, including seqnames, start, end, width, gc, mapp
#' @param cr.threshold bins with calling rate smaller than this threshold were removed.
#' @param GC.low lower bound of GC content
#' @param GC.up  upper bound of GC content
#' @param mapp.threshold bins with mappability samller than this threshold were removed.
#' @return RC the read counts matrix with high calling rate
#' @return ref_qc the reference genome information of the output RC matrix
#' 

Calling_rate_QC<-function(RC,ref_qc,
                          cr.threshold=0.8,
                          GC.low=0.1,
                          GC.up=0.4,
                          mapp.threshold=0.9){
  
  calling_rate <- apply(RC, 1, function(x){return(sum(x!=0)/length(x))})
  
  kept.bins       <- which(ref_qc$gc>GC.low & ref_qc$gc<GC.up 
                          & ref_qc$mapp>0.9
                          & calling_rate >= cr.threshold)
  
  RC           <- RC[kept.bins,]
  ref_qc       <- ref_qc[kept.bins,]
  
  return(list(RC=RC,ref_qc=ref_qc))
}


#'
#' This function remove bias introduced by GC content
#' @param RC the read count matrix
#' @param GCContent the vector of GC content
#' @param step define the the intervals of GC content where normalization is performed for read counts share the same interval 
#' @return RCNormList the normalized read counts where biases from GC content were removed
#' @export
#' 
CorrectGC<-function(RC,GCContent,step){
  stepseq<-seq(min(GCContent),max(GCContent),by=step)
  #stepseq<-seq(0,100,by=step)
  MasterMedian<-median(RC,na.rm=T)
  MedianGC<-rep(0,length(stepseq)-1)
  RCNormMedian<-RC
  for (i in 1:(length(stepseq)-1)){
    if (i==1){
      ind<-which(GCContent>=stepseq[i] & GCContent<=stepseq[i+1])
    }
    if (i!=1){
      ind<-which(GCContent>stepseq[i] & GCContent<=stepseq[i+1])
    }
    if (length(ind)>0){
      m<-median(RC[ind],na.rm=T)
      if (m>0){
        MedianGC[i]<-m
        RCNormMedian[ind]<-RC[ind]*MasterMedian/m
      }
    }
  }
  RCNormList<-list()
  RCNormList$Median<-MedianGC
  RCNormList$StepGC<-stepseq[1:(length(stepseq)-1)]
  RCNormList$RCNorm<-RCNormMedian
  return(RCNormList)
}


#' This function remove bias introduced by mappability
#' @param RC the read count matrix
#' @param MAPContent the vector of mappability
#' @param step define the the intervals of mappability where normalization is performed for read counts share the same interval 
#' @return RCNormList the normalized read counts where biases from mappability were removed
#' @export
#' 
CorrectMAP<-function(RC,MAPContent,step){
  stepseq<-seq(min(MAPContent),max(MAPContent),by=step)
  #stepseq<-seq(0,100,by=step)
  MasterMedian<-median(RC,na.rm=T)
  MedianMAP<-rep(0,length(stepseq)-1)
  RCNormMedian<-RC
  for (i in 1:(length(stepseq)-1)){
    if (i==1){
      ind<-which(MAPContent>=stepseq[i] & MAPContent<=stepseq[i+1])
    }
    if (i!=1){
      ind<-which(MAPContent>stepseq[i] & MAPContent<=stepseq[i+1])
    }
    
    if (length(ind)>0){
      m<-median(RC[ind],na.rm=T)
      if (m>0){
        MedianMAP[i]<-m
        RCNormMedian[ind]<-RC[ind]*MasterMedian/m
      }
    }
  }
  RCNormList<-list()
  RCNormList$Median<-MedianMAP
  RCNormList$StepMAP<-stepseq[1:(length(stepseq)-1)]
  RCNormList$RCNorm<-RCNormMedian
  return(RCNormList)
}

#' This function remove bias introduced by variability in bin size
#' @param RC the read count matrix
#' @param L the vector of bin size
#' @param step define the the intervals of bin size where normalization is performed for read counts share the same interval 
#' @return RCNormList the normalized read counts where biases from bin size were removed
#' @export
#' 
CorrectSize<-function(RC,L,step){
  stepseq<-seq(min(L),max(L),by=step)
  #stepseq<-seq(0,max(L),by=step)
  MasterMedian<-median(RC,na.rm=T)
  MedianL<-rep(0,length(stepseq)-1)
  RCNormMedian<-RC
  for (i in 1:(length(stepseq)-1)){
    if (i==1){
      ind<-which(L>=stepseq[i] & L<=stepseq[i+1])
    }
    if (i!=1){
      ind<-which(L>stepseq[i] & L<=stepseq[i+1])
    }
    if (length(ind)>0){
      m<-median(RC[ind],na.rm=T)
      if (m>0){
        MedianL[i]<-m
        RCNormMedian[ind]<-RC[ind]*MasterMedian/m
      }
    }
  }
  RCNormList<-list()
  RCNormList$Median<-MedianL
  RCNormList$RCNorm<-RCNormMedian
  return(RCNormList)
}


#' This function remove bias introduced by variabilities in bin size, GC content and mappbility using functions defined above
#' @param RC the read count matrix
#' @param ref_qc the reference genomoe data frame, including seqnames, start, end, width, gc, mapp
#' @param step define the the intervals of mappability where normalization is performed for read counts share the same interval 
#' @return RC_norm the normalized read counts where biases were removed
#' @return log2Rdata the normlaized log2 ratio between RC_norm and cell mean
#' @export
#'

# RC=count.mat
# ref_qc=ref_qc
# cr.threshold=0.8
# GC.low=0.1
# GC.up=0.4
# mapp.threshold=0.9

GC.MAP.normalization<-function(RC, ref_qc){
  
  # start=ref_qc$start
  # end=ref_qc$end
  # exon_size=end-start
  
  start=ref_qc@ranges@start
  end=ref_qc@ranges@start + ref_qc@ranges@width
  exon_size=ref_qc@ranges@width
  
  gc=ref_qc$gc
  map=ref_qc$mapp  
  
  RC_norm=matrix(data=NA,nrow = dim(RC)[1],ncol = dim(RC)[2])
  if (!is.null(rownames(RC))){ rownames(RC_norm)=rownames(RC) }
  if (!is.null(colnames(RC))){ colnames(RC_norm)=colnames(RC) }
  sub=1
  for(sub in 1:dim(RC)[2]){
    ### bin size normalization only when bin size are different###
    if (length(unique(exon_size))>1){
      RCTL <- RC/exon_size
      step <- 5
      RCLNormListIn <- CorrectSize(RCTL[,sub],L=exon_size,step)
      RCLNormIn <- RCLNormListIn$RCNorm
      const <- 1/mean(exon_size)
    }else if(length(unique(exon_size))==1){
      RCLNormIn <- RC[,sub]
      const <- 1
    }
    
    
    ### Mappability normalization ###
    step <- 0.01
    RCMAPNormListIn <- CorrectMAP(RCLNormIn,MAPContent=map,step)
    RCMAPNormIn <- RCMAPNormListIn$RCNorm
    
    ### GC-content Normalization ###
    step <- 0.05
    RCGCNormListIn <- CorrectGC(RCMAPNormIn,GCContent=gc,step)
    RCGCNormIn <- RCGCNormListIn$RCNorm
    
    RC_norm[,sub]=RCGCNormIn
  }
  
  lrr=matrix(data=NA,nrow = dim(RC_norm)[1],ncol = dim(RC_norm)[2])
  rownames(lrr)=rownames(RC_norm)
  colnames(lrr)=colnames(RC_norm)
  for(i in 1:dim(RC_norm)[2]){
    mean=median(RC_norm[,i])
    # mean=mean(RC_norm[,i])
    rc=RC_norm[,i]
    lrr[,i]=log2((rc+const)/mean)
  }
  log2Rdata=lrr
  # log2Rdata<-cbind(pos,log2Rdata)
  # RC_norm  <-cbind(pos,RC_norm)
  return(list(RC_norm=RC_norm,log2Rdata=log2Rdata))
  
}

# norl.res2<-GC.MAP.normalization(RC=count.mat2,ref_qc=ref_qc)
