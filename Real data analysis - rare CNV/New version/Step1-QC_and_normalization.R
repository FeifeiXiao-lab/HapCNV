# -------------------------------------------------------------------------------------------- #
# below are normalization code for normalization for GC content and mappability



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
  RCNormList
}

################ Correzione dei RC Per il Bin Size ####################
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
  RCNormList
}

################ Correzione dei RC dalla Mappability ####################
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
  RCNormList
}





GC.MAP.normalization<-function(Y, #raw matrix
                               ref_qc, #info including GC, mappability positions, etc.
                               pos
){
  
  sampname <- colnames(Y)
  start=ref_qc$start
  end=ref_qc$end
  exon_size=end-start
  gc=ref_qc$gc
  map=ref_qc$mapp  
  
  RCTL <- Y/exon_size
  # RCTL <- Y
  RC_norm=matrix(data=NA,nrow = dim(RCTL)[1],ncol = dim(RCTL)[2])
  rownames(RC_norm)=rownames(RCTL)
  colnames(RC_norm)=colnames(RCTL)
  # sub<-1
  for(sub in 1:dim(RCTL)[2]){
    ### Exon size normalization only if InTarget###
    step <- 5
    RCLNormListIn <- CorrectSize(RCTL[,sub],L=exon_size,step)
    RCLNormIn <- RCLNormListIn$RCNorm
    
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
    mean=mean(RC_norm[,i])
    rc=RC_norm[,i]
    lrr[,i]=log2((rc+0.01)/mean)
  }
  log2Rdata=lrr
  log2Rdata<-cbind(pos,log2Rdata)
  RC_norm  <-cbind(pos,RC_norm)
  return(list(RC_norm=RC_norm,log2Rdata=log2Rdata))
  
}

# ----------------------------------------------------------------------------------------- #
setwd("/Data/downsampling/NewDSdata")
# low quality samples have already been removed.
files<-list.files()
files<-files[-c(1:3)]
cell.names<-substr(files,1,nchar(files)-48)

down.sampling.mat<-read.delim(files[1],header = FALSE)
for (i in 2:length(cell.names)){
  print(i)
  tmp.data<-read.delim(files[i],header = FALSE)
  down.sampling.mat<-cbind(down.sampling.mat,tmp.data[,4])
}
colnames(down.sampling.mat)<-c(c("chr","start","end"),cell.names)
down.sampling.mat<-down.sampling.mat[,c(c(1,2,3),which(substr(colnames(down.sampling.mat),1,1) %in% c("t","u")))]

# sample cells annotation

setwd("/Data/downsampling")

sampleinfo<-read.csv("2023 SC Samples Table.csv")
sampleinfo<-sampleinfo[1:93,]
sampleinfo$sample_ID<-paste0(sampleinfo$sample.name,"_")
sampleinfo$sample_ID<-substr(sampleinfo$sample_ID,1,5)

match1<-na.omit(match(substr(colnames(down.sampling.mat),1,5),sampleinfo$sample_ID))
sampleinfo<-sampleinfo[match1,]
na.omit(match(substr(colnames(down.sampling.mat),1,5),sampleinfo$sample_ID))

keep.2cell<-c(c(1,2,3),3+which(sampleinfo$sample.type=="2 cell"))
down.sampling.mat_2cell<-down.sampling.mat[,as.numeric(keep.2cell)]
keep.10cell<-c(c(1,2,3),3+which(sampleinfo$sample.type=="10 cell"))
down.sampling.mat_10cell<-down.sampling.mat[,keep.10cell]


calling.rate<-rowSums(count.mat2>0)/ncol(count.mat2)

# ------------------------------------------------------------------------- #
# QC, remove bins of poor quality based on GC content and Mappability
# 2-cell smaples
setwd("Data/gc_mappability_files")
GC.1kb<-read.delim("pf3d7.mapability_gc.bin1kb.bed")
GC.names       <- paste0(substr(GC.1kb$chr, 7, 8), ":" ,GC.1kb$pos1)
kept.seq       <- GC.names[which(GC.1kb$GC>0.10 & GC.1kb$GC<0.40 & GC.1kb$mappability>0.9 & calling.rate < 0.8)]
rownames(down.sampling.mat_2cell) <- paste0(substr(down.sampling.mat_2cell$chr, 7, 8), ":" ,down.sampling.mat_2cell$start)
count.mat2     <- down.sampling.mat_2cell[match(kept.seq,rownames(down.sampling.mat_2cell)),-c(1:3)]

GC.1kb2 <- GC.1kb[match(kept.seq,GC.names),]
GC.names2 <- paste0(substr(GC.1kb2$chr, 7, 8), ":" ,GC.1kb2$pos1)
ref_qc <- data.frame(seqnames=GC.1kb2$chr, start=GC.1kb2$pos1, end=GC.1kb2$pos2,
                     width=GC.1kb2$pos2-GC.1kb2$pos1,
                     gc=GC.1kb2$GC, mapp=GC.1kb2$mappability)

norl.res.2_cell<-GC.MAP.normalization(Y=count.mat2,ref_qc=ref_qc,pos=ref_qc[,c(1:3)])
norl.res.2_cell[["sampleinfo"]]<-sampleinfo[which(sampleinfo$sample.type=="2 cell"),]

# ------------------------------------------------------------------------- #
# QC, remove bins of poor quality based on GC content and Mappability
# 10-cell smaples
setwd("/Data/gc_mappability_files")
GC.1kb<-read.delim("pf3d7.mapability_gc.bin1kb.bed")
dim(GC.1kb)

GC.names       <- paste0(substr(GC.1kb$chr, 7, 8), ":" ,GC.1kb$pos1)
kept.seq       <- GC.names[which(GC.1kb$GC>0.10 & GC.1kb$GC<0.40 & GC.1kb$mappability>0.9)]
rownames(down.sampling.mat_10cell) <- paste0(substr(down.sampling.mat_10cell$chr, 7, 8), ":" ,down.sampling.mat_10cell$start)
count.mat2     <- down.sampling.mat_10cell[match(kept.seq,rownames(down.sampling.mat_10cell)),-c(1:3)]

GC.1kb2 <- GC.1kb[match(kept.seq,GC.names),]
GC.names2 <- paste0(substr(GC.1kb2$chr, 7, 8), ":" ,GC.1kb2$pos1)
ref_qc <- data.frame(seqnames=GC.1kb2$chr, start=GC.1kb2$pos1, end=GC.1kb2$pos2,
                     width=GC.1kb2$pos2-GC.1kb2$pos1,
                     gc=GC.1kb2$GC, mapp=GC.1kb2$mappability)

norl.res.10_cell<-GC.MAP.normalization(Y=count.mat2,ref_qc=ref_qc,pos=ref_qc[,c(1:3)])
norl.res.10_cell[["sampleinfo"]]<-sampleinfo[which(sampleinfo$sample.type=="10 cell"),]


