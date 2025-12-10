# library("modSaRa2")
# library("DNAcopy")
# library("mclust")
# library("FLCNA")

#'
#' CNV detection: segmentation was performed by using CBS method, the priors for 3 CNV states: deletion, normal and duplication were generated using mclust package, and the posterior were updated using GMM. CNV states for segments were then determined
#' this function perform CNV calling
#' @param data.mat the matrix of normalized log2 Ratio
#' @param sample.name the sample/cell names of the samples/cells
#' @param alpha threshold for determining the number of segments in CBS
#' @param chr chromosome information for all bins
#' @param maploc  the location information for all bins
#' @return CNV.res the CNV profiling data frame
#'
#'@export
#'

CBS.function<-function(data.mat,sample.name=NULL,alpha=0.01,chr,maploc){
  if (dim(data.mat)[1]>= dim(data.mat)[2]){
    warning("rows denote cells and columns denote bins, transposing the matrix")
    data.mat <- t(data.mat)
  }
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


#'
#' calculate the prior information for 3 CNV states
#' @param mean.vector the segmental means from all cells/samples
#' @param G number of CNV states
#' @return priors the proportion, mean, and standard deviation of the means for each state
#'
#'@export
#'
mean.init_GapCNV<-function (mean.vector, G=3){

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



#'
#' update the distribution of each CNV state using Gaussian mixture model based on
#' @param CBS.seg.res the data frame containing segments information, including segmental means, segmental sum of square, number of markders
#' @param priors prior information for the CNV states
#' @param L upper bound of iteration time
#' @return CBS.seg.res the data frame contains the states for each segment
#'
#'@export
#'

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
# p.in=para.new$p.new
# mu.in=para.new$mu.new
# sigma.in=para.new$sigma.new
# means
# sum.x.sq
# N
# len
# st
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



