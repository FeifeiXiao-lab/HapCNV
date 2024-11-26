# library("FLCNA")
#' @title Shared CNA output
#'
#' @description This function clusters the identified change-points to make final CNA calling. The potential CNA segments between two neighbor candidate change-points are assigned to different copy number states according to the estimated mean matrix from FLCNA R function. We use three clusters including duplication, normal state and deletion. A Gaussisan Mixture Model based clustering strategy was applied to assign each segment to the most likely cluster/state.
#'
#' @param mean.matrix The cluster mean matrix estimated from FLCNA R function.
#' @param cutoff Cutoff value to further control the number of CNAs, the larger value of cutoff, the smaller number of CNAs. The default is 0.8.
#' @param L Repeat times in the EM algorithm, defaults to 100.
#'
#'
#' @return The return is the clustered CNA segments with start position, end position and copy number states.
#' \item{state}{The CNA states assigned.}
#' \item{start}{The start point of CNAs.}
#' \item{end}{The end point of CNAs.}
#' \item{chr}{Chromosome of CNAs.}
#' \item{width}{The width of CNAs.}
#'
#' @export
#'

CNA.out<-function (mean.matrix, ref, cutoff = 0.35, L = 100,max.cp=50000) {
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
    if ("try-error" %in% class(x.inv)){
      next
    }
    if (length(x.inv$CNA.end) == 0) {
      next
    }
    tmp.data <- data.frame(state = res$CNA.state, start = ref@ranges@start[res$CNA.start],
                           end = (ref@ranges@start)[res$CNA.end],
                           # end = (ref@ranges@start + ref@ranges@width-1)[res$CNA.end],
                           chr = rep(ref@seqnames@values, ref@seqnames@lengths)[res$CNA.start],
                           width_bins = (res$CNA.end - res$CNA.start + 1))
    tmp.data <- tmp.data[which(tmp.data$start<tmp.data$end),]
    CNAdata[[g]] <- tmp.data

  }
  return(CNAdata)
}


#' @title CNAcluster
#'
#' @description This function clusters CNAs into different states using a Gaussian Mixed Model based clustering strategy.
#' @param Y The numeric vector of the intensities of markers, which is the estimated mean vector in our study.
#' @param cp The numeric vector of the position index for the identified change-points.
#' @param L Repeat times in the EM algorithm, defaults to 100.
#' @param priors prior information for the CNV states
#'
#'
#' @return The return is the clustered CNA segments with the start position and end position, length of the CNA and the copy number states (duplication or deletion). It also returns a vector of final candidates of change-points.
#' \item{newcp}{The final list of change-points.}
#' \item{h}{The bandwidth used for the identification of change-points.}
#' \item{CNA.state}{Copy number state for each CNA.}
#' \item{CNA.start}{Start position of each CNA.}
#' \item{CNA.end}{End position of each CNA.}
#' @export
#'

CNAcluster<-function (Y, cp, L, priors,st = st)
{
  st=length(priors$mu)
  EM = gausianMixture(x = Y, cp, priors = priors, L, st = st)

  newcp = EM$cp.final
  h = EM$index.final
  CNA.state <- getState(EM = EM)
  return(list(newcp = newcp, h = h, CNA.state = CNA.state$CNA.state,
              CNA.start = CNA.state$CNA.start, CNA.end = CNA.state$CNA.end))
}


#' @title CNA states
#'
#' @description This function uses output of Gaussian Mixture Model to obtain different CNA states.
#' @param EM The output of Gaussian Mixture Model for clustering.
#'
#'
#' @return The return is the estimated CNA information.
#' \item{CNA.state}{Copy number state for each CNA.}
#' \item{CNA.start}{Start position of each CNA.}
#' \item{CNA.end}{End position of each CNA.}
#' @export
#'
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
#' @title CNA states
#'
#' @description This function uses output of Gaussian Mixture Model to obtain different CNA states.
#' @param mean.matrix Tmean.matrix is from the FLCNVA output
#' @param ref the genome reference in GRanges format
#'
#' @return The return is the estimated CNA information.
#' \item{CNA.state}{Copy number state for each CNA.}
#' \item{CNA.start}{Start position of each CNA.}
#' \item{CNA.end}{End position of each CNA.}
#' @export
#'

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




#' @title Esimation of mean matrix
#'
#' @description Esimation of mean matrix based on local quadratic approximation (LQA).
#'
#' @param k k-th cluster index.
#' @param K1 K1 cluster in total.
#' @param index.max Initial cluster index assigned accroding to posterior probability.
#' @param y n by p data matrix.
#' @param mu.t.all K by p mean matrix from previous EM-step.
#' @param mu.no.penal  K by p mean matrix of unpenalized estimates (lambda=0).
#' @param sigma.all p by p diagnal covariance matrix.
#' @param alpha n by K posterior probability matrix.
#' @param lambda Tuning parameter.
#' @param iter.num Max iterations in local quadratic approximation (LQA).
#' @param eps.LQA LQA stop criterion.
#' @param eps.diff Lower bound of mean difference.
#'
#'
#' @return Estimated mu hat for k-th cluster with totally K1 cluster.
#'
# -------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------- #
FLCNA_LQA_mod<-function (k, K1, index.max, y, mu.t.all, mu.no.penal, sigma.all,
                         alpha, lambda, iter.num = 20, eps.LQA = 1e-05, eps.diff = 1e-05)
{
  index.k <- which(index.max == k)

  if (length(index.k) <= 1) {
    out <- mu.t.all[k, ]
  }
  else {
    mu.t.k <- mu.t.all[k, ]
    mu.no.k <- mu.no.penal[k, ]
    y.k <- y[index.k, ]
    alpha.k <- as.matrix(alpha[index.k, k])
    n.k <- dim(y.k)[1]
    P <- ncol(y)
    AA <- data.frame(value = c(rep(alpha.k, P)/(2 * rep(sigma.all, n.k))), group = rep(1:P, each = n.k))
    # mat1 contain postive values
    mat1 <- diag(aggregate(AA$value, by = list(Category = AA$group),  FUN = sum)$x)
    BB <- data.frame(value = AA$value * as.vector(y.k), group = rep(1:P,each = n.k))
    vec1 <- as.matrix(aggregate(BB$value, by = list(Category = BB$group),
                                FUN = sum)$x)
    combn.mat <- matrix(0, P - 1, P)
    for (s.i in 1:nrow(combn.mat)) {
      combn.mat[s.i, s.i] <- -1
      combn.mat[s.i, s.i + 1] <- 1
    }
    dist.mu.no.k <- abs(combn.mat %*% mu.no.k)
    dist.mu.no.k[dist.mu.no.k < eps.diff] <- eps.diff
    # log.tau <- c(-log(abs(combn.mat %*% mu.no.k)))
    # !!! avoid 0s in log function
    log.tau <- abs(combn.mat %*% mu.no.k)
    log.tau[log.tau < eps.diff] <- eps.diff*0.1
    log.tau <- c(-log(log.tau))
    mu.k.hat <- matrix(0, nrow = iter.num + 1, ncol = P)
    mu.k.hat[1, ] <- mu.t.k

    # s<-1
    # t1<-Sys.time()
    for (s in 1:iter.num) {
      # print(paste0("LQA iter:",s))
      dist.old.tmp0 <- abs(combn.mat %*% mu.k.hat[s, ])
      dist.old.tmp0[dist.old.tmp0 < eps.diff] <- eps.diff
      c.kk <- as.vector(sqrt(exp(log.tau - log(2 * dist.old.tmp0))))
      # B is a symetric matrix
      B <- matrix(0, P, P)
      B[1, 1] <- c(c.kk^2)[1]
      B[length(c.kk), length(c.kk) + 1] <- -c(c.kk^2)[length(c.kk)]
      B[length(c.kk) + 1, length(c.kk)] <- -c(c.kk^2)[length(c.kk)]
      B[length(c.kk) + 1, length(c.kk) + 1] <- c(c.kk^2)[length(c.kk)]
      for (i in 1:(length(c.kk) - 1)) {
        B[i, i + 1] <- -c(c.kk^2)[i]
        B[i + 1, i] <- -c(c.kk^2)[i]
        B[i + 1, i + 1] <- c(c.kk^2)[i] + c(c.kk^2)[i + 1]
      }

      # system.time(C <- flowCalcCpp(mat1 + lambda * B, vec1))
      system.time(C<-batch_calcCpp(B=B, mat1=mat1, vec1=vec1,lambda=lambda))
      mu.k.hat[(s + 1), ] <- C$Imp
      if (sum(abs(mu.k.hat[(s + 1), ] - mu.k.hat[s, ]))/(sum(abs(mu.k.hat[s, ])) + 1) < eps.LQA) {
        out <- mu.k.hat[(s + 1), ]
        out[which(abs(out) < eps.diff)] <- 0
        break
      }
      else if (s == iter.num) {
        out <- mu.k.hat[(s + 1), ]
        out[which(abs(out) < eps.diff)] <- 0
      }
    }
    # Sys.time()-t1
  }
  return(out)
}


#' modify the functions of solving linear functions:
#' 1, infinite values in B will be replaced by mean diagonal values (infinite values lead to NAs)
#' 2, participate the whole matrix into batches to increase computational time.
#' @param chr.seq the chr information for the bins
#' @param B matrix calculated in the FLCNA_LQA function
#' @param mat1 matrix calculated in the FLCNA_LQA function
#' @param vec1 vector calculated in the FLCNA_LQA function
#' @param lambda the penalty term
#' @param batch.size the number of bins in each batch
#' @return Imp the solution of the linear functions
#'
#'@export
#'
#'
# compensate function i FLCNA_LQA
batch_calcCpp<-function(chr.seq=NULL, #sequence of start and end position for all chrs
                        B,
                        mat1,
                        vec1,
                        lambda,
                        batch.size=500){
  # if chr.seq is null, use default batch size as 500
  if(is.null(chr.seq)){
    chr.seq<-seq(1,ncol(B),by=batch.size)
    if(chr.seq[length(chr.seq)]!=ncol(B)){chr.seq[length(chr.seq)]=ncol(B)}
  }
  if (chr.seq[1]!=1){chr.seq<-c(1,chr.seq)}
  if (ncol(B)-chr.seq[length(chr.seq)]<batch.size){
    chr.seq[length(chr.seq)]<-ncol(B)
  }else{chr.seq<-c(chr.seq,ncol(B))}
  if(length(which(is.infinite(B)))>0){
    B[which(is.infinite(B))]<-mean(diag(B)[which(! is.infinite(diag(B)) )])
  }


  # print(chr.seq)
  chr.seq[1]<-0 # to make chr.seq[i]+1 feasible for 1st postion
  C<-c()
  for (i in 1:(length(chr.seq)-1)){
    test.start<-chr.seq[i]+1
    test.end  <-chr.seq[i+1]
    B2<-B[test.start:test.end,test.start:test.end]
    B2[1,2]<- -B2[1,1]
    B2[2,1]<- -B2[1,1]
    B2[ncol(B2)-1,ncol(B2)]<- -B2[ncol(B2),ncol(B2)]
    B2[ncol(B2),ncol(B2)-1]<- -B2[ncol(B2),ncol(B2)]
    mat2<-mat1[test.start:test.end,test.start:test.end]
    vec2<-as.matrix(vec1[test.start:test.end])
    tmp.C2 <- flowCalcCpp(mat2 + lambda * B2, vec2)
    C<-c(C,as.vector(tmp.C2$Imp))
  }

  return(list(Imp=as.matrix(C)))
}






#' @title FLCNA analysis
#'
#' @description Simultaneous CNA detection and subclone identification using single cell DNA sequencing data.
#'
#' @param tuning A 2-dimensional vector or a matrix with 2 columns, the first column is the number of clusters \eqn{K} and the second column is the tuning parameter \eqn{\lambda} in the penalty term. If this is missing, then \code{K} and \code{lambda} must be provided.
#' @param K The number of clusters \eqn{K}.
#' @param lambda The tuning parameter \eqn{\lambda} in the penalty term. The default is 1.5.
#' @param Y A p-dimensional data matrix. Each row is an observation.
#' @param N The maximum number of iterations in the EM algorithm. The default value is 100.
#' @param kms.iter The maximum number of iterations in kmeans algorithm for generating the starting value for the EM algorithm.
#' @param kms.nstart The number of starting values in K-means.
#' @param adapt.kms A indicator of using the cluster means estimated by K-means to calculate the adaptive parameters. The default value is FALSE.
#' @param eps.diff The lower bound of pairwise difference of two mean values. Any value lower than it is treated as 0.
#' @param eps.em The lower bound for the stopping criterion in the EM algorithm.
#' @param iter.LQA The number of iterations in the estimation of cluster means by using the local quadratic approximation (LQA).
#' @param eps.LQA The lower bound for the stopping criterion in the estimation of cluster means.
#' @param cutoff Cutoff value to further control the number of CNAs besed on mean matrix from FL model. Larger cutoff value, less CNAs.
#' @param L Repeat times in the EM algorithm while outputing CNA data, defaults to 100.
#' @param model.crit The criterion used to select the number of clusters \eqn{K}. It is either `bic' for Bayesian Information Criterion or `gic' for Generalized Information Criterion.
#'
#'
#'
#' @return This function returns the esimated parameters and some statistics of the optimal model within the given \eqn{K} and \eqn{\lambda}, which is selected by BIC when \code{model.crit = 'bic'} or GIC when \code{model.crit = 'gic'}.
#' \item{K.best}{The optimal number of clusters.}
#' \item{mu.hat.best}{The estimated cluster means in the optimal model.}
#' \item{sigma.hat.best}{The estimated covariance in the optimal model.}
#' \item{alpha.hat.best}{posterior probabilities in the optimal model.}
#' \item{p.hat.best}{The estimated cluster proportions in the optimal model.}
#' \item{s.hat.best}{The clustering assignments using the optimal model.}
#' \item{lambda.best}{The value of tuning hyperparameter lambda that provide the optimal model.}
#' \item{gic.best}{The GIC of the optimal model.}
#' \item{bic.best}{The BIC of the optimal model.}
#' \item{llh.best}{The log-likelihood of the optimal model.}
#' \item{ct.mu.best}{The degrees of freedom in the cluster means of the optimal model.}
#' \item{K}{The input k values.}
#' \item{lambda}{The input lambda values.}
#' \item{mu.hat}{The estimated cluster means for each parameter combination.}
#' \item{sigma.hat}{The estimated covariance for each parameter combination.}
#' \item{p.hat}{The estimated cluster proportions for each parameter combination.}
#' \item{s.hat = s.hat}{The clustering assignments for each parameter combination.}
#' \item{gic}{The GIC values for each parameter combination.}
#' \item{bic}{The BIC values for each parameter combination.}
#' \item{llh}{The log-likelihood values for each parameter combination.}
#' \item{ct.mu}{The degrees of freedom in the cluster means for each parameter combination.}
#'
#' @examples
#' y <- matrix(rnorm(10000, 0, 0.5),10, 1000)
#' output <- FLCNA(K = c(1:2), lambda = c(2,3), y=y)
#' output
#'
#' @export
FLCNA_mod <- function(tuning=NULL, K=NULL, lambda = c(1.5), y, N = 100, kms.iter = 100, kms.nstart = 100,
                      adapt.kms = FALSE, eps.diff = 1e-5, eps.em = 1e-5, iter.LQA = 20, eps.LQA = 1e-5,
                      cutoff=0.5, L=100, model.crit = 'bic'){
  ##  tuning: a matrix with 2 columns;
  ##  1st column = K (number of clusters), positive integer;
  ##  2nd column = lambda, nonnegative real number.

  ## ----- check invalid numbers
  #if (dim(y)[1]>= dim(y)[2]){
  #  warning("rows denote cells and columns denote bins, transposing the matrix")
  #  y = t(as.matrix(y))
  #}
  
  y = t(as.matrix(y))
  if (is.null(tuning)){
    if(is.null(K) | is.null(lambda)){
      stop("Require a matrix of tuning parameters for 'tuning' or vectors for 'K' and 'lambda' ")
    }
    else{
      tuning = as.matrix(expand.grid(K, lambda))
      colnames(tuning)=c()
    }
  }

  if(is.vector(tuning) == TRUE){
    tuning = unlist(tuning)
    if(length(tuning) != 2){
      stop(" 'tuning' must be a vector with length 2 or a matrix with two columns")
    }
    else if(dim(y)[1] < tuning[1]){
      stop(" The number of clusters 'K' is greater than the number of observations in 'y' ")
    }
    else if( tuning[1] - round(tuning[1]) > .Machine$double.eps^0.5 | tuning[1] <= 0 | tuning[2] < 0){
      stop(" 'K' must be a positive interger and 'lambda' must be nonnegative ")
    }
  }

  if(is.matrix(tuning) == TRUE){
    if(dim(tuning)[2] != 2){
      stop(" 'tuning' must be a vector with length 2 or a matrix with two columns")
    }
    else if(any(dim(y)[1] < tuning[,1])){
      stop(" The number of clusters 'K' must be greater than the number of observations in 'y' ")
    }
    else if(any(tuning[,1] - round(tuning[,1]) > .Machine$double.eps^0.5 | tuning[,1] <= 0 | tuning[,2] < 0)){
      stop(" 'K' must be a positive interger and 'lambda' must be nonnegative ")
    }
  }

  if(is.matrix(tuning) == FALSE){
    stop(" 'tuning' must be a vector with length 2 or a matrix with two columns")
  }


  ### --- data dimension ---
  n = nrow(y)
  d = ncol(y)

  ### --- outputs ---
  s.hat =list()	            # clustering labels
  mu.hat = list()						# cluster means
  # CNAdata = list()          # CNA output
  sigma.hat = list()				# cluster variance
  alpha.hat = list()        # Probability for each observation to be assigned into each cluster
  p.hat = list() 				    # cluster proportion
  not.cvg = rep(0, dim(tuning)[1])       # convergence status

  ## ------------ BIC and GIC ----
  llh=c()
  ct.mu=c()
  gic=c()
  bic=c()

  for(j.tune in 1:dim(tuning)[1]){

    K1 = tuning[j.tune,1]
    lambda1 = tuning[j.tune,2]

    message(paste0("Start estimating of K=",K1,", lambda=",lambda1))
    if(K1 == 1) {
      mu.hat[[j.tune]] <- colMeans(y)
      sigma.hat[[j.tune]] <- apply(y, 2, var)
      p.hat[[j.tune]] <- 1
      s.hat[[j.tune]] <- rep(1,n)
      llh[j.tune]   <- sum(dmvnorm(x = y, mean = mu.hat[[j.tune]], sigma = diag(sigma.hat[[j.tune]]), log=TRUE))
      ct.mu[j.tune] <- sum(abs(mu.hat[[j.tune]]) > eps.diff)
      gic[j.tune]   <- -2*llh[j.tune] + log(log(n))*log(d)*(d + ct.mu[j.tune])
      bic[j.tune]   <- -2*llh[j.tune] + log(n)*(d + ct.mu[j.tune])
    }

    else if(lambda1 == 0){
      temp.out <- nopenalty(K=K1, y=y, N=N, kms.iter=kms.iter, kms.nstart=kms.nstart, eps.diff=eps.diff,
                            eps.em=eps.em, model.crit=model.crit)
      mu.hat[[j.tune]]    <- temp.out$mu.hat.best
      sigma.hat[[j.tune]] <- temp.out$sigma.hat.best
      p.hat[[j.tune]]     <- temp.out$p.hat.best
      s.hat[[j.tune]]     <- temp.out$s.hat.best
      llh[j.tune]         <- temp.out$llh.best
      ct.mu[j.tune]       <- temp.out$ct.mu.best
      gic[j.tune]         <- temp.out$gic.best
      bic[j.tune]         <- temp.out$bic.best
    }

    else{
      sigma.iter <-  matrix(0, ncol = d, nrow = N)		#each row is variances (d dimensions)
      p <- matrix(0, ncol = K1, nrow = N)		# each row is clusters proportions
      mu <- array(0, dim = c(N, K1, d))		# N=iteration, K1=each cluster , d=dimension

      ### --- start from kmeans with multiple starting values	-----
      set.seed(11)
      kms1 <- kmeans(y, centers = K1, nstart= kms.nstart, iter.max = kms.iter)
      kms1.class <- kms1$cluster

      mean0.fn <- function(k){
        if(length(y[kms1.class == k, 1]) == 1) {
          ## if there is only one data in a cluster,
          out <- y[kms1.class == k,]
        }
        else {
          out <- colMeans(y[kms1.class == k,])
        }
        return(out)
      }
      mu[1,,] <- t(sapply(1:K1, mean0.fn))
      sigma.iter[1,] <- apply(y, 2, var)
      p[1, ] <- as.numeric(table(kms1.class))/n

      ## use kmeans' results as the adaptive parameter
      if(adapt.kms == FALSE){
        temp.out <- nopenalty(K=K1, y=y, N=N, kms.iter=kms.iter, kms.nstart=kms.nstart, eps.diff=eps.diff,
                              eps.em=eps.em, model.crit=model.crit)
        mu.no.penal <- temp.out$mu.hat.best
      }
      else {
        mu.no.penal <- mu[1,,]
      }

      ## array of posterior probability computed in E-step
      alpha.temp = array(0, dim=c(N,n,K1))

      for (t in 1:(N-1)) {
        # E-step: compute all the cluster-wise density
        # temp.normal = sapply(c(1:K1), dmvnorm_log, y=y, mu=mu[t,,], sigma = diag(sigma.iter[t,]))
        temp.normal = dmvnorm_log_sapply(seq.max    = K1,
                                         y          = y,
                                         mu         = mu[t, , ],
                                         sigma.iter = sigma.iter[t, ],batch.size=50)

        alpha.fn = function(k, log.dens = temp.normal, p.temp = p[t,]){
          if(p.temp[k] ==0){out.alpha = rep(0,dim(log.dens)[1])}
          else {
            log.mix1 = sweep(log.dens,2, log(p.temp), '+')
            log.ratio1 = sweep(log.mix1, 1, log.mix1[,k], '-')
            out.alpha = 1/rowSums(exp(log.ratio1))
            out.alpha[which(rowSums(exp(log.ratio1)) == 0)] = 0
          }
          return(out.alpha)
        }

        # E-step: update posterior probabilities alpha
        alpha.temp[(t+1),, ] = sapply(c(1:K1), alpha.fn, log.dens = temp.normal, p.temp = p[t,])

        # M-step  part 1: update cluster proportions p_k
        p[(t+1), ] = colMeans(alpha.temp[(t+1),,])
        index.max <- apply(alpha.temp[(t+1),, ] , 1, which.max)

        # M-step  part 2: update cluster sigma_j, c() read matrix by columns
        sig.fn = function(index, all.combined, alpha) {  #alpha = alpha[(t+1),,]
          combined = all.combined[,index]
          y.j = combined[1:n]		# n observation's l-th dimension
          mu.j = combined[(n+1):length(combined)]			# K means' l-th dimension
          return(sum((rep(y.j, times= length(mu.j)) - rep(mu.j, each=n))^2*c(alpha))/n)
        }
        all.combined = rbind(y, mu[t,,])   # (n+K1) by d matrix
        sigma.iter[(t+1),] = sapply(c(1:d), sig.fn, all.combined=all.combined, alpha=alpha.temp[(t+1),,])

        mu[(t+1),,] = t(sapply(1:K1, FLCNA_LQA_mod, K1, index.max=index.max, mu.t.all=mu[t,,], mu.no.penal=mu.no.penal, y=y, alpha=alpha.temp[(t+1),,],
                               sigma.all=sigma.iter[(t+1),], iter.num = iter.LQA, eps.LQA=eps.LQA,
                               eps.diff=eps.diff, lambda=lambda1))

        if(sum(abs(mu[(t+1),,]-mu[t,,]))/(sum(abs(mu[t,,]))+1)+ sum(abs(sigma.iter[(t+1),]-sigma.iter[t,]))/sum(abs(sigma.iter[t,])) +
           sum(abs(p[(t+1),] - p[t,]))/sum(p[t,]) < eps.em){
          s.hat[[j.tune]] = apply(alpha.temp[(t+1),,], 1, which.max)
          mu.hat[[j.tune]] = mu[(t+1),,]							# cluster means
          sigma.hat[[j.tune]] = sigma.iter[(t+1),]					# cluster variance
          alpha.hat[[j.tune]] <- alpha.temp[(t+1),,]
          p.hat[[j.tune]] = p[(t+1),]
          break
        }
        if(t==N-1) {
          s.hat[[j.tune]] = apply(alpha.temp[(t+1),,],1,which.max)		# clustering labels
          mu.hat[[j.tune]] = mu[(t+1),,]							# cluster means
          sigma.hat[[j.tune]] = sigma.iter[(t+1),]					# cluster variance
          alpha.hat[[j.tune]] <- alpha.temp[(t+1),,]
          p.hat[[j.tune]] = p[(t+1),]
          warning(paste('not converge when lambda = ', lambda1, 'and K = ', K1, sep=''))
        }
      }      #ends of each estimates

      #CNA output
      # message(paste0("Start output CNA data for K=",K1,", lambda=",lambda1))
      # CNAdata[[j.tune]] <- CNA.out(mean.matrix=mu.hat[[j.tune]], cutoff=cutoff, L=L)

      ## ------------ BIC and GIC ----
      label.s = sort(unique(s.hat[[j.tune]]))

      ## with empty clusters
      if(length(label.s) < K1){
        llh[j.tune] = sum(table(s.hat[[j.tune]])*log(p.hat[[j.tune]][label.s]))
        for(k in label.s){
          ## llh[j.tune] = llh[j.tune] + sum(dmvnorm(y[s.hat[[j.tune]]==k,], mean=mu.hat[[j.tune]][k,], sigma = diag(sigma.hat[[j.tune]]), log=TRUE))
          llh[j.tune] = llh[j.tune] + sum(dmvnorm_sapply(y1     = y[s.hat[[j.tune]] ==  k, ],
                                                         mean1  = mu.hat[[j.tune]][k, ],
                                                         sigma1 = diag(sigma.hat[[j.tune]]),
                                                         batch.size=50) )
        }
      }
      ## no empty clusters
      else {
        llh[j.tune] = sum(table(s.hat[[j.tune]])*log(p.hat[[j.tune]]))
        for(k in 1:K1){
          ## llh[j.tune] = llh[j.tune] + sum(dmvnorm(y[s.hat[[j.tune]]==k,], mean=mu.hat[[j.tune]][k,], sigma = diag(sigma.hat[[j.tune]]), log=TRUE))
          llh[j.tune] = llh[j.tune] + sum(dmvnorm_sapply(y1     = y[s.hat[[j.tune]] ==  k, ],
                                                         mean1  = mu.hat[[j.tune]][k, ],
                                                         sigma1 = diag(sigma.hat[[j.tune]]),
                                                         batch.size=50) )
        }
      }
      ct.mu[j.tune] = sum(apply(mu.hat[[j.tune]], 2, count.mu, eps.diff = eps.diff))
      gic[j.tune] = -2*llh[j.tune] + log(log(n))*log(d)*(length(label.s)-1+d+ct.mu[j.tune])
      bic[j.tune] = -2*llh[j.tune] + log(n)*(length(label.s)-1+d+ct.mu[j.tune])
    }
  } # end of multiple lambda

  if(model.crit == 'bic'){
    index.best = which.min(bic)
  }
  else {
    index.best = which.min(gic)
  }
  gic.best = gic[index.best]
  bic.best = bic[index.best]
  llh.best = llh[index.best]
  ct.mu.best = ct.mu[index.best]

  alpha.hat.best = alpha.hat[[index.best]]
  p.hat.best = p.hat[[index.best]]
  s.hat.best = s.hat[[index.best]]
  mu.hat.best = mu.hat[[index.best]]
  sigma.hat.best = sigma.hat[[index.best]]

  K.best = tuning[index.best,1]
  lambda.best = tuning[index.best,2]

  output = structure(list(K.best = K.best, mu.hat.best=mu.hat.best, sigma.hat.best=sigma.hat.best,
                          alpha.hat.best= alpha.hat.best, p.hat.best=p.hat.best, s.hat.best=s.hat.best,
                          lambda.best=lambda.best, gic.best = gic.best, bic.best=bic.best,
                          llh.best = llh.best, ct.mu.best = ct.mu.best,
                          K = tuning[,1], lambda = tuning[,2],
                          mu.hat = mu.hat, sigma.hat = sigma.hat, p.hat = p.hat,
                          s.hat = s.hat, gic = gic, bic = bic, llh = llh, ct.mu = ct.mu))
  return(output)
}
