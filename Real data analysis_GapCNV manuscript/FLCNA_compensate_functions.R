# ------------------------------------------------------------------------------------------ #
dmvnorm_log_sapply<-function(seq.max,
                             y,
                             mu,
                             sigma.iter,
                             batch.size=50){
  nlength<-ncol(y)
  
  times<-nlength%/%batch.size
  redu <-nlength%%batch.size
  
  tmp.normal.density = sapply(c(1:seq.max), dmvnorm_log, 
                              y     = y[,1:batch.size], 
                              mu    = mu[, 1:batch.size], 
                              sigma = diag(sigma.iter[1:batch.size]))
  if (times>1){
    for (i in 2:times) {
      tmp.normal.density1 = sapply(c(1:seq.max), dmvnorm_log, 
                                   y     = y[,((i-1)*batch.size+1):(i*batch.size)], 
                                   mu    = mu[, ((i-1)*batch.size+1):(i*batch.size)], 
                                   sigma = diag(sigma.iter[((i-1)*batch.size+1):(i*batch.size)]))
      tmp.normal.density<-tmp.normal.density+tmp.normal.density1
    }
  }
  
  if (redu>0){
    # calculate dmvnorm_log for the remaining sequence
    tmp.normal.density1 = sapply(c(1:seq.max), dmvnorm_log, 
                                 y     = y[,(times*batch.size+1):nlength], 
                                 mu    = mu[, (times*batch.size+1):nlength], 
                                 sigma = diag(sigma.iter[(times*batch.size+1):nlength]))
    tmp.normal.density<-tmp.normal.density+tmp.normal.density1
  }
  
  
  return(tmp.normal.density)
  
}



# ------------------------------------------------------------------------------------------ #
dmvnorm_sapply<-function(y1,
                         mean1,
                         sigma1,
                         batch.size=50){
  
  if (is.matrix(y1)){nlength<-ncol(y1)}
  
  if(is.vector(y1)){
    nlength<-length(y1)
    y1<-matrix(y1,nrow = 1)
    mean1<-matrix(mean1,nrow = 1)
  }
  # dim(y1)
  # dim(mean1)
  times<-nlength%/%batch.size
  redu <-nlength%%batch.size
  
  tmp.normal.density = dmvnorm(x= y1[,1:batch.size], 
                               mean = mean1[1:batch.size], 
                               sigma = sigma1[1:batch.size,1:batch.size], 
                               log = TRUE)
  if (times>1){
    for (i in 2:times) {
      tmp.normal.density1 = dmvnorm(x= y1[,((i-1)*batch.size+1):(i*batch.size)], 
                                    mean = mean1[((i-1)*batch.size+1):(i*batch.size)], 
                                    sigma = sigma1[((i-1)*batch.size+1):(i*batch.size),((i-1)*batch.size+1):(i*batch.size)], 
                                    log = TRUE)
      tmp.normal.density<-tmp.normal.density+tmp.normal.density1
    }
  }
  
  if (redu>0){
    # calculate dmvnorm for the remaining sequence
    tmp.normal.density1 = dmvnorm(x= y1[,(times*batch.size+1):nlength], 
                                  mean = mean1[(times*batch.size+1):nlength], 
                                  sigma = sigma1[(times*batch.size+1):nlength,(times*batch.size+1):nlength], 
                                  log = TRUE)
    tmp.normal.density<-tmp.normal.density+tmp.normal.density1
  }
  
  return(tmp.normal.density)
  
}



# ------------------------------------------------------------------------------------------ #
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
  
  
  
  print(chr.seq)
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






# ------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------ #

FLCNA<-function (tuning = NULL, K = NULL, lambda = c(5), Y, N = 100,
                 kms.iter = 100, kms.nstart = 100, adapt.kms = FALSE, eps.diff = 1e-05, 
                 eps.em = 1e-05, iter.LQA = 20, eps.LQA = 1e-05, cutoff = 0.5, 
                 L = 100, model.crit = "bic") 
{
  y = t(as.matrix(Y))
  if (is.null(tuning)) {
    if (is.null(K) | is.null(lambda)) {
      stop("Require a matrix of tuning parameters for 'tuning' or vectors for 'K' and 'lambda' ")
    }
    else {
      tuning = as.matrix(expand.grid(K, lambda))
      colnames(tuning) = c()
    }
  }
  if (is.vector(tuning) == TRUE) {
    tuning = unlist(tuning)
    if (length(tuning) != 2) {
      stop(" 'tuning' must be a vector with length 2 or a matrix with two columns")
    }
    else if (dim(y)[1] < tuning[1]) {
      stop(" The number of clusters 'K' is greater than the number of observations in 'y' ")
    }
    else if (tuning[1] - round(tuning[1]) > .Machine$double.eps^0.5 | 
             tuning[1] <= 0 | tuning[2] < 0) {
      stop(" 'K' must be a positive interger and 'lambda' must be nonnegative ")
    }
  }
  if (is.matrix(tuning) == TRUE) {
    if (dim(tuning)[2] != 2) {
      stop(" 'tuning' must be a vector with length 2 or a matrix with two columns")
    }
    else if (any(dim(y)[1] < tuning[, 1])) {
      stop(" The number of clusters 'K' must be greater than the number of observations in 'y' ")
    }
    else if (any(tuning[, 1] - round(tuning[, 1]) > .Machine$double.eps^0.5 | 
                 tuning[, 1] <= 0 | tuning[, 2] < 0)) {
      stop(" 'K' must be a positive interger and 'lambda' must be nonnegative ")
    }
  }
  if (is.matrix(tuning) == FALSE) {
    stop(" 'tuning' must be a vector with length 2 or a matrix with two columns")
  }
  n = nrow(y)
  d = ncol(y)
  s.hat = list()
  mu.hat = list()
  sigma.hat = list()
  alpha.hat = list()
  p.hat = list()
  not.cvg = rep(0, dim(tuning)[1])
  llh = c()
  ct.mu = c()
  gic = c()
  bic = c()
  # j.tune<-1
  for (j.tune in 1:dim(tuning)[1]) {
    K1 = tuning[j.tune, 1]
    lambda1 = tuning[j.tune, 2]
    message(paste0("Start estimating of K=", K1, ", lambda=", 
                   lambda1))
    if (K1 == 1) {
      mu.hat[[j.tune]] <- colMeans(y)
      sigma.hat[[j.tune]] <- apply(y, 2, var)
      p.hat[[j.tune]] <- 1
      s.hat[[j.tune]] <- rep(1, n)
      llh[j.tune] <- sum(dmvnorm(x = y, mean = mu.hat[[j.tune]], 
                                 sigma = diag(sigma.hat[[j.tune]]), log = TRUE))
      ct.mu[j.tune] <- sum(abs(mu.hat[[j.tune]]) > eps.diff)
      gic[j.tune] <- -2 * llh[j.tune] + log(log(n)) * log(d) * 
        (d + ct.mu[j.tune])
      bic[j.tune] <- -2 * llh[j.tune] + log(n) * (d + ct.mu[j.tune])
    }
    else if (lambda1 == 0) {
      temp.out <- nopenalty(K = K1, y = y, N = N, kms.iter = kms.iter, 
                            kms.nstart = kms.nstart, eps.diff = eps.diff, 
                            eps.em = eps.em, model.crit = model.crit)
      mu.hat[[j.tune]] <- temp.out$mu.hat.best
      sigma.hat[[j.tune]] <- temp.out$sigma.hat.best
      p.hat[[j.tune]] <- temp.out$p.hat.best
      s.hat[[j.tune]] <- temp.out$s.hat.best
      llh[j.tune] <- temp.out$llh.best
      ct.mu[j.tune] <- temp.out$ct.mu.best
      gic[j.tune] <- temp.out$gic.best
      bic[j.tune] <- temp.out$bic.best
    }
    else {
      sigma.iter <- matrix(0, ncol = d, nrow = N)
      p <- matrix(0, ncol = K1, nrow = N)
      mu <- array(0, dim = c(N, K1, d))
      set.seed(11)
      kms1 <- kmeans(y, centers = K1, nstart = kms.nstart, 
                     iter.max = kms.iter)
      kms1.class <- kms1$cluster
      mean0.fn <- function(k) {
        if (length(y[kms1.class == k, 1]) == 1) {
          out <- y[kms1.class == k, ]
        }
        else {
          out <- colMeans(y[kms1.class == k, ])
        }
        return(out)
      }
      mu[1, , ] <- t(sapply(1:K1, mean0.fn))
      sigma.iter[1, ] <- apply(y, 2, var)
      p[1, ] <- as.numeric(table(kms1.class))/n
      
      # CHECK1 Temporally end here
      if (adapt.kms == FALSE) {
        temp.out <- nopenalty(K=K1,
                              y=y,
                              N=N,
                              kms.iter=kms.iter,
                              kms.nstart=kms.nstart,
                              eps.diff=eps.diff,
                              eps.em=eps.em,
                              model.crit=model.crit)
        mu.no.penal <- temp.out$mu.hat.best
      } else {
        mu.no.penal <- mu[1, , ]
      }
      alpha.temp = array(0, dim = c(N, n, K1))
      # t = 1
      for (t in 1:(N - 1)) {
        print(paste0("FLCNA iter: ",t))
        # temp.normal = sapply(c(1:K1), dmvnorm_log, y = y, 
        #                      mu = mu[t, , ], sigma = diag(sigma.iter[t,]))
        temp.normal = dmvnorm_log_sapply(seq.max    = K1,
                                         y          = y,
                                         mu         = mu[t, , ],
                                         sigma.iter = sigma.iter[t, ],batch.size=50)
        
        alpha.fn = function(k, log.dens = temp.normal, 
                            p.temp = p[t, ]) {
          if (p.temp[k] == 0) {
            out.alpha = rep(0, dim(log.dens)[1])
          } else {
            log.mix1 = sweep(log.dens, 2, log(p.temp), 
                             "+")
            log.ratio1 = sweep(log.mix1, 1, log.mix1[, 
                                                     k], "-")
            out.alpha = 1/rowSums(exp(log.ratio1))
            out.alpha[which(rowSums(exp(log.ratio1)) == 
                              0)] = 0
          }
          return(out.alpha)
        }
        alpha.temp[(t + 1), , ] = sapply(c(1:K1), alpha.fn, 
                                         log.dens = temp.normal, p.temp = p[t, ])
        p[(t + 1), ] = colMeans(alpha.temp[(t + 1), , 
                                           ])
        index.max <- apply(alpha.temp[(t + 1), , ], 1, 
                           which.max)
        sig.fn = function(index, all.combined, alpha) {
          combined = all.combined[, index]
          y.j = combined[1:n]
          mu.j = combined[(n + 1):length(combined)]
          return(sum((rep(y.j, times = length(mu.j)) - 
                        rep(mu.j, each = n))^2 * c(alpha))/n)
        }
        all.combined = rbind(y, mu[t, , ])
        sigma.iter[(t + 1), ] = sapply(c(1:d), sig.fn, 
                                       all.combined = all.combined, alpha = alpha.temp[(t + 
                                                                                          1), , ])
        mu[(t + 1), , ] = t(sapply(1:K1, FLCNA_LQA, K1, 
                                   index.max = index.max, mu.t.all = mu[t, , ], 
                                   mu.no.penal = mu.no.penal, y = y, alpha = alpha.temp[(t + 1), , ], 
                                   sigma.all = sigma.iter[(t + 1),], 
                                   iter.num = iter.LQA, 
                                   eps.LQA = eps.LQA, 
                                   eps.diff = eps.diff, 
                                   lambda = lambda1))
        
        if (sum(abs(mu[(t + 1), , ] - mu[t, , ]))/(sum(abs(mu[t, , ])) + 1) + sum(abs(sigma.iter[(t + 1), ] - 
                                                                                      sigma.iter[t, ]))/sum(abs(sigma.iter[t, ])) + 
            sum(abs(p[(t + 1), ] - p[t, ]))/sum(p[t, ]) < 
            eps.em) {
          s.hat[[j.tune]] = apply(alpha.temp[(t + 1), 
                                             , ], 1, which.max)
          mu.hat[[j.tune]] = mu[(t + 1), , ]
          sigma.hat[[j.tune]] = sigma.iter[(t + 1), ]
          alpha.hat[[j.tune]] <- alpha.temp[(t + 1), 
                                            , ]
          p.hat[[j.tune]] = p[(t + 1), ]
          break
        }
        if (t == N - 1) {
          s.hat[[j.tune]] = apply(alpha.temp[(t + 1), 
                                             , ], 1, which.max)
          mu.hat[[j.tune]] = mu[(t + 1), , ]
          sigma.hat[[j.tune]] = sigma.iter[(t + 1), ]
          alpha.hat[[j.tune]] <- alpha.temp[(t + 1), 
                                            , ]
          p.hat[[j.tune]] = p[(t + 1), ]
          warning(paste("not converge when lambda = ", 
                        lambda1, "and K = ", K1, sep = ""))
        }
      }
      label.s = sort(unique(s.hat[[j.tune]]))
      if (length(label.s) < K1) {
        llh[j.tune] = sum(table(s.hat[[j.tune]]) * log(p.hat[[j.tune]][label.s]))
        for (k in label.s) {
          # llh[j.tune] = llh[j.tune] + sum(dmvnorm(y[s.hat[[j.tune]] == k, ], 
          #                                         mean = mu.hat[[j.tune]][k, ], 
          #                                         sigma = diag(sigma.hat[[j.tune]]),
          #                                         log = TRUE))
          llh[j.tune] = llh[j.tune] + sum(dmvnorm_sapply(y1     = y[s.hat[[j.tune]] ==  k, ],
                                                         mean1  = mu.hat[[j.tune]][k, ],
                                                         sigma1 = diag(sigma.hat[[j.tune]]),
                                                         batch.size=50) )
        }
      }
      else {
        llh[j.tune] = sum(table(s.hat[[j.tune]]) * log(p.hat[[j.tune]]))
        for (k in 1:K1) {
          # llh[j.tune] = llh[j.tune] + sum(dmvnorm(y[s.hat[[j.tune]] ==  k, ],
          #                                         mean = mu.hat[[j.tune]][k, ], 
          #                                         sigma = diag(sigma.hat[[j.tune]]),
          #                                         log = TRUE))
          llh[j.tune] = llh[j.tune] + sum(dmvnorm_sapply(y1     = y[s.hat[[j.tune]] ==  k, ],
                                                         mean1  = mu.hat[[j.tune]][k, ],
                                                         sigma1 = diag(sigma.hat[[j.tune]]),
                                                         batch.size=50) )
        }
      }
      ct.mu[j.tune] = sum(apply(mu.hat[[j.tune]], 2, count.mu, 
                                eps.diff = eps.diff))
      gic[j.tune] = -2 * llh[j.tune] + log(log(n)) * log(d) * 
        (length(label.s) - 1 + d + ct.mu[j.tune])
      bic[j.tune] = -2 * llh[j.tune] + log(n) * (length(label.s) - 
                                                   1 + d + ct.mu[j.tune])
    }
  }
  if (model.crit == "bic") {
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
  K.best = tuning[index.best, 1]
  lambda.best = tuning[index.best, 2]
  output = structure(list(K.best = K.best, mu.hat.best = mu.hat.best, 
                          sigma.hat.best = sigma.hat.best, alpha.hat.best = alpha.hat.best, 
                          p.hat.best = p.hat.best, s.hat.best = s.hat.best, lambda.best = lambda.best, 
                          gic.best = gic.best, bic.best = bic.best, llh.best = llh.best, 
                          ct.mu.best = ct.mu.best, K = tuning[, 1], lambda = tuning[, 
                                                                                    2], mu.hat = mu.hat, sigma.hat = sigma.hat, p.hat = p.hat, 
                          s.hat = s.hat, gic = gic, bic = bic, llh = llh, ct.mu = ct.mu))
  return(output)
}



# -------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------- #
nopenalty<-function (K, y, N = 100, kms.iter = 100, kms.nstart = 100, eps.diff = 1e-05, 
                     eps.em = 1e-05, model.crit = "gic") 
{
  if (missing(K) | any(dim(y)[1] < K)) {
    stop(" The number of clusters 'K' must be greater than the number of observations in 'y' ")
  }
  n = nrow(y)
  D = ncol(y)
  s.hat = list()
  mu.hat = list()
  sigma.hat = list()
  p.hat = list()
  llh = c()
  ct.mu = c()
  gic = c()
  bic = c()
  # j.tune<-1
  for (j.tune in 1:length(K)) {
    if (K[j.tune] == 1) {
      mu.hat[[j.tune]] <- colMeans(y)
      sigma.hat[[j.tune]] <- apply(y, 2, var)
      p.hat[[j.tune]] <- 1
      s.hat[[j.tune]] <- rep(1, n)
      llh[j.tune] <- sum(dmvnorm(x = y, mean = mu.hat[[j.tune]], 
                                 sigma = diag(sigma.hat[[j.tune]]), log = TRUE))
      ct.mu[j.tune] <- sum(abs(mu.hat[[j.tune]]) > eps.diff)
      gic[j.tune] <- -2 * llh[j.tune] + log(log(n)) * log(D) * 
        (D + ct.mu[j.tune])
      bic[j.tune] <- -2 * llh[j.tune] + log(n) * (D + ct.mu[j.tune])
    }
    else {
      sigma.iter = matrix(0, ncol = D, nrow = N)
      p = matrix(0, ncol = K[j.tune], nrow = N)
      mu = array(0, dim = c(N, K[j.tune], D))
      set.seed(11)
      kms1 = kmeans(y, centers = K[j.tune], nstart = kms.nstart, 
                    iter.max = kms.iter)
      kms1.class = kms1$cluster
      mean0.fn = function(k) {
        if (length(y[kms1.class == k, 1]) == 1) {
          out = y[kms1.class == k, ]
        }
        else {
          out = colMeans(y[kms1.class == k, ])
        }
        return(out)
      }
      mu[1, , ] = t(sapply(1:K[j.tune], mean0.fn))
      sigma.iter[1, ] = apply(y, 2, var)
      p[1, ] = as.numeric(table(kms1.class))/n
      alpha.temp = array(0, dim = c(N, n, K[j.tune]))
      # t<-1
      # setwd("C:\\Users\\xuanxuan\\Dropbox\\2021_Researches\\CNV detection cooperation project\\codes\\20220824_merge_datasets\\check_FLCNA")
      # saveRDS(y,"y.chk3.RDS")
      # saveRDS(mu,"mu.chk3.RDS")
      # saveRDS(sigma.iter,"sigma.iter.chk3.RDS")
      # t<-1
      for (t in 1:(N - 1)) {
        # temp.normal = sapply(c(1:K[j.tune]), dmvnorm_log, 
        #                      y = y, mu = mu[t, , ], sigma = diag(sigma.iter[t, ]))
        temp.normal = dmvnorm_log_sapply(seq.max=K[j.tune],y = y, mu = mu[t, , ], sigma.iter = sigma.iter[t, ],batch.size=50)
        alpha.fn = function(k, log.dens = temp.normal, 
                            p.temp = p[t, ]) {
          if (p.temp[k] == 0) {
            out.alpha = rep(0, dim(log.dens)[1])
          }
          else {
            log.mix1 = sweep(log.dens, 2, log(p.temp), 
                             "+")
            log.ratio1 = sweep(log.mix1, 1, log.mix1[, 
                                                     k], "-")
            out.alpha = 1/rowSums(exp(log.ratio1))
            out.alpha[which(rowSums(exp(log.ratio1)) ==  0)] = 0
          }
          return(out.alpha)
        }
        alpha.temp[(t + 1), , ] = sapply(c(1:K[j.tune]),  alpha.fn, log.dens = temp.normal, p.temp = p[t, ])
        p[(t + 1), ] = colMeans(alpha.temp[(t + 1), ,  ])
        sig.est.fn = function(index, all.combined, alpha) {
          combined = all.combined[, index]
          y.l = combined[1:n]
          mu.l = combined[(n + 1):length(combined)]
          return(sum((rep(y.l, times = K[j.tune]) - rep(mu.l,  each = n))^2 * c(alpha))/n)
        }
        all.combined = rbind(y, mu[t, , ])
        sigma.iter[(t + 1), ] = sapply(c(1:D), sig.est.fn, 
                                       all.combined = all.combined, alpha = alpha.temp[(t + 1), , ])
        mu.est.fn = function(k.index, y, alpha) {
          mu.k1 = alpha[, k.index] %*% y/sum(alpha[, k.index])
          return(mu.k1)
        }
        mu[(t + 1), , ] = t(sapply(1:K[j.tune], mu.est.fn, y = y, alpha = alpha.temp[(t + 1), , ]))
        CRTERIA1 = sum(abs(mu[(t + 1),,]-mu[t,,]))
        CRTERIA2 = sum(abs(mu[t,,])) + 1
        CRTERIA3 = sum(abs(sigma.iter[(t+1),]-sigma.iter[t,]))
        CRTERIA4 = sum(abs(sigma.iter[t,]))
        CRTERIA5 = sum(abs(p[(t+1),]-p[t,]))
        CRTERIA6 = sum(p[t, ])
        CRTERIA  = CRTERIA1/CRTERIA2+CRTERIA3/CRTERIA4+CRTERIA5/CRTERIA6
        # CRTERIA = sum(abs(mu[(t + 1),,] - mu[t,,]))/(sum(abs(mu[t,,])) + 1) + sum(abs(sigma.iter[(t+1),]- sigma.iter[t, ]))/sum(abs(sigma.iter[t, ])) + sum(abs(p[(t + 1), ] - p[t, ]))/sum(p[t, ])
        if (CRTERIA<eps.em) {
          s.hat[[j.tune]] = apply(alpha.temp[(t + 1),  , ], 1, which.max)
          mu.hat[[j.tune]] = mu[(t + 1), , ]
          sigma.hat[[j.tune]] = sigma.iter[(t + 1), ]
          p.hat[[j.tune]] = p[(t + 1), ]
          break
        }
        if (t == (N-1)) {
          s.hat[[j.tune]] = apply(alpha.temp[(t + 1),  , ], 1, which.max)
          mu.hat[[j.tune]] = mu[(t + 1), , ]
          sigma.hat[[j.tune]] = sigma.iter[(t + 1), ]
          p.hat[[j.tune]] = p[(t + 1), ]
          warning(paste("not converge when lambda = 0 and K=", 
                        K[j.tune], sep = ""))
        }
      }
      label.s = sort(unique(s.hat[[j.tune]]))
      if (length(label.s) < K[j.tune]) {
        llh[j.tune] = sum(table(s.hat[[j.tune]]) * log(p.hat[[j.tune]][label.s]))
        for (k in label.s) {
          # llh[j.tune] = llh[j.tune] + sum(dmvnorm(y[s.hat[[j.tune]] ==  k, ],
          #                                         mean = mu.hat[[j.tune]][k, ], 
          #                                         sigma = diag(sigma.hat[[j.tune]]), 
          #                                         log = TRUE))
          llh[j.tune] = llh[j.tune] + sum(dmvnorm_sapply(y1     =y[s.hat[[j.tune]] == k, ],
                                                         mean1  = mu.hat[[j.tune]][k, ],
                                                         sigma1 = diag(sigma.hat[[j.tune]]),
                                                         batch.size=50) )
        }
      }
      else {
        llh[j.tune] = sum(table(s.hat[[j.tune]]) * log(p.hat[[j.tune]]))
        for (k in 1:K[j.tune]) {
          # llh[j.tune] = llh[j.tune] + sum(dmvnorm(y[s.hat[[j.tune]] == k, ],
          #                                         mean = mu.hat[[j.tune]][k, ],
          #                                         sigma = diag(sigma.hat[[j.tune]]),
          #                                         log = TRUE))
          
          llh[j.tune] = llh[j.tune] + sum(dmvnorm_sapply(y1     =y[s.hat[[j.tune]] == k, ],
                                                         mean1  = mu.hat[[j.tune]][k, ],
                                                         sigma1 = diag(sigma.hat[[j.tune]]),
                                                         batch.size=50) )
        }
      }
      ct.mu[j.tune] = sum(apply(mu.hat[[j.tune]], 2, count.mu, 
                                eps.diff = eps.diff))
      gic[j.tune] = -2 * llh[j.tune] + log(log(n)) * log(D) * 
        (length(label.s) - 1 + D + ct.mu[j.tune])
      bic[j.tune] = -2 * llh[j.tune] + log(n) * (length(label.s) - 1 + D + ct.mu[j.tune])
    }
  }
  if (model.crit == "bic") {
    index.best = which.min(bic)
  } else {
    index.best = which.min(gic)
  }
  gic.best = gic[index.best]
  bic.best = bic[index.best]
  llh.best = llh[index.best]
  ct.mu.best = ct.mu[index.best]
  p.hat.best = p.hat[[index.best]]
  s.hat.best = s.hat[[index.best]]
  mu.hat.best = mu.hat[[index.best]]
  sigma.hat.best = sigma.hat[[index.best]]
  K.best = K[index.best]
  output = structure(list(K.best = K.best, mu.hat.best = mu.hat.best, 
                          sigma.hat.best = sigma.hat.best, p.hat.best = p.hat.best, 
                          s.hat.best = s.hat.best, gic.best = gic.best, bic.best = bic.best, 
                          llh.best = llh.best, ct.mu.best = ct.mu.best), class = "parse_fit")
  return(output)
}





# -------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------- #
FLCNA_LQA<-function (k, K1, index.max, y, mu.t.all, mu.no.penal, sigma.all,
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
    # avoid 0s in log function
    log.tau <- abs(combn.mat %*% mu.no.k)
    log.tau[log.tau < eps.diff] <- eps.diff*0.1
    log.tau <- c(-log(log.tau))
    mu.k.hat <- matrix(0, nrow = iter.num + 1, ncol = P)
    mu.k.hat[1, ] <- mu.t.k
    
    # s<-1
    # t1<-Sys.time()
    for (s in 1:iter.num) {
      print(paste0("LQA iter:",s))
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




# =========================================================================== #
# unrestrict the maximum number of change points
CNA.out<-function (mean.matrix, ref, cutoff = 0.35, L = 100,max.cp=10000) 
{
  CNAdata <- vector(mode = "list", length = nrow(mean.matrix))
  # g<-1
  for (g in 1:nrow(mean.matrix)) {
    cp.index <- 1 + which(abs(diff(mean.matrix[g, ], 1)) >=  cutoff)
    if ((length(cp.index) < 1) | (length(cp.index) > max.cp)) {
      next
    }
    x.inv <- try(res <- CNAcluster(Y = mean.matrix[g, ], 
                                   cp = cp.index, L), silent = TRUE)
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




dmvnorm_log <- function(index, mu, sigma, y) {
  
  return(mvtnorm::dmvnorm(y, mu[index,], sigma, log=TRUE))
}


#' @title count.mu 
#' @description Computing the number of unique cluster means for each dimension, which was used in computing BIC or GIC.
#'
#' @param mu.j Mean vector.
#' @param eps.diff Lower bound of mean difference.
count.mu <- function(mu.j, eps.diff) {
  temp.dist <- as.matrix(dist(mu.j, method = 'manhattan'))
  ct <- length(mu.j[abs(mu.j)>eps.diff])
  ## initial counts (nonzero elements)
  ## --- if exists same means
  temp.dist[upper.tri(temp.dist, diag = T)] <- NA
  if(any(temp.dist < eps.diff, na.rm = T)) {
    temp.index <- which(temp.dist < eps.diff, arr.ind = TRUE)
    temp1 <- mu.j
    ## --- truncated distance so means are not exactly same, make them equal
    for(i in 1:dim(temp.index)[1]){
      temp1[temp.index[i,]] <- min(temp1[temp.index[i,]])
    }
    ct <- length(unique(temp1[abs(temp1)>eps.diff]))
  }
  return(ct)
}

