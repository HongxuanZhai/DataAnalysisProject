library(mixtools)
library(parallel)

acq.data.temp <- read.csv("~/Documents/DAnew/acquisition.csv", header = FALSE, sep = "")$V1
# acq.data.temp <- read.csv("~/Projects/DA/acquisition.csv", header = FALSE, sep = "")$V1
acq.data <- acq.data.temp[which(acq.data.temp > 0)]
remove(acq.data.temp)

#' Function to a list of hyperparameter list
#' in constrained prior model 
#'
#' @param a expanded grid of hyperparameters, should be in order
#' of (v0, ss0, a, A, g1, g2, c, d)
#' 
#' @return a list of list representing each different
#' hyperparameter configuration
c.gethp <- function(grid) {
  hplist <- list()
  nRow <- nrow(grid)
  for (i in seq(nRow)) {
    hplist[[i]] <- list(v0 = grid[i, 1],
                        ss0 = grid[i, 2],
                        a = grid[i, 3],
                        A = grid[i, 4],
                        g1 = grid[i, 5],
                        g2 = grid[i, 6],
                        c = grid[i, 7],
                        d = grid[i, 8])
  }
  return(hplist)
}


#' Function to a list of hyperparameter list
#' in flexible prior model 
#'
#' @param a expanded grid of hyperparameters, should be in order
#' of (v0, ss0, a, A, g1, g2, c, d, df, v1)
#' 
#' @return a list of list representing each different
#' hyperparameter configuration
h.gethp <- function(grid) {
  hplist <- list()
  nRow <- nrow(grid)
  for (i in seq(nRow)) {
    hplist[[i]] <- list(v0 = grid[i, 1],
                        ss0 = grid[i, 2],
                        a = grid[i, 3],
                        A = grid[i, 4],
                        g1 = grid[i, 5],
                        g2 = grid[i, 6],
                        c = grid[i, 7],
                        d = grid[i, 8],
                        df = grid[i, 9],
                        v1 = grid[i, 10])
  }
  return(hplist)
}

#' Function to get log-scaled training samples 
#'
#' @param rawdata untransformed data, default is acquisition dataset
#' @param seed random seed to control randomness of sampled training sets
#' @param size sample size of one training set, default is 500
#' 
#' @return a list with log-scaled training sample and its index in rawdata
get.logdata <- function(seed, rawdata = acq.data, size = 500) {
  set.seed(seed)
  ind <- sample(seq(rawdata), size)
  return(list(data = log(rawdata[ind]),
              index = ind))
}

#' Function to get truncated predictions from KDE estimates
#'
#' @param logdata log-scaled data serving as input to "density(.)" function
#' @param npred number of predictions drawn from kde estimate
#' @param maxval knowledge of maximum of population
#'
#' @return a set of predictions of length npred from KDE estimate on the
#' log-scaled data.
pred.kde <- function(logdata, npred, maxval = max(acq.data)) {
  dens <- density(logdata)
  samsize <- length(logdata)
  pred <- rep(NA, npred)
  for (i in seq(npred)) {
    acc <- 0
    while (!acc) {
      ind <- sample(seq(samsize), 1)
      ## first draw an index
      prop <- exp(rnorm(1, logdata[ind], dens$bw))
      ## propose a lognormal draw with sd being bw 
      if(prop <= maxval) {
        pred[i] <- prop
        acc <- 1
      } else {
        acc <- 0
      }
    }
  }
  return(pred)
}

#' Function to get estimation of the mean by C.L.T. 
#'
#' @param seed.seq random seed sequence to control randomness of sampled training sets
#' 
#' @return a vector of sample means w.r.t different samples
clt.est <- function(seed.seq){
  cltmean <- rep(NA, length(seed.seq))
  for (i in seq(length(seed.seq))) {
    sam.data <- get.logdata(seed.seq[i])
    cltmean[i] <- mean(exp(sam.data$data))
  }
  return(cltmean)
}

#' Function to get estimation of the mean by bootstrapping. 
#'
#' @param seed.seq random seed sequence to control randomness of sampled training sets
#' @param nbs sample size of bootstrapping, must be equal to sample size
#' 
#' @return a vector of bootstrapping mean for each sample
bs.est <- function(seed.seq, nbs = 500) {
  bsmean <- rep(NA, length(seed.seq))
  for (i in seq(length(seed.seq))) {
    sam.data <- get.logdata(seed.seq[i])
    bs <- sample(seq(length(sam.data$data)), nbs, replace = TRUE)
    bsmean[i] <- mean(exp(sam.data$data[bs]))
  }
  return(bsmean)
}

#' Function to get truncated predictions from fittings 
#' of mixtools package, which implements EM algorithm
#' 
#' @param mixfit an fitted object from mixtools
#' @param npred prespecified number of predictions
#' @param maxval knowledge of maximum of population
#' 
#' @return a scalar; mean of the predictions
pred.em <- function(mixfit, npred, maxval = max(acq.data)) {
  pred <- rep(NA, npred)
  for (i in seq(npred)) {
    acc <- 0
    while(!acc){
      ind <- sample(seq(length(mixfit$mu)), 1, prob = mixfit$lambda)
      prop <- exp(rnorm(1, mixfit$mu[ind], mixfit$sigma[ind]))
      if (prop <= maxval){
        pred[i] <- prop
        acc <- 1
      } else {
        acc <- 0
      }
    }
  }
  return(mean(pred))
}

#' Function to calculate multiple predicted mean using 
#' predictions from mixtools
#' 
#' @param seed.seq random seed sequence to control randomness 
#' of sampled training sets
#' @param npred prespecified number of predictions for each 
#' training sample
#' @param ncluster prespecified number of cluster in finite
#' mixture model in mixtools, default as 10 for fast computing speed
#' 
#' @return a vector of mean estimates for each different samples
multi.pred.em <- function(seed.seq, npred, ncluster = 10) {
  predmean <- rep(NA, length(seed.seq))
  for (i in seq(length(seed.seq))) {
    sam <- get.logdata(seed.seq[i])
    initial <- kmeans(sam$data, ncluster)$centers
    fit <- normalmixEM(sam$data, mu = initial, maxit = 55000)
    predmean[i] <- pred.em(fit, npred)
  }
  return(predmean)
}


#' Function to calculate the L2 prediction risks
#' 
#' @param parameters are truevalue of the population mean;
#' predictions from np (constrained or flexible) models; predictions
#' from KDE; predictions from C.L.T; predictions from finite
#' mixture model; predictions from bootstrapping
#' 
#' @return a list of risks
cal.risk <- function(trueval, nppred, kdepred, cltpred, fmm, bs) {
  emrisk <- sqrt(sum((trueval-fmm)^2))
  bsrisk <- sqrt(sum((trueval-bs)^2))
  nprisk <- sqrt(sum((trueval-nppred)^2))
  kderisk <- sqrt(sum((trueval-kdepred)^2))
  cltrisk <- sqrt(sum((trueval-cltpred)^2))
  return(list(np = mean(nprisk),kde = mean(kderisk),
              fmm = mean(emrisk),
              bs = mean(bsrisk),
              clt = mean(cltrisk)))
}

#' Function to fit NP mixture model with constrained prior
#' 
#' @param hyperpar: hyperparameter list
#' @param niter: number of MCMC iterations
#' @param dat: log-scaled training dataset
#' @param n.sample: number of posterior predictions to draw at each MCMC iteration
#' @param n.cluster: truncation level specified to approximate stick breaking process
#' @param burnin1: burnin of MCMC; number of steps for shufflings; default as 3000
#' @param burnin2: burnin after shufflings
#' @param thinin: thinin of MCMC 
#' 
#' @return a list of all posterior draws and posterior prediction matrix
cDP_mixture_lnorm_ro <- function(hyperpar, niter, dat, n.sample, 
                                 n.cluster, burnin1 = 1000, burnin2 = 1000,
                                 thinin = 1) {
  
  # number of mixing components
  N = n.cluster
  
  # hyperparameter specification
  
  ## hp for phi
  nu0 <- (hyperpar$v0) * 2
  ss0 <- (hyperpar$ss0) * 2
  
  ## hp for mu0
  a <- hyperpar$a
  A <- hyperpar$A
  
  ## hp for kappa0
  
  w <- (hyperpar$g1) * 2
  W <- (hyperpar$g2) * 2
  
  ## hp for alpha
  c <- hyperpar$c
  d <- hyperpar$d
  
  initial_val <- function() {
    km <- kmeans(dat, N)
    # sort the mean, variance and mixing weights in decreasing order
    km_ordered <- sort(km$centers)
    km_cluster_ordered <- rep(NA, length(km$cluster))
    km_phi_ordered <- rep(NA, length(km$withinss))
    km_size_ordered <- rep(NA, length(km$size))
    for (i in seq(length(km_ordered))) {
      index <- which(km$centers == km_ordered[i])
      km_size_ordered[i] <- km$size[index]
      km_phi_ordered[i] <- 1 / km$withinss[index]
      index2 <- which(km$cluster == index)
      km_cluster_ordered[index2] <- i
    }
    km_p_ordered <- km_size_ordered / length(dat)
    return(list(P = km_p_ordered,
                K = km_cluster_ordered,
                mu_k = km_ordered,
                phi_k = km_phi_ordered,
                alpha = rgamma(1, shape = c, rate = d),
                mu0 = rnorm(1, a, sqrt(A)),
                kappa0 = rgamma(1, shape = 0.5 * w, rate = 0.5 * W)))
  }
  
  # Marsaglia and Tsang's method for sampling log(beta)
  MT_loggam <- function(alpha){
    if (alpha >= 1){
      ua <- alpha
      add = FALSE
    } else {
      ua <- alpha + 1
      add = TRUE
    }
    flag <- 1
    while (flag) {
      d <- ua - 1/3
      c <- 1 / sqrt(9*d)
      Z <- rnorm(1)
      if (Z > -1 / c){
        V <- (1 + c*Z)^3
        U <- runif(1)
        flag = (log(U) > (0.5*Z^2+d-d*V+d*log(V)))
      }
    }
    if(add) {
      return(log(d) +log(V) +(1 / (alpha))*log(runif(1)))
    } else{
      return(log(d)+ log(V))
    }
  }
  
  # Marsaglia and Tsang's method for sampling log(beta)
  lbeta_from_lgam <- function(alpha, beta){
    A <- MT_loggam(alpha)
    B <- MT_loggam(beta)
    A2 <- A - max(A, B)
    B2 <- B - max(A, B)
    A3 <- A2 - log1p(exp(min(A2, B2)))
    B3 <- B2 - log1p(exp(min(A2, B2)))
    return(c(A3, B3))
  }
  
  # gibbs sampler update function
  
  update.Z <- function(K, mu0, kappa0) {
    unique.K <- unique(K)
    # update Z_k for k not in K* (unique.K)
    # and for k in K*
    
    # preallocation
    mu_k <- rep(NA, N)
    phi_k <- rep(NA, N)
    Z_matrix <- matrix(NA, nrow = N, ncol = 2)
    
    for (i in seq(N)) {
      if(sum(i == unique.K)) {
        # statistics
        n_j.star <- sum(i == K)
        i.index <- which(i == K)
        Mpost <- (kappa0*mu0 + sum(dat[i.index])) / (kappa0 + n_j.star)
        Cpost <- kappa0 + n_j.star
        Bpost <- 0.5*(ss0 + kappa0*mu0^2 - Cpost*Mpost^2 + sum(dat[i.index]^2))
        # update
        phi_k[i] <- rgamma(1, (nu0 + n_j.star) / 2,
                           rate = Bpost)
        mu_k[i] <- rnorm(1, Mpost, sqrt(1 / Cpost / phi_k[i]))
      } else {
        phi_k[i] <- rgamma(1, 0.5 * nu0, rate = 0.5 * ss0)
        mu_k[i] <- rnorm(1, mu0, sqrt(1 / kappa0 / phi_k[i]))
      }
    }
    Z_matrix[,1] <- mu_k
    Z_matrix[,2] <- phi_k
    return(Z_matrix)
  }
  
  update.K <- function(mu_k, phi_k, P) {
    # preallocation
    K <- rep(NA, length(dat))
    for (i in seq(dat)) {
      prob <- P * dnorm(dat[i], mu_k, sqrt(1 / phi_k))
      K[i] <- sample(seq(N), size = 1, prob = prob)
    }
    return(K)
  }
  
  update.P <- function(K, alpha) {
    # preallocation
    V_k <- rep(NA, N)
    V_k_minus <- rep(NA, N-1)
    P <- rep(NA, N)
    # calculate M_k's
    M_vector <- sapply(c(1 : N), function(i, x) 
    {sum(i == x)}, x = K)
    for (i in seq(N-1)) {
      # statistic
      M_k <- M_vector[i]
      M_k_sum <- sum(M_vector[(i + 1): N])
      temp <- lbeta_from_lgam(1 + M_k, alpha + M_k_sum)
      # V_k[i] <- rbeta(1, 1 + M_k, alpha + M_k_sum)
      V_k[i] <- temp[1]
      V_k_minus[i] <- temp[2]
    }
    V_k[N] <- log(1) # equivalent to setting P_N = 1 - sum(P_i)
    # stick-breaking update
    P[1] <- exp(V_k[1])
    for (j in 2:N) {
      # logsum <- sum(log(1  - V_k[1:(j-1)])) + log(V_k[j])
      logsum <- sum(V_k_minus[1:(j-1)]) + V_k[j]
      P[j] <- exp(logsum)
    }
    return(list(P = P, V = V_k_minus))
  }
  
  update.mu0 <- function(mu_k, phi_k, kappa0) {
    # statistic
    a_n <- (a + A * kappa0 * (sum(mu_k * phi_k))) /
      (1 + A * kappa0 * sum(phi_k))
    A_n <- A / (1 + A * kappa0 * sum(phi_k))
    return(rnorm(1, a_n, sqrt(A_n)))
  }
  
  update.kappa0 <- function(mu_k, phi_k, mu0) {
    # statistic
    sh <- 0.5 * (w + N)
    ra <- 0.5 * (W + sum(phi_k * (mu_k - mu0)^2))
    return(rgamma(1, shape = sh,
                  rate = ra))
  }
  
  update.alpha <- function(log_V_minus) {
    return(rgamma(1, shape = c + N - 1,
                  rate = d - sum(log_V_minus[1 : (N-1)])))
  }
  
  # posterior predictives before burnin (untruncated)
  pred.post <- function(mu_k, phi_k, P) {
    # preallocation
    pred <- rep(NA, n.sample)
    # sample component
    component.index <- sample(c(1:N), size = n.sample,
                              replace = TRUE,
                              prob = P)
    for (i in seq(n.sample)) {
      pred[i] <- rlnorm(1, mu_k[component.index[i]],
                        sqrt(1 / phi_k[component.index[i]]))
    }
    return(pred)
  }
  
  # TRUNCATED posterior predictives after burin
  pred.post.after <- function(mu_k, phi_k, P, ns = n.sample) {
    # preallocation
    pred <- rep(NA, ns)
    # sample component
    for (i in seq(ns)) {
      acc <- FALSE
      while(!acc){
        temp_com <- sample(seq(N), size = 1,
                           prob = P)
        temppred <- rlnorm(1, mu_k[temp_com],
                           sqrt(1 / phi_k[temp_com]))
        if (temppred <= max((acq.data))) {
          pred[i] <- temppred 
          acc <- 1
        } else
        {acc <- 0}
      }
    }
    return(pred)
  }
  
  out_len <- floor((niter - burnin1 - burnin2)/thinin)
  # MCMC output preallocation
  # mean.out <- matrix(nrow = niter, ncol = N)
  # var.out <- matrix(nrow = niter, ncol = N)
  K.out <- matrix(nrow = out_len, ncol = length(dat))
  # P.out <- matrix(nrow = niter, ncol = N)
  alpha.out <- rep(NA, out_len)
  # mu0.out <- rep(NA, niter)
  kappa0.out <- rep(NA, out_len)
  spred.out <- rep(NA, out_len)
  
  # initialization
  G <- initial_val()
  
  # MCMC
  for (i in seq(niter)) {
    Z.update <- update.Z(G$K, G$mu0, G$kappa0)
    G$mu_k <- Z.update[,1]
    G$phi_k <- Z.update[,2]
    G$K <- update.K(G$mu_k, G$phi_k, G$P)
    P.V.update <- update.P(G$K, G$alpha)
    G$P <- P.V.update$P
    G$mu0 <- update.mu0(G$mu_k, G$phi_k, G$kappa0)
    G$kappa0 <- update.kappa0(G$mu_k, G$phi_k, G$mu0)
    G$alpha <- update.alpha(P.V.update$V)
    
    # reorder the component by its occupancy
    if(i == burnin1){
      # order the occupied clusters by number of members within clusters
      occ_cluster <- as.numeric(names(summary(as.factor(G$K))))[order(unname(summary(as.factor(G$K))), decreasing = TRUE)]
      un_occ_cluster <- setdiff(seq(N), occ_cluster)
      reindex <- 1
      temp <- G$K
      for (item in occ_cluster) {
        G$K[which(temp == item)] <- reindex
        G$mu_k[reindex] <- G$mu_k[item] 
        G$phi_k[reindex] <- G$phi_k[item]
        reindex = reindex + 1
      }
      l = 0
      for (unuse in un_occ_cluster) {
        G$mu_k[reindex+l] <- G$mu_k[unuse]
        G$phi_k[reindex+l] <- G$phi_k[unuse]
        l = l + 1
      }
    }
    
    if (i > (burnin1 + burnin2) & !i %% thinin) {
      # record all the update
      j = (i-burnin1-burnin2)/thinin
      # record all the update
      # mean.out[i,] <- G$mu_k
      # var.out[i,] <- 1 / G$phi_k
      K.out[j,] <- G$K
      # P.out[i,] <- G$P
      alpha.out[j] <- G$alpha
      # mu0.out[i] <- G$mu0
      kappa0.out[j] <- G$kappa0
      spred.out[j] <- pred.post.after(G$mu_k, G$phi_k, G$P, ns = 1)
    }
  }
  
  return(list(K = K.out,
              alpha = alpha.out,
              kappa = kappa0.out,
              spred = spred.out
              ))
}

#' Function to fit compare predictions of population mean
#' under constrained prior
#' 
#' @param hyperpar: hyperparameter list
#' @param seed.seq: random seed sequence to control randomness 
#' of sampled training sets
#' @param niter: number of MCMC iterations
#' @param n.sample: number of posterior predictions to draw at each MCMC iteration
#' @param n.cluster: truncation level specified to approximate stick breaking process
#' 
#' @return a list of true population mean;
#' predicted mean of constrained model;
#' and predicted mean of KDE.
cDP_pred_mean_compare <- function(hyperpar, seed.seq, niter, n.sample, 
                                  n.cluster, data = acq.data) {
  trueMean <- rep(NA, length(seed.seq))
  cpredMean <- rep(NA, length(seed.seq))
  kdeMean <- rep(NA, length(seed.seq))
  hp <- hyperpar
  for (i in seq(length(seed.seq))) {
    sam.data <- get.logdata(seed.seq[i])
    trueMean[i] <- mean(data[-sam.data$index])
    tr.fit <- cDP_mixture_lnorm_ro(hp, niter, sam.data$data, n.sample, n.cluster)
    kdepred <- pred.kde(sam.data$data, 8983)
    cpredMean[i] <- mean(tr.fit$spred)
    kdeMean[i] <- mean(kdepred)
  }
  return(list(true = trueMean,
              cpred = cpredMean,
              kde = kdeMean))
}

#' Function to get fitting results of constrained model
#' 
#' @param hyperpar: hyperparameter list
#' @param seed.seq: random seed sequence to control randomness 
#' of sampled training sets
#' @param niter: number of MCMC iterations
#' @param n.sample: number of posterior predictions to draw at each MCMC iteration
#' @param n.cluster: truncation level specified to approximate stick breaking process
#' 
#' 
#' @return a list consisting of fitting results under different 
#' training sample and values for the holdout dataset.
cDP_get_fit <- function(hyperpar, seed.seq, niter, n.sample, 
                        n.cluster, data = acq.data) {
  fit <- list()
  popdata <- list()
  for (i in seq(length(seed.seq))) {
    dat <- get.logdata(seed.seq[i])
    fit[[i]] <- cDP_mixture_lnorm_ro(hyperpar, niter, dat$data, n.sample, 
                                     n.cluster)
    popdata[[i]] <- data[-dat$index]
  }
  return(list(fit = fit,
              popdata = popdata))
}


#' Function to fit NP mixture model with flexible prior
#' 
#' @param hyperpar: hyperparameter list
#' @param niter: number of total MCMC iterations
#' @param dat: log-scaled training dataset
#' @param n.sample: number of posterior predictions to draw at each MCMC iteration
#' @param n.cluster: truncation level specified to approximate stick breaking process
#' @param burnin1: burnin of MCMC; number of steps for shufflings; default as 3000
#' @param burnin2: burnin after shufflings
#' @param thinin: thinin of MCMC
#' 
#' @return a list of all posterior draws and posterior prediction matrix
hDP_mixture_lnorm_ro <- function(hyperpar, niter, dat, n.sample, 
                                 n.cluster, burnin1 = 1000, burnin2 = 1000,
                                 thinin = 1) {
  
  N = n.cluster
  
  # hyperparameter
  g1 <- hyperpar$g1
  g2 <- hyperpar$g2
  v1 <- hyperpar$v1
  ss0 <- hyperpar$ss0
  v0 <- hyperpar$v0
  df <- hyperpar$df
  c <- hyperpar$c
  d <- hyperpar$d
  a <- hyperpar$a
  A <- hyperpar$A
  

  initial_val <- function() {
    km <- kmeans(dat, N)
    # sort the mean, variance and mixing weights in decreasing order
    km_ordered <- sort(km$centers)
    km_cluster_ordered <- rep(NA, length(km$cluster))
    km_phi_ordered <- rep(NA, length(km$withinss))
    km_size_ordered <- rep(NA, length(km$size))
    for (i in seq(length(km_ordered))) {
      index <- which(km$centers == km_ordered[i])
      km_size_ordered[i] <- km$size[index]
      km_phi_ordered[i] <- 1 / (km$withinss[index]/km$size[index]) 
      index2 <- which(km$cluster == index)
      km_cluster_ordered[index2] <- i
    }
    km_p_ordered <- km_size_ordered / length(dat)
    return(list(P = km_p_ordered,
                K = km_cluster_ordered,
                mu_k = km_ordered,
                phi_k = km_phi_ordered,
                h_k = rgamma(N, v1, rate = v1),
                alpha = rgamma(1, c, rate = d),
                mu0 = mean(dat),
                r_k = rgamma(N, df/2, rate = df/2),
                kappa_0 = 1))
  }
  
  # Marsaglia and Tsang's method for sampling log(beta)
  MT_loggam <- function(alpha){
    if (alpha >= 1){
      ua <- alpha
      add = FALSE
    } else {
      ua <- alpha + 1
      add = TRUE
    }
    flag <- 1
    while (flag) {
      d <- ua - 1/3
      c <- 1 / sqrt(9*d)
      Z <- rnorm(1)
      if (Z > -1 / c){
        V <- (1 + c*Z)^3
        U <- runif(1)
        flag = (log(U) > (0.5*Z^2+d-d*V+d*log(V)))
      }
    }
    if(add) {
      return(log(d) +log(V) +(1 / (alpha))*log(runif(1)))
    } else{
      return(log(d)+ log(V))
    }
  }
  
  # Marsaglia and Tsang's method for sampling log(beta)
  lbeta_from_lgam <- function(alpha, beta){
    A <- MT_loggam(alpha)
    B <- MT_loggam(beta)
    A2 <- A - max(A, B)
    B2 <- B - max(A, B)
    A3 <- A2 - log1p(exp(min(A2, B2)))
    B3 <- B2 - log1p(exp(min(A2, B2)))
    return(c(A3, B3))
  }
  
  kappa.update <- function(mu0, Mu, Phi, r) {
    return(rgamma(1, 0.5*N + g1,
                  rate = 0.5*(sum(Phi*r*(Mu - mu0)^2)) + g2))
  }
  
  Z.update <- function(K, mu0, Mu, kappa0, h, r) {
    unique.K <- unique(K)
    
    # preallocation
    mu_k <- rep(NA, N)
    phi_k <- rep(NA, N)
    h_k <- rep(NA, N)
    r_k <- rep(NA, N)
    Z_matrix <- matrix(NA, nrow = N, ncol = 4)
    
    for (i in seq(N)) {
      if(sum(i == unique.K)){
        # some statistics
        i.index <- which(i == K)
        n_j.star <- length(i.index)
        Apost <- v0 + 0.5*(n_j.star + 1)
        Bpost <- ss0*h[i] + 0.5*(kappa0*r[i]*(Mu[i]-mu0)^2) +
          0.5*(sum((dat[i.index]-Mu[i])^2))
        # update
        phi_k[i] <- rgamma(1, Apost,
                           rate = Bpost)
        h_k[i] <- rgamma(1, v0+v1, rate = v1+ss0*phi_k[i])
        Mpost <- (phi_k[i]*n_j.star*mean(dat[i.index]) + kappa0*phi_k[i]*r[i]*mu0)/
          (phi_k[i]*n_j.star + kappa0*phi_k[i]*r[i])
        Vpost <- 1/(phi_k[i]*n_j.star + kappa0*phi_k[i]*r[i])
        mu_k[i] <- rnorm(1, Mpost, sqrt(Vpost))
        r_k[i] <- rgamma(1, 0.5*(df + 1), rate = 0.5*(df+kappa0*phi_k[i]*(Mu[i]-mu0)^2))
      } else{
        h_k[i] <- rgamma(1, v1, rate = v1)
        phi_k[i] <- rgamma(1, v0, rate = ss0*h_k[i])
        r_k[i] <- rgamma(1, df/2, rate = df/2)
        mu_k[i] <- rnorm(1, mu0, sqrt(1/phi_k[i]/r_k[i]/kappa0))
      }
    }
    Z_matrix[,1] <- mu_k
    Z_matrix[,2] <- phi_k
    Z_matrix[,3] <- h_k
    Z_matrix[,4] <- r_k
    return(Z_matrix)
  }
  
  K.update <- function(Mu, Phi, P) {
    # preallocation
    K <- rep(NA, length(dat))
    for (i in seq(dat)) {
      prob <- P * dnorm(dat[i], Mu, sqrt(1 / Phi))
      K[i] <- sample(seq(N), size = 1, prob = prob)
    }
    return(K)
  }
  
  P.update <- function(K, alpha) {
    # preallocation
    V_k <- rep(NA, N)
    V_k_minus <- rep(NA, N-1) # log(rbeta)
    P <- rep(NA, N)
    # calculate M_k's
    M_vector <- sapply(c(1 : N), function(i, x) 
    {sum(i == x)}, x = K)
    for (i in seq(N-1)) {
      # statistic
      M_k <- M_vector[i]
      M_k_sum <- sum(M_vector[(i + 1): N])
      temp <- lbeta_from_lgam(1 + M_k, alpha + M_k_sum)
      # V_k[i] <- rbeta(1, 1 + M_k, alpha + M_k_sum)
      V_k[i] <- temp[1]
      V_k_minus[i] <- temp[2]
    }
    V_k[N] <- log(1) # equivalent to setting P_N = 1 - sum(P_i)
    # stick-breaking update
    P[1] <- exp(V_k[1])
    for (j in 2:N) {
      # logsum <- sum(log(1  - V_k[1:(j-1)])) + log(V_k[j])
      logsum <- sum(V_k_minus[1:(j-1)]) + V_k[j]
      P[j] <- exp(logsum)
    }
    return(list(P = P, V = V_k_minus))
  }
  
  mu0.update <- function(Mu, Phi, Kappa, r) {
    # some statistics
    mpost <- (A*Kappa*sum(Mu*Phi*r) + a) / (1 + A*Kappa*sum(Phi*r))
    vpost <- A / (1 + A*Kappa*sum(Phi*r))
    return(rnorm(1, mpost, sqrt(vpost)))
  }
  
  alpha.update <- function(log_V_minus) {
    return(rgamma(1, shape = c + N - 1,
                  rate = d - sum(log_V_minus[1 : (N-1)])))
  }
  
  # posterior predictives before burnin
  pred.post <- function(Mu, Phi, P) {
    # preallocation
    pred <- rep(NA, n.sample)
    # sample component
    component.index <- sample(c(1:N), size = n.sample,
                              replace = TRUE,
                              prob = P)
    for (i in seq(n.sample)) {
      pred[i] <- rlnorm(1, Mu[component.index[i]],
                        sqrt(1 / Phi[component.index[i]]))
    }
    return(pred)
  }
  
  # TRUNCATED posterior predictives after burin
  pred.post.after <- function(mu_k, phi_k, P, ns = n.sample) {
    # preallocation
    pred <- rep(NA, ns)
    # sample component
    for (i in seq(ns)) {
      acc <- FALSE
      while(!acc){
        temp_com <- sample(seq(N), size = 1,
                           prob = P)
        temppred <- rlnorm(1, mu_k[temp_com],
                           sqrt(1 / phi_k[temp_com]))
        if (temppred <= max((acq.data))) {
          pred[i] <- temppred 
          acc <- 1
        } else
        {acc <- 0}
      }
    }
    return(pred)
  }
  
  
  # MCMC output preallocation
  # mu.out <- matrix(nrow = niter, ncol = N)
  # phi.out <- matrix(nrow = niter, ncol = N)
  out_len <- floor((niter - burnin1 - burnin2)/thinin)
  K.out <- matrix(nrow = out_len, ncol = length(dat))
  # P.out <- matrix(nrow = niter, ncol = N)
  # h.out <- matrix(nrow = niter, ncol = N)
  alpha.out <- rep(NA, out_len)
  # mu0.out <- rep(NA, niter)
  # r.out <- matrix(nrow = niter, ncol = N)
  kappa.out <- rep(NA, out_len)
  spred.out <- rep(NA, out_len)
  
  # initialization
  G <- initial_val()
  
  for (i in seq(niter)) {
    G$kappa_0 <- kappa.update(G$mu0, G$mu_k, G$phi_k, G$r_k) 
    update.Z <- Z.update(G$K, G$mu0, G$mu_k, G$kappa_0, G$h_k, G$r_k)
    G$mu_k <- update.Z[,1]
    G$phi_k <- update.Z[,2]
    G$h_k <- update.Z[,3] 
    G$r_k <- update.Z[,4]
    G$K <- K.update(G$mu_k, G$phi_k, G$P)
    P.V.update <- P.update(G$K, G$alpha) 
    G$P <- P.V.update$P 
    G$mu0 <- mu0.update(G$mu_k, G$phi_k, G$kappa_0, G$r_k) 
    G$alpha <- alpha.update(P.V.update$V)
    
    # reorder the component by its occupancy
    if(i == burnin1){
      # order the occupied clusters by number of members within clusters
      occ_cluster <- as.numeric(names(summary(as.factor(G$K))))[order(unname(summary(as.factor(G$K))), decreasing = TRUE)]
      un_occ_cluster <- setdiff(seq(N), occ_cluster)
      reindex <- 1
      temp <- G$K
      for (item in occ_cluster) {
        G$K[which(temp == item)] <- reindex
        G$mu_k[reindex] <- G$mu_k[item] 
        G$phi_k[reindex] <- G$phi_k[item]
        G$r_k[reindex] <- G$r_k[item]
        G$h_k[reindex] <- G$h_k[item]
        reindex = reindex + 1
      }
      l = 0
      for (unuse in un_occ_cluster) {
        G$mu_k[reindex+l] <- G$mu_k[unuse]
        G$phi_k[reindex+l] <- G$phi_k[unuse]
        G$r_k[reindex+l] <- G$r_k[unuse]
        G$h_k[reindex+l] <- G$h_k[unuse]
        l = l + 1
      }
    }
    
    
    if (i > (burnin1 + burnin2) & !i %% thinin) {
      # record all the update
      j = (i-burnin1-burnin2)/thinin
      # mu.out[i,] <- G$mu_k
      # phi.out[i,] <- G$phi_k
      K.out[j,] <- G$K
      # P.out[i,] <- G$P
      # h.out[i,] <- G$h_k
      alpha.out[j] <- G$alpha
      # mu0.out[i] <- G$mu0
      kappa.out[j] <- G$kappa_0
      # r.out[i,] <- G$r_k
      spred.out[j] <- pred.post.after(G$mu_k, G$phi_k, G$P, ns = 1)
    }
  }
  return(list(K = K.out,
              alpha = alpha.out,
              kappa = kappa.out,
              spred = spred.out))
  
}

#' Function to fit compare predictions of population mean
#' under flexible prior
#' 
#' @param hyperpar: hyperparameter list
#' @param seed.seq: random seed sequence to control randomness 
#' of sampled training sets
#' @param niter: number of total MCMC iterations
#' @param n.sample: number of posterior predictions to draw at each MCMC iteration
#' @param n.cluster: truncation level specified to approximate stick breaking process
#' 
#' @return a list of true population mean;
#' predicted mean of flexible model;
#' and predicted mean of KDE.
hDP_pred_mean_compare <- function(hyperpar, seed.seq, niter, n.sample, 
                                  n.cluster, data = acq.data) {
  trueMean <- rep(NA, length(seed.seq))
  cpredMean <- rep(NA, length(seed.seq))
  kdeMean <- rep(NA, length(seed.seq))
  hp <- hyperpar
  for (i in seq(length(seed.seq))) {
    sam.data <- get.logdata(seed.seq[i])
    trueMean[i] <- mean(data[-sam.data$index])
    tr.fit <- hDP_mixture_lnorm_ro(hp, niter, sam.data$data, n.sample, n.cluster)
    kdepred <- pred.kde(sam.data$data, 8983)
    cpredMean[i] <- mean(tr.fit$spred)
    kdeMean[i] <- mean(kdepred)
  }
  return(list(true = trueMean,
              cpred = cpredMean,
              kde = kdeMean))
}

#' Function to get fitting results of flexible model
#' 
#' @param hyperpar: hyperparameter list
#' @param seed.seq: random seed sequence to control randomness 
#' of sampled training sets
#' @param niter: number of MCMC iterations
#' @param n.sample: number of posterior predictions to draw at each MCMC iteration
#' @param n.cluster: truncation level specified to approximate stick breaking process
#' 
#' @return a list consisting of fitting results under different 
#' training sample and values for the holdout dataset.
hDP_get_fit <- function(hyperpar, seed.seq, niter, n.sample, 
                        n.cluster, data = acq.data) {
  fit <- list()
  popdata <- list()
  for (i in seq(length(seed.seq))) {
    dat <- get.logdata(seed.seq[i])
    fit[[i]] <- hDP_mixture_lnorm_ro(hyperpar, niter, dat$data, n.sample, 
                                     n.cluster)
    popdata[[i]] <- data[-dat$index]
  }
  return(list(fit = fit,
              popdata = popdata))
}

# list(v0 = 10, ss0 = 3, a = 0, A = 100, g1 = 2, g2 = 8, c = 8, d = 1)
# list(v0 = 10, ss0 = 3, a = 0, A = 100, g1 = 2, g2 = 8, c = 8,d=1, df = 8, v1 = 5)
