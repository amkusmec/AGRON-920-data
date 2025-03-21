require(RSpectra)

### N.B.: Optimization of the Gaussian and deep kernels should be done using 
### the *training* phenotypic data to reduce information leakage.

### Marginal likelihood for the bandwidth parameter of a Gaussian kernel following
### Perez-Elizalde et al. (2015) J. Agric. Biol. Environ. Stat.
# https://github.com/gcostaneto/KernelMethods/blob/master/Gaussian_Kernel.R
margh.fun <- function(theta, y, D, q = 1, nu = 0.0001, Sc = 0.0001, nuh = NULL, 
                      Sch = NULL, prior = NULL) {
  h <- theta[1]
  phi <- theta[2]
  
  Kh <- exp(-h*D/q)
  
  eigenKh <- eigs_sym(Kh, k = nrow(Kh) - 1, which = "LM")
  nr <- length(which(eigenKh$values > 1e-10))
  d <- t(eigenKh$vectors[, 1:nr]) %*% matrix(y, ncol = 1)
  
  Iden <- -1/2*sum(log(1 + phi*eigenKh$values[1:nr])) - 
    (nu + nr - 1)/2*log(Sc + sum(d^2/(1 + phi*eigenKh$values[1:nr])))
  if(!is.null(prior)) lprior <- dgamma(h, nuh, Sch, log = TRUE) else lprior <- 0
  
  Iden <- -(Iden + lprior)
  
  return(Iden)
}

### Marginal likelihood for the first-order arc-cosine kernel at given number of 
### layers following Cuevas et al. (2019) G3
# https://github.com/gcostaneto/KernelMethods/blob/master/DeepKernels.R
margl.fun <- function(y, Kh) {
  lden.fun <- function(theta, nr, Sh, d) {
    phi <- theta[1]
    lden  <- -1/2*sum(log((1 + phi*Sh))) - (nr - 1)/2*log(sum(d^2/((1 + phi*Sh))))
    lden <- -(lden)
    return(lden)
  }
  
  eigenKh <- eigs_sym(Kh, k = nrow(Kh) - 1, which = "LM")
  nr <- length(which(eigenKh$values > 1e-10))
  d <- t(eigenKh$vectors[, 1:nr]) %*% matrix(y, ncol = 1)
  
  sol <- optim(c(1), lden.fun, nr = nr, Sh = eigenKh$values[1:nr], d = d, 
               method = "L-BFGS-B", lower = c(0.0005), upper = c(200))
  phi <- sol$par[1]
  
  log.vero <- -1/2*sum(log((1 + phi*eigenKh$values[1:nr]))) - 
    (nr - 1)/2*log(sum(d^2/((1 + phi*eigenKh$values[1:nr]))))
  
  return(log.vero)
}

### Normalized by the average of the diagonal
linear_kernel <- function(X, normalize = TRUE) {
  nr <- nrow(X)
  K <- tcrossprod(X)
  
  if (normalize) {
    K <- K/(sum(diag(K))/nr)
  }
  
  rownames(K) <- colnames(K) <- rownames(X)
  return(K)
}

gaussian_kernel <- function(X, h = 1, prob = NULL, y = NULL, verbose = FALSE) {
  # Preserves row names
  if (verbose) cat("Distance matrix...\n")
  D <- as.matrix(dist(X, method = "euclidean", diag = TRUE, upper = TRUE))^2
  
  # Scaling factor for numerical stability
  if (is.null(prob)) {
    q <- 1
  } else {
    q <- quantile(D, probs = prob)
  }
  
  # Optimize the bandwidth parameter?
  if (!is.null(y)) {
    if (verbose) cat("Optimizing bandwidth...\n")
    
    sol <- optim(c(1, 3), margh.fun, y = y, D = D, q = q, method = "L-BFGS-B", 
                 lower = c(0.01, 0.05), upper = c(6, 30), prior = "gamma", nuh = 3, 
                 Sch = 1.5, control = list(trace = 4L))
    h <- sol$par[1]
    if(verbose) cat("Optimal h =", h, "; fn evals =", sol$counts[1], "\n")
  }
  
  return(list(Kernel = exp(-h*D/q), 
              h = h, q = q))
}

# This version normalizes the kernel by its median entry at each level. This is 
# actually how Cuevas et al. (2019) compute the kernel. It no longer has the 
# theoretical properties that the un-normalized version has. In practice, 
# it seems to work *much* better for predicting phenotypes.
deep_kernel <- function(X, nl = 1, y = NULL, verbose = FALSE) {
  ### Single-layer kernel
  # (Constant) magnitude dependence
  LK <- tcrossprod(X)
  norms <- outer(sqrt(diag(LK)), sqrt(diag(LK)))
  
  # Angular dependence
  c0 <- LK/norms
  c0[c0 > 1] <- 1   # Rounding errors can lead to values >1
  theta <- acos(c0)
  DK1 <- (1/pi)*norms*(sin(theta) + (pi - theta)*cos(theta))
  mm <- median(DK1)
  DK1 <- DK1/mm
  
  # Clean-up
  rm(LK, c0, norms); gc()
  
  # Layer optimization interprets `nl` as a maximum number of layers
  if (!is.null(y)) {
    l <- 1
    ll0 <- margl.fun(y, DK1)
    if (verbose) cat("l =", l, "; LL =", ll0, "\n")
    
    DK2 <- DK1
    repeat {
      l <- l + 1
      if (l > nl) break   # Reached the maximum number of layers
      
      # Compute the marginal likelihood for the next layer
      norms <- outer(sqrt(diag(DK1)), sqrt(diag(DK1)))
      c0 <- DK1/norms
      c0[c0 > 1] <- 1
      theta <- acos(c0)
      DK2 <- (1/pi)*norms*(sin(theta) + (pi - theta)*cos(theta))
      mm <- c(mm, median(DK2))
      DK2 <- DK2/mm[l]
      ll1 <- margl.fun(y, DK2)
      if (verbose) cat("l =", l, "; LL =", ll1, "\n")
      
      if (ll1 <= ll0) break   # No increase in the marginal likelihood
      
      ll0 <- ll1
      DK1 <- DK2
    }
    if(verbose) cat("Optimal l =", l - 1, "\n")
    nl <- l - 1
  } else if (nl > 1) {
    for (i in 2:nl) {
      norms <- outer(sqrt(diag(DK1)), sqrt(diag(DK1)))
      c0 <- DK1/norms
      c0[c0 > 1] <- 1
      theta <- acos(c0)
      DK2 <- (1/pi)*norms*(sin(theta) + (pi - theta)*cos(theta))
      mm <- c(mm, median(DK2))
      DK1 <- DK2/mm[i]
    }
  }
  
  # Final kernel
  return(list("Kernel" = DK1, 
              "Medians" = mm[1:nl], 
              "l" = nl))
}
