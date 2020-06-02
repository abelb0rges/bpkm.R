# --------------------------------------
# binary lp-based k-means implementation
# author: Abel Borges
# --------------------------------------

# r: matrix N-by-K classifying columns of x
# mu: matrix K-by-p of centroids of clusters
# X: matrix N-by-p with vectors to be classified

# euclidian distance
ed <- function(x, y=rep(0,length(x))) {
  if(length(x) != length(y))
    stop('x and y must have same length')
  sqrt(crossprod(x-y))
}

# weighted euclidian distance
wed <- function(x, y=rep(0,length(x)), w) {
  if(length(x) != length(y))
    stop('x and y must have same length')
  sqrt(w * (x - y)^2)
}

# euclidian distance-based cost function
cost_ed <- function(r, mu, X) {
  mu_sum <- as.numeric(apply(mu, 1, crossprod) %*% apply(r, 2, sum))
  cross_sum <- sum(X %*% t(mu) %*% t(r))
  mu_sum - 2 * cross_sum
}

# mahalanobis distance-based cost function
cost_md <- function(r, mu, X) {
  require(MASS)
  N <- dim(r)[1]
  R <- dim(r)[2]
  K <- dim(mu)[1]
  J <- 0
  
  for(k in 1:K) {
    S <- cov(X[r[,k] == 1,])
    Sinv <- ginv(S)
    for(n in 1:N) {
      J <- J + mahalanobis(X[n,], mu[k,], Sinv, inverted=TRUE)
    }
  }
  
  J
}

bpkm <- function(X, K, init_iter=10, min_it=2, MAHALANOBIS=FALSE, RBASED=TRUE) {
  # X: (N,p)-matrix with n R^p points
  # K: number of groups
  
  # dependencies
  require(MASS)
  if(!require(lpSolve)) stop("run install.packages('lpSolve') please")
  
  COST <- if(MAHALANOBIS) cost_md else cost_ed
  
  N <- NROW(X)
  p <- NCOL(X)
  # initial partition
  guess <- function() {
    r <- matrix(0, N, K)
    init_part <- list(r = r, centers = NULL)
    cost <- Inf
    
    for(i in 1:init_iter) {
      centers <- X[sample(1:N,K),]
      r <- matrix(0, N, K)
      for(j in 1:N) {
        xx <- as.numeric(X[j,])
        best <- which.min(apply(centers, 1, function(m) ed(xx, m)))
        r[j, best] = 1
      }
                                
      newcost <- COST(r, centers, X)
      if(newcost < cost) {
        init_part$r <- r
        init_part$centers <- centers
        cost <- newcost
      }
    }
    
    init_part
  }
  
  # initial partition based on randomisation of r directly
  rbased_guess <- function() {
    r <- matrix(0, N, K)
    init_part <- list(r = r, centers = NULL)
    cost <- Inf
    
    for(i in 1:init_iter) {
      for(j in 1:N) r[j, sample(1:K, 1)] <- 1
      nk <- apply(r, 2, sum)
      centers <- apply(t(r) %*% X, 2, function(m) m/nk)
      newcost <- COST(r, centers, X)
      if(newcost < cost) {
        init_part$r <- r
        init_part$centers <- centers
        cost <- newcost
      }
    }

    init_part
  }
  
  # ----------
  # K-means!!!
  # ----------
  
  # for binary optimisation with lpSolve:
  
  # matrix of constraints
  A <- matrix(0, N, N*K)
  for(n in 1:N) A[n, ((n-1)*K + 1):(n*K)] <- 1
  for(k in 1:K) {
    kth <- numeric(K)
    kth[k] <- 1
    kth <- rep(kth, N)
    A <- rbind(A, kth)
  }
  # directions of constraints
  dir <- c(rep("==", N), rep(">=", K))
  # rhs of constraints
  rhs <- c(rep(1, N), rep(1, K))
  
  niter <- 0
  while(niter <= min_it) {
    # initial values for k-means iterations
    init_part <- if(RBASED) rbased_guess() else guess()
    oldr <- init_part$r
    mu <- init_part$centers
    
    # stop criteria
    niter <- 0
    bad <- TRUE
    while(bad) {
      niter <- niter + 1
      cat(niter,'...\n')
      
      # step 1: r optimisation
      obj <- numeric(N*K)
      ind <- 0
      for(n in 1:N) {
        for(k in 1:K) {
          ind <- ind + 1
          if(MAHALANOBIS) {
            S <- cov(X[oldr[k,] == 1,])
            Sinv <- ginv(S)
            obj[ind] <- mahalanobis(X[n,], mu[k,], Sinv, inverted = TRUE)
          } else {
            obj[ind] <- ed(as.numeric((X[n,] - mu[k,])))
          }
        }
      }
      lp_res <- lp(objective.in = obj, const.mat = A, const.dir = dir,
              const.rhs = rhs, all.bin = TRUE)
      r <- matrix(lp_res$solution, ncol = K, byrow = TRUE)
      # step 2: mu optimisation (update cluster's mean)
      for(k in 1:K) {
        kth <- (r[,k] == 1) # rows of k-th cluster
        mu[k,] <- apply(as.matrix(X[kth,], nrow = length(kth)), 2, mean)
      }
      
      bad <- any(!(oldr == r))
      oldr <- r
    }
  }
  
  sol <- numeric(N)
  for(i in 1:N) sol[i] <- which(r[i,] == 1)
  
  list(X=X, centers=mu, sol=sol, rsol=r, niter=niter)
}

failrate <- function(labs, sol) {
  # labs: labels in same order of the
  #+database rows entry for reach solution 'sol'
  
  N <- length(labs)
  errors <- logical(N*(N-1)/2)
  ind <- 0
  for(i in 1:(N-1)) {
    for(j in (i+1):N) {
      ind <- ind + 1
      flag1 <- (labs[i] == labs[j])
      flag2 <- (sol[i] == sol[j])
      errors[ind] <- (flag1 != flag2)
    }
  }
  
  mean(errors)
}
