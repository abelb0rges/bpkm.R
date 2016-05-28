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
  crossprod(x-y) %>% sqrt
}

# weighted euclidian distance
wed <- function(x, y=rep(0,length(x)), w) {
  if(length(x) != length(y))
    stop('x and y must have same length')
  res <- 0
  for(i in seq_along(x))
    res <- res + w[i] * (x[i] - y[i])^2
  
  sqrt(res)
}

# euclidian distance-based cost function
cost_ed <- function(r, mu, X) {
  require(magrittr)
  
  mu_sum <- (apply(mu, 1, crossprod) %*% apply(r, 2, sum)) %>%
    as.numeric
  cross_sum <- (X %*% t(mu) %*% t(r)) %>% sum
  
  mu_sum - 2 * cross_sum
}

# mahalanobis distance-based cost function
cost_md <- function(r, mu, X) {
  N <- dim(r)[1]; R <- dim(r)[2]
  K <- dim(mu)[1]
  
  J <- 0
  for(k in 1:K) {
    kth <- (r[,k] == 1) # rows of k-th cluster
    S <- cov(X[kth,])
    # TO-DO: include this check on S
    # if(S is not singular) {
    for(n in 1:N) {
      J <- J + mahalanobis(X[n,], mu[k,], S)
    }
    #}
  }
  
  J
}

bpkm <- function(X, K, init_iter=10, min_it=2, MAHALANOBIS=FALSE, RBASED=TRUE) {
  # X: (N,p)-matrix with n R^p points
  # K: number of groups
  
  # dependencies
  if(!require(magrittr))
    stop("run install.packages('magrittr') please")
  if(!require(lpSolve))
    stop("run install.packages('lpSolve') please")
  
  if(MAHALANOBIS)
    COST <- cost_md
  else
    COST <- cost_ed
  
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
        best <- apply(centers, 1, function(m) ed(xx, m)) %>%
          which.min
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
      centers <- (t(r) %*% X) %>% apply(2, function(m) m/nk)
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
    if(RBASED)
      init_part <- rbased_guess()
    else
      init_part <- guess()
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
          ind %<>% { . + 1 }
          if(MAHALANOBIS) {
            kth <- (oldr[k,] == 1)
            S <- cov(X[kth,])
            obj[ind] <- mahalanobis(X[n,],mu[k,],S)
          } else {
            obj[ind] <- (X[n,] - mu[k,]) %>% as.numeric %>% ed
          }
        }
      }
      r <- lp(objective.in = obj, const.mat = A, const.dir = dir,
              const.rhs = rhs, all.bin = TRUE) %$%
        solution %>% matrix(ncol=K, byrow=TRUE)
      # step 2: mu optimisation (update cluster's mean)
      mu %<>% (function(m) {
        for(k in 1:K) {
          kth <- (r[,k] == 1) # rows of k-th cluster
          m[k,] <- as.matrix(X[kth,], nrow=length(kth)) %>%
            apply(2,mean)
        }
        return(m)
      })
      
      bad <- any(!(oldr == r))
      oldr <- r
    }
  }
  
  sol <- numeric(N)
  for(i in 1:N)
    sol[i] <- which(r[i,] == 1)
  
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

