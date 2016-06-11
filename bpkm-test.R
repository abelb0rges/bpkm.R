# benchmarking with stats::kmeans
# iris dataset
dataset <- as.matrix(iris[,1:4])
test <- kmeans(dataset, 3)
failrate(labs=iris$Species, sol=test$cluster)
sols <- list()
length(sols) <- 10
for(i in 1:10) {
  sol <- bpkm(X = dataset, K = 3, min_it = 4)$sol
  sols[[i]] <- sol
}
fails <- lapply(sols, function(ll) {
  failrate(iris$Species, ll)
}) %>% as.numeric
meanfail <- mean(fails)

# simulation study
bpkm_sim <- function(nreps=100, MAHALANOBIS, RBASED) {
  require(MASS)
  mu1 <- c(0,0)
  mu2 <- c(2,2)
  S1 <- matrix(c(1, .5,
                 .5, 1), 2, byrow=TRUE)
  S2 <- matrix(c(1.5, -.2,
                 -.2, 1.5), 2, byrow=TRUE)
  N1 <- 40
  N2 <- 60
  
  # list of solutions
  sols <- list()
  length(sols) <- nreps
  for(i in 1:nreps) {
    X1 <- mvrnorm(N1, mu1, S1)
    X2 <- mvrnorm(N2, mu2, S2)
    X <- cbind(rbind(X1, X2), c(rep(1, N1), rep(2, N2)))
    
    sols[[i]] <- bpkm(X=X, K=2, min_it=2,
                                MAHALANOBIS=MAHALANOBIS,
                                RBASED=RBASED)
  }
  
  # fail rate
  frate <- lapply(sols, function(ll) {
    labs <- ll$X[,3]
    sol <- ll$sol
    failrate(labs, sol)
  }) %>% as.numeric
  
  frate
}

# comparison according to:
# 1. use mahalanobis or euclidian distance
# 2. use r matrix-based or centers-based first guess
bpkm_compare <- function() {
  # returns the mean fail rate for each combination
  # of distance function and first guess approach
  frate <- matrix(0, 2, 2)
  rownames(frate) <- c('mahalanobis', 'euclides')
  colnames(frate) <- c('r_based','m_based')
  i <- 0
  for(m in c(TRUE, FALSE)) {
    i <- i + 1
    j <- 0
    for(rb in c(TRUE, FALSE)) {
      j <- j + 1
      frate[i,j] <- bpkm_sim(nreps=100, MAHALANOBIS=m,
                                       RBASED=rb) %>% mean
    }
  }
  
  frate
}
