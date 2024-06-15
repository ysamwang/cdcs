n <- 100
ss <- 1000
rec <- matrix(0, nrow = ss, ncol = 4)
for(i in 1:ss){
  # X <- rmutil::rlaplace(n, 0, 1)
  # Y <- 1.5 * X + rmutil::rlaplace(n, 0, 1)
  
  # X <- runif(n, -sqrt(3), sqrt(3))
  # Y <- -.8 * X + runif(n, -sqrt(3), sqrt(3))
  
  X <- rgamma(n, 1, 1) - 1
  Y <- .2 * X + rgamma(n, 1, 1) - 1
  
  X <- scale(X)
  Y <- scale(Y)
  
  G <- array(0, dim = c(n, 10, 1))
  
  
  G[, 1, 1] <- scale(sin(Y))
  G[, 2, 1] <- scale(cos(Y))
  temp <- 2
  G[, 3, 1] <- scale(sin(temp * Y))
  G[, 4, 1] <- scale(cos(temp * Y))
  temp <- 3
  G[, 5, 1] <- scale(sin(temp * Y))
  G[, 6, 1] <- scale(cos(temp * Y))
  G[, 7, 1] <- Y^2
  G[, 8, 1] <- Y^3
  G[, 9, 1] <- scale(ifelse(abs(Y) < 1, Y^2, Y^3))
  G[, 10, 1] <- scale(ifelse(abs(Y) < 1, Y^3, Y^2))

  out <- cdcs::gofTest(covariates = matrix(Y, ncol = 1), Y = X, G =G[,c(1,2, 7, 8) , ,drop = F], bs = 200,
                withinAgg = 2, intercept = 1, sampleSplit = 0)
  
  rec[i, 1] <- mean(out$testStatTwo < out$nullDistTwo )
  
  out <- cdcs::gofTest(covariates = matrix(Y, ncol = 1), Y = X, G =G[, 1:6 , , drop = F], bs = 200,
                       withinAgg = 2, intercept = 1, sampleSplit = 0)
  
  rec[i, 2] <- mean(out$testStatTwo < out$nullDistTwo )
  
  out <- cdcs::gofTest(covariates = matrix(Y, ncol = 1), Y = X, G =G[, 7:8 , , drop = F], bs = 200,
                       withinAgg = 2, intercept = 1, sampleSplit = 0)
  
  rec[i, 3] <- mean(out$testStatTwo < out$nullDistTwo )
  
  out <- cdcs::gofTest(covariates = matrix(Y, ncol = 1), Y = X, G =G[, 9:10 , , drop = F], bs = 200,
                       withinAgg = 2, intercept = 1, sampleSplit = 0)
  
  rec[i, 4] <- mean(out$testStatTwo < out$nullDistTwo )
  
  if(i %%100 == 0){
    cat(i)
  }
}

colMeans(rec < .1)

###############

p <- 5
n <- 500
distro <- "gamma"
ss <- 1000
K <- 10
rec <- matrix(0, ss, 2)
pa <- 2
ch <- 1
M <- 5
bs <- 200
for(i in 1:ss){
  dat <- cdcs::rDAG(p, n, parent_prob = 1/3, lowScale = .8,
                    highScale = 1, edgeVar = n^(-1/8),
                    dist = distro, uniqueTop = T)
  
  Y <- scale(dat$Y)
  
  G <- array(0, dim = c(n, 2 * K, p))
  
  for(j in 1:p){
    for(k in 1:K){
      temp <- rnorm(1)
      G[, 2*(k-1) + 1, j] <- scale(sin(temp* Y[, j]))
      G[, 2*(k), j] <- scale(cos(temp * Y[,j]))
    }
  }

  
  rec[i,1] <- bnbHelperanm(poly(Y[,pa, drop = F], M), Y[,ch, drop = F],
                         G[, , pa, drop = F], withinAgg = 2, aggType = 2, bs = bs,
                         intercept = 1)$pVals[1]
  out <- gofTest(Y[,pa, drop = F], Y[,ch, drop = F],
                 G[, , pa, drop = F], withinAgg = 2, bs = bs,
                 intercept = 1, sampleSplit = F)
  
  rec[i,2] <- (sum(out$nullDistTwo > out$testStatTwo) + 1) / (bs + 1)

  if(i %% 100 == 0 ){
    cat(i)
  }
}
par(mfrow = c(2,1))
hist(rec[,1])
hist(rec[,2])
colMeans(rec < .1)

p <- 8
n <- 5000
K <- 1
dat <- cdcs::rDAG(p, n, parent_prob = 1/3, lowScale = .8,
                  highScale = 1, edgeVar = .5,
                  dist = distro, uniqueTop = T)

Y <- scale(dat$Y)

G <- array(0, dim = c(n, 2 * K, p))

for(j in 1:p){
  for(k in 1:K){
    temp <- rnorm(1)
    # G[, 2*(k-1) + 1, j] <- sin(temp* Y[, j])
    # G[, 2*(k), j] <- cos(temp * Y[,j])
    G[, 2*(k-1) + 1, j] <- sqrt(abs(Y[,j])) * sign(Y[,j])
    G[, 2*(k), j] <- sqrt(abs(Y[,j]))
  }
}

out <- brandAndBound_anm(Y, G, bs =200, withinAgg = 2, 
                  aggType = 2, alpha = .1,
                  pValueAgg = "tippett",
                  intercept = 1, verbose = T, basis = "poly", M = 4)

out[1,]
dim(out)[1] / factorial(p)


# ancestMat[, c(11,12,13,14,15)]

