runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)



run.once <- function(p, n, distro, error.sd = 1, corrParam = .8, intercept = T, bs = 200, K = c(2,3)){

  if(distro == "gamma"){
    
    ## Independent Gamma's
    X <- matrix(rgamma(n * p, 1, 1) - 1, nrow = n, ncol = p)
    
    if(intercept){
      X <- cbind(rep(1, n), X)
      beta <- rnorm(p + 1)
    } else {
      beta <- rnorm(p)
    }
    
    Y <- X %*% beta + (rgamma(n, 1, 1) - 1) * error.sd
    
  } else if(distro == "gammaCor"){
    
    ## Correlated Gamma's
    M <- toeplitz(corrParam^(0:(p-1)))
    X <- lcmix::rmvgamma(n, corr = M) - 1
    
    if(intercept){
      X <- cbind(rep(1, n), X)
      beta <- rnorm(p + 1)
    } else {
      beta <- rnorm(p)
    }
    
    Y <- X %*% beta + (rgamma(n, 1, 1) - 1) * error.sd
    
  } else if(distro == "unif"){
    
    X <- matrix(runif(n * p, -sqrt(3), sqrt(3)), nrow = n, ncol = p)
    
    if(intercept){
      X <- cbind(rep(1, n), X)
      beta <- rnorm(p + 1)
    } else {
      beta <- rnorm(p)
    }
    
    
    Y <- X %*% beta + runif(n, -sqrt(3), sqrt(3)) * error.sd
    
  } else if(distro == "unifCor"){
    
    M <- toeplitz(corrParam^(0:(p-1)))
    X <- (MultiRNG::draw.d.variate.uniform(n, p, M) - .5) * 2  * sqrt(3)
    
    if(intercept){
      X <- cbind(rep(1, n), X)
      beta <- rnorm(p + 1)
    } else {
      beta <- rnorm(p)
    }
    
    Y <- X %*% beta + runif(n, -sqrt(3), sqrt(3)) * error.sd
    
  } else if (distro == "laplace"){
    
    X <- matrix(rmutil::rlaplace(n * p, 0, 1/sqrt(2)), nrow = n, ncol = p)
    
    if(intercept){
      X <- cbind(rep(1, n), X)
      beta <- rnorm(p + 1)
    } else {
      beta <- rnorm(p)
    }
    
    Y <- X %*% beta + rmutil::rlaplace(n, 0, 1/sqrt(2)) * error.sd
    
  } else if (distro == "laplaceCor"){
    
    M <- toeplitz(corrParam^(0:(p-1)))
    X <- LaplacesDemon::rmvl(n, mu = rep(0, p), Sigma = M)
    
    if(intercept){
      X <- cbind(rep(1, n), X)
      beta <- rnorm(p + 1)
    } else {
      beta <- rnorm(p)
    }
    
    Y <- X %*% beta + rmutil::rlaplace(n, 0, 1/sqrt(2)) * error.sd
    
  }
  
  
  rec <- rep(0, 6)
  
  if(p < 150){
    
    rec[1:2] <- cdcs::dkw2020(X, Y, bs = bs, statType = "moment", K = c(2,3))
    rec[3] <- cdcs::dkw2020(X, Y, bs = bs, statType = "hsic")
    
    rec[4:5] <- cdcs::senSen2014(X, Y, bs = bs, statType = "moment", K = c(2,3))
    rec[6] <- cdcs::senSen2014(X, Y, bs = bs, statType = "hsic")
    
    rec <- cdcs::convToCont(rec, bs)
    
  } else {
    
    rec[1:2] <- cdcs::dkw2020(X, Y, bs = bs, statType = "moment", K = c(2, 3))
    rec[4:5] <- cdcs::senSen2014(X, Y, bs = bs, statType = "moment", K = c(2, 3))
    rec[c(1,2,4, 5)] <- cdcs::convToCont(rec[c(1,2,4, 5)], bs)
    
  }
  
  
  names(rec) <- c("dkwMoment1", "dkwMoment2", "dkwHSIC", "ssMoment1", "ssMoment2", "ssHSIC")

  return(rec)
}


##################
library(parallel)
library(cdcs)

sample.size <- 1000
rep.runs <- 100
p.list <- c(10, 25, 50, 100, 150, 200)
d.list <- c("gamma", "gammaCor", "laplace", "laplaceCor", "unif", "unifCor")
param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list)
## Param Grid 360

p <- param.grid[runInd, 1]
n <- round(p^(15/8))
distro <- param.grid[runInd, 2]

cl <- makeCluster(3)
clusterExport(cl, ls())

out <- t(parSapply(cl, 1:rep.runs, function(x){run.once(p, n, distro, error.sd = 1, corrParam = .8, intercept = T, bs = 200)}))

outTab <- data.frame(p, n, distro, out)

colnames(outTab) <- c("p", "n", "distro", "dkwMoment", "dkwHSIC", "ssMoment", "ssHSIC")

write.csv(t(out), paste("../lmgdRes/lmgdRes_", runInd, ".csv", sep = ""))


stopCluster(cl)


