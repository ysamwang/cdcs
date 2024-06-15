runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)






run.once <- function(p, n, distro, error.sd = 1, corrParam = .6,
                     intercept = F, bs = 200, K = c(2,3), sampleSplit = F){

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
    
  } else if (distro == "weibull"){
    a <- 3/4
    
    X <- matrix( (rweibull(n * p, shape = a, 1) - gamma(1 + 1/a)) / sqrt(gamma(1 + 2/a) - gamma(1 + 1/a)^2),
                 nrow = n, ncol = p, byrow = T)
    
    if(intercept){
      X <- cbind(rep(1, n), X)
      beta <- rnorm(p + 1)
    } else {
      beta <- rnorm(p)
    }
    
    Y <- X %*% beta + (rweibull(n, shape = a, 1) - gamma(1 + 1/a)) / sqrt(gamma(1 + 2/a) - gamma(1 + 1/a)^2) * error.sd
    
  } else if (distro == "weibullCor") {
    a <- 3/4
    decay <- 1^(-a)
    
    M <- toeplitz(corrParam^(0:(p-1)))
    X <- (lcmix::rmvweisd(n, shape=a, decay=decay, corr= M) - gamma(1 + 1/a)) / sqrt(gamma(1 + 2/a) - gamma(1 + 1/a)^2)
    
    if(intercept){
      X <- cbind(rep(1, n), X)
      beta <- rnorm(p + 1)
    } else {
      beta <- rnorm(p)
    }
    
    Y <- X %*% beta + (rweibull(n, shape = a, 1) - gamma(1 + 1/a)) / sqrt(gamma(1 + 2/a) - gamma(1 + 1/a)^2) * error.sd
    
  } else if (distro == "bimodalGauss"){
    
    X <- matrix(rnorm(n * p) + sample(c(-2, 2), size = n * p, replace = T), nrow = n, ncol = p)
    
    if(intercept){
      X <- cbind(rep(1, n), X)
      beta <- rnorm(p + 1)
    } else {
      beta <- rnorm(p)
    }
    
    Y <- X %*% beta + (rnorm(n) + sample(c(-2, 2), size = n, replace = T)) * error.sd
    
  }
  
  ### Scale all obs and covariates
  Y <- scale(Y)
  X <- scale(X)
  
  rec <- rep(0, 6)
  

  
  if(p <= 200){
    
    rec[1:3] <- cdcs::dkw2020(X, Y, bs = bs, statType = "moment", K = c(2,3), sampleSplit = sampleSplit, cutoff = log(n) * 3)
    rec[4] <- cdcs::dkw2020(X, Y, bs = bs, statType = "hsic", sampleSplit = sampleSplit, kernel = "gaussian")
    
    rec[5] <- cdcs::senSen2014(X, Y, bs = bs, statType = "moment", K = c(2,3), sampleSplit = sampleSplit)
    rec[6] <- cdcs::senSen2014(X, Y, bs = bs, statType = "hsic", sampleSplit = sampleSplit)
    
    regOut <- RcppArmadillo::fastLm(X = X, y = Y)

    rec[7] <- IndepTest::MINTregression(X,Y, max(round(sqrt(n))/2,3),  max(round(sqrt(n)), 5) ,w=FALSE, sample(regOut$res, n * bs, replace = T))
    rec[8] <- IndepTest::MINTregression(X,Y,  max(round(sqrt(n))/2,3),  max(round(sqrt(n)), 5), w=FALSE, rnorm(n * bs))
    
    # rec <- cdcs::convToCont(rec, bs)
    
  } else {
    
    rec[1:2] <- cdcs::dkw2020(X, Y, bs = bs, statType = "moment", K = c(2, 3), sampleSplit = sampleSplit)
    rec[4:5] <- cdcs::senSen2014(X, Y, bs = bs, statType = "moment", K = c(2, 3), intercept = 0, sampleSplit = sampleSplit)
    # rec[c(1,2,4, 5)] <- cdcs::convToCont(rec[c(1,2,4, 5)], bs)
    
  }
  
  
  names(rec) <- c("dkwMoment1", "dkwMoment2", "dkwMoment3", "dkwHSIC", "ssMoment3", "ssHSIC", "MintRes", "MintParam")

  return(rec)
}

run.once(p = 5, n = 100, "gamma")


##################
library(parallel)
library(cdcs)

sample.size <- 1000
rep.runs <- 100
p.list <- c(10, 25, 50, 100)
d.list <- c("gamma", "gammaCor", "laplace", "laplaceCor", "unif", "unifCor", "weibull", "weibullCor")
param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list)
## Param Grid 480

p <- param.grid[runInd, 1]
n <- round(p^(3/2))
distro <- param.grid[runInd, 2]

cl <- makeCluster(3)
clusterExport(cl, ls())

out <- t(parSapply(cl, 1:rep.runs, function(x){run.once(p, n, distro, error.sd = 1,
                                                        corrParam = .6, intercept = F, bs = 200)}))



outTab <- data.frame(p, n, distro, out)

colnames(outTab) <- c("p", "n", "distro", "dkwMoment1", "dkwMoment2", "dkwHSIC", "ssMoment1", "ssMoment2", "ssHSIC", "MintRes", "MintParam")

write.csv(t(out), paste("../results/lmgd/lmgd_", runInd, ".csv", sep = ""))

stopCluster(cl)


#### Analysis ####
sample.size <- 1000
rep.runs <- 100
p.list <- c(10, 25, 50, 100, 150, 200)
d.list <- c("gamma", "gammaCor", "laplace", "laplaceCor", "unif", "unifCor", "weibull", "weibullCor")
param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list)

i <- 1
dat <- data.frame(p = rep(param.grid[i, 1], rep.runs),
                  n = rep(round(param.grid[i, 1]^(3/2)), rep.runs),
                  dist = rep(param.grid[i, 2], rep.runs),
                  t(read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/lmgdResMint/lmgd_", i, ".csv", sep = ""))[,-1]))

for(i in 2:nrow(param.grid)){
  temp <- data.frame(p = rep(param.grid[i, 1], rep.runs),
                     n = rep(round(param.grid[i, 1]^(3/2)), rep.runs),
                     dist = rep(param.grid[i, 2], rep.runs),
                     t(read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/lmgdResMint/lmgd_", i, ".csv", sep = ""))[,-1]))
  dat <- rbind(dat, temp)
}

colnames(dat) <- c("p", "n", "distro", "dkwMoment1", "dkwMoment2", "dkwHSIC", "ssMoment1", "ssMoment2", "ssHSIC", "MintRes", "MintParam")

dat[is.na(dat)] <- 0

res <- cbind(aggregate(dkwMoment1 < .1 ~ p + n + distro, data = dat, FUN = mean),
             aggregate(dkwMoment2 < .1 ~ p + n + distro, data = dat, FUN = mean)[,4],
              aggregate(dkwHSIC < .1 ~ p + n + distro, data = dat, FUN = mean)[,4],
              aggregate(ssMoment1 < .1 ~ p + n + distro, data = dat, FUN = mean)[,4],
              aggregate(ssMoment2 < .1 ~ p + n + distro, data = dat, FUN = mean)[,4],
              aggregate(ssHSIC < .1 ~ p + n + distro, data = dat, FUN = mean)[,4],
              aggregate(MintRes < .1 ~ p + n + distro, data = dat, FUN = mean)[,4],
              aggregate(MintParam < .1 ~ p + n + distro, data = dat, FUN = mean)[,4])


colnames(res) <- c("p", "n", "distro", "dkwMoment1", "dkwMoment2", "dkwHSIC", "ssMoment1", "ssMoment2", "ssHSIC", "MintRes", "MintParam")
