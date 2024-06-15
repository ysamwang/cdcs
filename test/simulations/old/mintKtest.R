runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)



run.once <- function(p, n, distro, error.sd = 1, corrParam = .6, intercept = T, bs = 200, K = c(2,3), adj = 2, sampleSplit = F){
  
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
  
  regOut <- RcppArmadillo::fastLm(X = X, y = Y)
  
  rec <- rep(0, 8)
  
  k <- max(round(sqrt(n))/2,3)
  keps <- max(round(sqrt(n)) ,5)
  
  rec[1] <- IndepTest::MINTregression(X,Y, k,  keps, w=FALSE, sample(regOut$res, n * bs, replace = T))
  rec[2] <- IndepTest::MINTregression(X,Y, k,  keps, w=FALSE, rnorm(n * bs))
  
  
  k <- max(round(sqrt(n))/4 ,3)
  keps <- max(round(sqrt(n))/2 ,5)
  
  rec[3] <- IndepTest::MINTregression(X,Y, k,  keps, w=FALSE, sample(regOut$res, n * bs, replace = T))
  rec[4] <- IndepTest::MINTregression(X,Y, k,  keps, w=FALSE, rnorm(n * bs))
  

  k <- max(n / 5 ,3)
  keps <- max(n/10 ,5)
  rec[5] <- IndepTest::MINTregression(X,Y, k,  keps, w=FALSE, sample(regOut$res, n * bs, replace = T))
  rec[6] <- IndepTest::MINTregression(X,Y, k,  keps, w=FALSE, rnorm(n * bs))
    
  k <- max(n / 10 ,3)
  keps <- max(n/20 ,5)
  rec[7] <- IndepTest::MINTregression(X,Y, k,  keps, w=FALSE, sample(regOut$res, n * bs, replace = T))
  rec[8] <- IndepTest::MINTregression(X,Y, k,  keps, w=FALSE, rnorm(n * bs))
  
  
  
  names(rec) <- c("sq2res", "sq2param",
                  "sq4res", "sq4param",
                  "d5res", "d5param",
                  "d10res", "d10param")
  
  return(rec)
}


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
                                                        corrParam = .6, intercept = T, bs = 200)}))

outTab <- data.frame(p, n, distro, out)

colnames(outTab) <- c("p", "n", "distro", "sq2res", "sq2param",
                                            "sq4res", "sq4param",
                                            "d5res", "d5param",
                                            "d10res", "d10param")

write.csv(outTab, paste("../results/mint/mint_", runInd, ".csv", sep = ""))

stopCluster(cl)


#### Analysis ####

sample.size <- 1000
rep.runs <- 100
p.list <- c(10, 25, 50, 100)
d.list <- c("gamma", "gammaCor", "laplace", "laplaceCor", "unif", "unifCor", "weibull", "weibullCor")
param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list)

i <- 1
dat <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/mint/mint_", i, ".csv", sep = ""))

for(i in 2:nrow(param.grid)){
  temp <-read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/mint/mint_", i, ".csv", sep = ""))
  dat <- rbind(dat, temp)
}
res <- cbind(aggregate(sq2res < .1 ~ p + n + distro, data = dat, FUN = mean),
             aggregate(sq2param < .1 ~ p + n + distro, data = dat, FUN = mean)[,4],
             aggregate(sq4res < .1 ~ p + n + distro, data = dat, FUN = mean)[,4],
             aggregate(sq4param < .1 ~ p + n + distro, data = dat, FUN = mean)[,4],
             aggregate(d5res < .1 ~ p + n + distro, data = dat, FUN = mean)[,4],
             aggregate(d5param < .1 ~ p + n + distro, data = dat, FUN = mean)[,4],
             aggregate(d10res < .1 ~ p + n + distro, data = dat, FUN = mean)[,4],
             aggregate(d10param < .1 ~ p + n + distro, data = dat, FUN = mean)[,4])


colnames(res) <- c("p", "n", "distro", "sq2res", "sq2param", "sq4res", "sq4param", "d5res", "d5param", "d10res", "d10param")

dis <- "gamma"
p <- 50

hist(dat[which(dat$distro == dis & dat$p == p), 5])
