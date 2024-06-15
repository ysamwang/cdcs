runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)




run.once <- function(p, n, distro, bs = 200, parent_prob = 1/2, verbose = F,
                     cutoff = NULL){
  
  k.list <- c(1:3)
  
  gFunc <- function(x, k){
    if(k == 1){
      return(scale(x^2))
    } else if (k ==2 ){
      return(scale(x^3))
    } else if (k ==3 ){
      return(scale(abs(x)^(2.5) * sign(x)))
    }
    
  }
  
  
  dat <- cdcs::rDAG(p, n, parent_prob = parent_prob, lowScale = .8,
                    highScale = 1, lowEdge = .1, highEdge = 1,
                    dist = distro, uniqueTop = T)
  
  Y <- scale(dat$Y)
  
  ### Null Hyp

  G <- array(0, dim = c(n , length(k.list), p))
  for(j in 1:p){
    for(k in 1:length(k.list)){
      G[ , k, j] <- gFunc(Y[, j], k.list[k])  
    }
  }
  
  
  pvals <- sapply(2:p, FUN = function(z){cdcs::bnbHelperanm(Y[, 1:(z-1),
                                                     drop = F], Y[, z, drop = F], G = G[ , , 1:(z-1), drop = F],
                                                   withinAgg = 2, aggType = 1, bs = bs, intercept = 1)$pVals[1]})
  
  
  
  return(pbeta(min(pvals), 1, p-1))
  
}



##################
library(parallel)
library(cdcs)
sample.size <- 1000
rep.runs <- 200
p.list <- c(10, 15, 20, 30, 45)
d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace", "mixed")
pow.list <- c(2, 5/4)
pp.list <- c(1/2)
param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list, pow.list, pp.list)
nrow(param.grid)
## Param Grid 480


p <- param.grid[runInd, 1]
distro <- param.grid[runInd, 2]
pow.i <- param.grid[runInd, 3]
parent_prob <- param.grid[runInd, 4]
n <- round(p^pow.i)



cl <- makeCluster(7)
# cl <- makeCluster(3)
clusterExport(cl, ls())

out <- t(parSapply(cl, 1:rep.runs, function(x){run.once(p, n, distro, parent_prob = parent_prob, bs = 2000)}))

# out <- t(parSapply(cl, 1:3, function(x){run.once(p, n, distro, parent_prob = parent_prob, bs = 250)}))



outTab <- data.frame(p= rep(p,rep.runs), n, distro, parent_prob, pow.i, t(out))

colnames(outTab) <- c("p", "n", "distro", "parent_prob", "pow","dkwPow1")

write.csv(outTab, paste("../results/singleMod/singleMod_", runInd, ".csv", sep = ""), row.names = F)

stopCluster(cl)

####
library(parallel)
library(cdcs)
sample.size <- 1000
rep.runs <- 200
p.list <- c(10, 15, 20, 30, 45)
d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace", "mixed")
pow.list <- c(2, 5/4)
pp.list <- c(1/2)
param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list, pow.list, pp.list)
nrow(param.grid)


runInd <- 1
outTab <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/singleMod/singleMod_", runInd, ".csv", sep = ""))

colnames(outTab) <- c("p", "n", "distro", "parent_prob", "pow","dkwPow1")

for(runInd in 2:nrow(param.grid)){
  temp <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/singleMod/singleMod_", runInd, ".csv", sep = ""))
  colnames(temp) <- c("p", "n", "distro", "parent_prob", "pow","dkwPow1")

  outTab <- rbind(outTab, temp)
}

out <- aggregate(dkwPow1 <= .1 ~ p + distro + pow, data = outTab, FUN = mean)

