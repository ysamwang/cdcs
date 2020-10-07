runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)



run.onceBnb <- function(p, n, distro, bs = 200, K = 3, aggType = "sum", parent_prob = .5, verbose = F, uni = T){

  dat <- cdcs::rDAG(p, n, parent_prob = parent_prob, lowScale = .6,
                highScale = .95, lowEdge = .3, highEdge = 1.2, dist = distro, uniqueTop = uni)

  
  outTime <- system.time(out <- cdcs::brandAndBound(dat$Y, K = 3, bs = bs, alpha = .8, aggType = aggType, verbose = verbose))
  
  
  rec <- c(sum(out$pValue > .05), sum(out$pValue > .1),
           all(out[1, -1] == 1:p),
           all(out[1, -1] == 1:p) * out[1, 1],
           outTime[3]) 
  
  return(rec)
}


##################
library(parallel)
library(cdcs)

sample.size <- 1000
n.list <- c(100, 500, 1000, 2500)
d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace")

agg.list <- c("sum", "min")
param.grid <- expand.grid(n.list, d.list, agg.list)
### Param grid Size 48

cl <- makeCluster(31)


p <- 10
n <- param.grid[runInd, 1]
distro <- param.grid[runInd, 2]
aggType <- param.grid[runInd, 3]

clusterExport(cl, ls())
clusterEvalQ(cl = cl, library(cdcs))


out <- parSapply(cl, 1:sample.size, function(X){run.onceBnb(p, n, distro, aggType = aggType)})

outTab <- t(out)
colnames(outTab) <- c("setSize05", "setSize10", "trueFirst", "pValTrue", "time")


write.csv(outTab, paste("../bnb/bnbRes_",p, "p_", n,"n_",aggType, "_", distro, ".csv", sep = ""))
#
stopCluster(cl)


##################
# runInd <- 5
# library(parallel)
# library(cdcs)
# 
# sample.size <- 20
# n.list <- c(100, 500, 1000, 2500)
# d.list <- c("gauss", "unif", "lognorm", "gamma", "weibull", "laplace")
# 
# agg.list <- c("sum", "min")
# param.grid <- expand.grid(n.list, d.list, agg.list)
# ### Param grid Size 48
# 
# cl <- makeCluster(2)
# 
# 
# p <- 4
# n <- param.grid[runInd, 1]
# distro <- param.grid[runInd, 2]
# aggType <- param.grid[runInd, 3]
# 
# clusterExport(cl, ls())
# clusterEvalQ(cl = cl, library(cdcs))
# 
# 
# out <- parSapply(cl, 1:sample.size, function(X){run.onceBnb(p, n, distro, aggType = aggType)})
# 
# outTab <- t(out)
# colnames(outTab) <- c("setSize05", "setSize10", "trueFirst", "pValTrue", "time")
# 
# 
# # write.csv(outTab, paste("../bnb/bnbRes_",p, "p_", n,"n_",aggType, "_", distro, ".csv", sep = ""))
# 
# stopCluster(cl)

