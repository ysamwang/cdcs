runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)



run.onceLvl <- function(p, n, distro, bs = 500, K = 3, parent_prob = NULL){
  
  if(is.null(parent_prob)){
    parent_prob <- min(.5, 5/p)
  }
  
  dat <- cdcs::rDAG(p, n, parent_prob = parent_prob, lowScale = .6,
              highScale = .95, lowEdge = 0, highEdge = .9, dist = distro, uniqueTop = T)
  
  outTime <- system.time(pval <- cdcs::fullTest(dat$Y, 1:p, K = K, bs = bs))
  
  ret <- c(pval, outTime[3])
  names(ret) <- c("totalpVal", "sumpVal", "minpVal", "time")
  return(ret)
}


##################
library(parallel)
library(cdcs)


sample.size <- 1000
rep.runs <- 100
p.list <- c(10, 25, 50, 100, 150, 200)
d.list <- c("gauss", "unif", "lognorm", "gamma", "weibull", "laplace")
param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list)
## Param grid size 360

p <- param.grid[runInd, 1]
n <- round(p^(15/8))
distro <- param.grid[runInd, 2]

out <- t(replicate(rep.runs, run.onceLvl(p, n, distro)))

outTab <- data.frame(p, n, distro, out)
names(outTab) <- c("p", "n", "distro", "totalpVal", "sumpVal", "minpVal", "time")


write.csv(outTab, paste("../levelRes/levelRes_", runInd, ".csv", sep = ""))

