runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)



run.onceLvl <- function(p, n, distro, bs = 200, K = c(2,3), parent_prob = 1/3){
  
  dat <- cdcs::rDAG(p, n, parent_prob = 5/p, lowScale = .8,
                    highScale = 1, lowEdge = .2, highEdge = .8,
                    dist = distro, uniqueTop = T)
  
  max(abs(dat$Y))

  
  
  outTime <- system.time(pval <- cdcs::testOrdering(dat$Y, 1:p, K = K, bs = bs, aggType = 3))
  
  ret <- c(fisherInfPval = pval["fisherInfPval"],
           fisherOnePval = pval["fisherOnePval"],
           fisherTwoPval = pval["fisherTwoPval"],
           tipettInfPval = pval["tipettInfPval"],
           tipettOnePval = pval["tipettOnePval"],
           tipettTwoPval = pval["tipettTwoPval"], time = outTime[3])

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

cl <- makeCluster(3)
clusterExport(cl, ls())

out <- t(parSapply(cl, 1:rep.runs, function(x){run.onceLvl(p, n, distro)}))

outTab <- data.frame(p, n, distro, out)



write.csv(outTab, paste("../levelRes/levelRes_", runInd, ".csv", sep = ""))

stopCluster(cl)


run.onceLvl(p, n, distro)
