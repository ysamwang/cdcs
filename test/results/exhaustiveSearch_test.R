runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)



run.onceEx <- function(p, n, distro, bs = 200, K = 3, parent_prob = .5){

  dat <- cdcs::rDAG(p, n, parent_prob = parent_prob, lowScale = .6,
              highScale = .95, lowEdge = .3, highEdge = 1.2, dist = distro, uniqueTop = T)
  
  outTime <- system.time(pvals <- cdcs::exhaustiveTest(dat$Y,K = K, bs = bs))
  
  rec <- c(sum(pvals[, 1] > .05), sum(pvals[, 1] > .1), pvals[1, 1],
           sum(pvals[, 2] > .05), sum(pvals[, 2] > .1), pvals[1, 2],
           sum(pvals[, 3] > .05), sum(pvals[, 3] > .1), pvals[1, 3],
           outTime[3])
  names(rec) <- c("totalSetSize05", "totalSetSize10", "totalpValTrue",
                  "sumSetSize05", "sumSetSize10", "sumpValTrue",
                  "minSetSize05", "minSetSize10", "minpValTrue", "time")
  
  return(rec)
}


##################
library(parallel)
library(cdcs)

sample.size <- 1000
rep.runs <- 100
n.list <- c(50, 100, 500, 1000, 2500)
d.list <- c("gauss", "unif", "lognorm", "gamma", "weibull", "laplace")
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), d.list)
## Param grid size 300


p <- 7
n <- param.grid[runInd, 1]
distro <- param.grid[runInd, 2]


out <- t(replicate(rep.runs, run.onceEx(p, n, distro)))

outTab <- data.frame(p, n, distro, out)
names(outTab) <- c("p", "n", "distro", "totalSetSize05", "totalSetSize10", "totalpValTrue",
                   "sumSetSize05", "sumSetSize10", "sumpValTrue",
                   "minSetSize05", "minSetSize10", "minpValTrue", "time")

write.csv(outTab, paste("../exhaustive/exhaustiveRes_", runInd, ".csv", sep = ""))


  

