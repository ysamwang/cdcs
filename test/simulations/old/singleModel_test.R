runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)


run.once <- function(p, n, distro, bs = 500){
  dat <- cdcs::rDAG(p, n, parent_prob = parent_prob, lowScale = .8,
                    highScale = 1, lowEdge = .1, highEdge = .95,
                    dist = distro, uniqueTop = T)
  
  Y <- scale(dat$Y)
  
}

run.once(p = 5, n = 100, "gamma")


##################
library(parallel)
library(cdcs)

sample.size <- 500
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

colnames(outTab) <- c("p", "n", "distro", "true_dkwMoment3", "true_senSen", "true_rpTest", "false_dkwMoment3", "false_senSen", "false_rpTest")

write.csv(t(out), paste("../results/singleMod/singleMod_", runInd, ".csv", sep = ""))

stopCluster(cl)


# #### Analysis ####
# sample.size <- 500
# rep.runs <- 100
# p.list <- c(10, 25, 50, 100, 150, 200)
# d.list <- c("gamma", "gammaCor", "laplace", "laplaceCor", "unif", "unifCor", "weibull", "weibullCor")
# param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list)
# 
# i <- 1
# dat <- data.frame(p = rep(param.grid[i, 1], rep.runs),
#                   n = rep(round(param.grid[i, 1]^(3/2)), rep.runs),
#                   dist = rep(param.grid[i, 2], rep.runs),
#                   t(read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/lmgdResMint/lmgd_", i, ".csv", sep = ""))[,-1]))
# 
# for(i in 2:nrow(param.grid)){
#   temp <- data.frame(p = rep(param.grid[i, 1], rep.runs),
#                      n = rep(round(param.grid[i, 1]^(3/2)), rep.runs),
#                      dist = rep(param.grid[i, 2], rep.runs),
#                      t(read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/lmgdResMint/lmgd_", i, ".csv", sep = ""))[,-1]))
#   dat <- rbind(dat, temp)
# }
# 
# colnames(dat) <- c("p", "n", "distro", "dkwMoment1", "dkwMoment2", "dkwHSIC", "ssMoment1", "ssMoment2", "ssHSIC", "MintRes", "MintParam")
# 
# dat[is.na(dat)] <- 0
# 
# res <- cbind(aggregate(dkwMoment1 < .1 ~ p + n + distro, data = dat, FUN = mean),
#              aggregate(dkwMoment2 < .1 ~ p + n + distro, data = dat, FUN = mean)[,4],
#              aggregate(dkwHSIC < .1 ~ p + n + distro, data = dat, FUN = mean)[,4],
#              aggregate(ssMoment1 < .1 ~ p + n + distro, data = dat, FUN = mean)[,4],
#              aggregate(ssMoment2 < .1 ~ p + n + distro, data = dat, FUN = mean)[,4],
#              aggregate(ssHSIC < .1 ~ p + n + distro, data = dat, FUN = mean)[,4],
#              aggregate(MintRes < .1 ~ p + n + distro, data = dat, FUN = mean)[,4],
#              aggregate(MintParam < .1 ~ p + n + distro, data = dat, FUN = mean)[,4])
# 
# 
# colnames(res) <- c("p", "n", "distro", "dkwMoment1", "dkwMoment2", "dkwHSIC", "ssMoment1", "ssMoment2", "ssHSIC", "MintRes", "MintParam")
