runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)









run.onceBnb <- function(p, n, distro, bs = 200, parent_prob = 1/3, verbose = F){
  
  .getAncest <- function(tab){
    
    .estOrder <- function(ind){
      check <- apply(tab, MAR = 1, function(x){which(x == ind[1]) < which(x == ind[2])})
      
      return(mean(check))
    }
    
    p <- dim(tab)[2]
    
    pairs <- combn(1:p, 2)
    checkRes <- apply(pairs, MAR = 2, .estOrder)
    
    A <- matrix(0, p, p)
    A[cbind(pairs[2,], pairs[1,])] <- checkRes
    
    return(A)
  }
  
  
  
  dat <- cdcs::rDAG(p, n, parent_prob = parent_prob, lowScale = .8,
                    highScale = 1, edgeVar = n^(-1/8), dist = distro, uniqueTop = T)
  
  outlingamDirect <- cdcs::directLiNGAM(dat$Y, verbose = verbose, metric = "dhsic")
  

  outTime_1_T_23 <- system.time(out_1_T_23 <- cdcs::brandAndBound(dat$Y, K = c(2, 3), bs = bs, alpha = .1, aggType = 1, pValueAgg = "tippet", verbose = verbose))
  
  rec <- c(sum(out_1_T_23$pValue > .1), all(out_1_T_23[1, -1] == 1:p), outTime_1_T_23[3], mean(.getAncest(out_1_T_23[,-1])[lower.tri(dat$B)] == 1 ) )
  
  names(rec) <- c(paste("oneTippet23", rep(c("size", "truth", "time", "ancest"), times = 8), sep = "_"))
  
  return(rec)
}


##################
library(parallel)
library(cdcs)

sample.size <- 10
rep.runs <- 10

p.list <- c(10, 12, 15, 18)
d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace")
param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list)
### Param grid Size 400



n <- 2500
p <- param.grid[runInd, 1]
distro <- param.grid[runInd, 2]

cl <- makeCluster(3)
clusterExport(cl, ls())
out <- t(parSapply(cl, 1:rep.runs, function(x){run.onceBnb(p, n, distro)}))


outTab <- data.frame(p, n, distro, out)

colnames(outTab) <- c("p","n", "distro",c(paste("oneTippet23", rep(c("size", "truth", "time", "ancest"), times = 8), sep = "_")))


write.csv(outTab, paste("../cl/clRes_",runInd, ".csv", sep = ""))

stopCluster(cl)




