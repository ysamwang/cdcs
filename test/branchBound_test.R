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
    
    
    outTime_1_F_3 <- system.time(out_1_F_3 <- cdcs::brandAndBound(dat$Y, K = c(3), bs = bs, alpha = .1, aggType = 1, pValueAgg = "fisher", verbose = verbose))
    outTime_1_T_3 <- system.time(out_1_T_3 <- cdcs::brandAndBound(dat$Y, K = c(3), bs = bs, alpha = .1, aggType = 1, pValueAgg = "tippet", verbose = verbose))
    
    outTime_1_F_23 <- system.time(out_1_F_23 <- cdcs::brandAndBound(dat$Y, K = c(2, 3), bs = bs, alpha = .1, aggType = 1, pValueAgg = "fisher", verbose = verbose))
    outTime_1_T_23 <- system.time(out_1_T_23 <- cdcs::brandAndBound(dat$Y, K = c(2, 3), bs = bs, alpha = .1, aggType = 1, pValueAgg = "tippet", verbose = verbose))
    
    outTime_2_F_3 <- system.time(out_2_F_3 <- cdcs::brandAndBound(dat$Y, K = c(3), bs = bs, alpha = .1, aggType = 2, pValueAgg = "fisher", verbose = verbose))
    outTime_2_T_3 <- system.time(out_2_T_3 <- cdcs::brandAndBound(dat$Y, K = c(3), bs = bs, alpha = .1, aggType = 2, pValueAgg = "tippet", verbose = verbose))
    
    outTime_2_F_23 <- system.time(out_2_F_23 <- cdcs::brandAndBound(dat$Y, K = c(2, 3), bs = bs, alpha = .1, aggType = 2, pValueAgg = "fisher", verbose = verbose))
    outTime_2_T_23 <- system.time(out_2_T_23 <- cdcs::brandAndBound(dat$Y, K = c(2, 3), bs = bs, alpha = .1, aggType = 2, pValueAgg = "tippet", verbose = verbose))
    
  
  
  rec <- c(sum(out_1_F_3$pValue > .1), all(out_1_F_3[1, -1] == 1:p), outTime_1_F_3[3], mean(.getAncest(out_1_F_3[,-1])[lower.tri(dat$B)] == 1 ),
           sum(out_1_T_3$pValue > .1), all(out_1_T_3[1, -1] == 1:p), outTime_1_T_3[3], mean(.getAncest(out_1_T_3[,-1])[lower.tri(dat$B)] == 1 ),
           sum(out_1_F_23$pValue > .1), all(out_1_F_23[1, -1] == 1:p), outTime_1_F_23[3], mean(.getAncest(out_1_F_23[,-1])[lower.tri(dat$B)] == 1 ),
           sum(out_1_T_23$pValue > .1), all(out_1_T_23[1, -1] == 1:p), outTime_1_T_23[3], mean(.getAncest(out_1_T_23[,-1])[lower.tri(dat$B)] == 1 ),
           
           sum(out_2_F_3$pValue > .1), all(out_2_F_3[1, -1] == 1:p), outTime_2_F_3[3], mean(.getAncest(out_2_F_3[,-1])[lower.tri(dat$B)] == 1 ),
           sum(out_2_T_3$pValue > .1), all(out_2_T_3[1, -1] == 1:p), outTime_2_T_3[3], mean(.getAncest(out_2_T_3[,-1])[lower.tri(dat$B)] == 1 ),
           sum(out_2_F_23$pValue > .1), all(out_2_F_23[1, -1] == 1:p), outTime_2_F_23[3], mean(.getAncest(out_2_F_23[,-1])[lower.tri(dat$B)] == 1 ),
           sum(out_2_T_23$pValue > .1), all(out_2_T_23[1, -1] == 1:p), outTime_2_T_23[3], mean(.getAncest(out_2_T_23[,-1])[lower.tri(dat$B)] == 1 ),
           
           all(outlingamDirect == 1:p))
  
  names(rec) <- c(paste(rep(c("oneFisher3", "oneTippet3", "oneFisher23", "oneTippet23",
                  "twoFisher3", "twoTippet3", "twoFisher23", "twoTippet23"), each = 3), rep(c("size", "truth", "time", "ancest"), times = 8), sep = "_"),
                  "pointEst")
  
  return(rec)
}


##################
library(parallel)
library(cdcs)

sample.size <- 1000
rep.runs <- 50

n.list <- c(500, 1000, 1500, 2500)
d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace")
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), d.list)
### Param grid Size 400


p <- 8

n <- param.grid[runInd, 1]
distro <- param.grid[runInd, 2]

cl <- makeCluster(3)
clusterExport(cl, ls())
out <- t(parSapply(cl, 1:rep.runs, function(x){run.onceBnb(p, n, distro)}))


outTab <- data.frame(p, n, distro, out)

colnames(outTab) <- c("p","n", "distro",
                      c(paste(rep(c("oneFisher3", "oneTippet3", "oneFisher23", "oneTippet23",
                                    "twoFisher3", "twoTippet3", "twoFisher23", "twoTippet23"), each = 3), rep(c("size", "truth", "time", "ancest"), times = 8), sep = "_"),
                        "pointEst"))


write.csv(outTab, paste("../bnb/bnbRes_",runInd, ".csv", sep = ""))

stopCluster(cl)




