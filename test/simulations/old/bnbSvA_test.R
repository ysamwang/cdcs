runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)



run.onceSize <- function(p, n, distro, bs = 500, parent_prob = 1/3, verbose = F, uni = T){
  
  .reRunSvA <- function(BInput, scalesInput){
    dat <- cdcs::rDAG(p, n, parent_prob = 1/3, dist = distro, BInput = BInput, scalesInput = scalesInput)
    outlingamDirect <- cdcs::directLiNGAM(dat$Y, verbose = F, metric = "tauStar")
    return(all(outlingamDirect == 1:p))
  }
  
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
  
  
  datOrig <- cdcs::rDAG(p, n, parent_prob = parent_prob, lowScale = .6,
                        highScale = 1, edgeVar = n^(-1/8), dist = distro, uniqueTop = uni)
  
  recoveryProp <- mean(replicate(100, .reRunSvA(datOrig$B, datOrig$scales)))
  
  

  outTime_1_T_23 <- system.time(out_1_T_23 <- cdcs::brandAndBound(datOrig$Y, K = c(2, 3), bs = bs,
                                                                  alpha = .1, aggType = 1, pValueAgg = "tippet", verbose = verbose))
  outTime_2_T_23 <- system.time(out_2_T_23 <- cdcs::brandAndBound(datOrig$Y, K = c(2, 3),
                                                                  bs = bs, alpha = .1, aggType = 2, pValueAgg = "tippet", verbose = verbose))
  
  
  
  rec <- c(sum(out_1_T_23$pValue > .1), all(out_1_T_23[1, -1] == 1:p), outTime_1_T_23[3], mean(.getAncest(out_1_T_23[,-1])[lower.tri(datOrig$B)] == 1 ),
           sum(out_2_T_23$pValue > .1), all(out_2_T_23[1, -1] == 1:p), outTime_2_T_23[3], mean(.getAncest(out_2_T_23[,-1])[lower.tri(datOrig$B)] == 1 ),
           recoveryProp)
  
  names(rec) <- c(paste(rep(c("oneTippet23", "twoTippet23"), each = 4), rep(c("size", "truth", "time", "ancest"), times = 2), sep = "_"),
                  "recoveryProp")
  
  return(rec)
}


##################
library(parallel)
library(cdcs)

sample.size <- 500
rep.runs <- 50

n.list <- c(500, 1000, 1500, 2500)
d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace")
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), d.list)
param.grid <- rbind(param.grid, param.grid)
### Param grid Size 400


p <- 8


n <- param.grid[runInd, 1]
distro <- param.grid[runInd, 2]

cl <- makeCluster(3)
clusterExport(cl, ls())
out <- t(parSapply(cl, 1:rep.runs, function(x){run.onceSize(p, n, distro)}))

outTab <- data.frame(p, n, distro, out)

colnames(outTab) <- c("p","n", "distro", paste(rep(c("oneTippet23", "twoTippet23"), each = 4),
                                               rep(c("size", "truth", "time", "ancest"), times = 2), sep = "_"),
                        "recoveryProp")


write.csv(outTab, paste("../setSize/size_",runInd, ".csv", sep = ""))

stopCluster(cl)



# system.time(out <- run.onceSize(8, n, distro))
# Rchestra::play_song("hello adele")

