runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)



run.onceSampleProp <- function(p, n, distro, bs = 200,
                               parent_prob = 1/3, verbose = F,
                               aggType = 3, propBsCutoff = 500, alpha = .1){
  
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
  
  
  out <- cdcs::bsProp(dat$Y, K = c(2, 3), propBsCutoff = propBsCutoff,
                recaptureCutoff = .1, proposalBs = 250, bs = 200, aggType = aggType)
  
  sizes <- c(sum(out[, "fisherOnePval"] > alpha),
            sum(out[, "fisherTwoPval"] > alpha),
            sum(out[, "tipettOnePval"] > alpha),
            sum(out[, "tipettTwoPval"] > alpha))
  
  trueInProp <- which(apply(out[, -c(1:6)], MAR = 1, function(x){all(x == 1:p)}))
  if(length(trueInProp) == 1){
    inc <- out[trueInProp, c("fisherOnePval", "fisherTwoPval", "tipettOnePval", "tipettTwoPval") ] > alpha
  } else {
    inc <- rep(0, 4)
  } 
               
  rec <- c(all(outlingamDirect == 1:p), 
           length(trueInProp) == 1, outTime[3], sizes, inc)
  
  names(rec) <- c("pointEst", "trueInProp", "time",
                  "f1Size", "f2Size", "t1Size", "t2Size",
                  "f1Inc", "f2Inc", "t1Inc", "t2Inc")
  
  return(rec)
}


##################
library(parallel)
library(cdcs)

sample.size <- 1000
rep.runs <- 100

p.list <- c(15, 20, 25)
n <- 1000
d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace")
param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list)
### Param grid Size 400


p <- param.grid[runInd, 1]
distro <- param.grid[runInd, 2]

cl <- makeCluster(3)
clusterExport(cl, ls())
out <- t(parSapply(cl, 1:rep.runs, function(x){run.onceSampleProp(p, n, distro, propBsCutoff = 1000)}))

run.onceSampleProp(p = 12, n = 1000, distro = "gamma", propBsCutoff = 200)

outTab <- data.frame(p, n, distro, out)

colnames(outTab) <- c("p", "n", "distro",
                      "trueInProp", "time",
                      "f1Size", "f2Size", "t1Size", "t2Size",
                      "f1Inc", "f2Inc", "t1Inc", "t2Inc")


write.csv(outTab, paste("../sp/spRes_",runInd, ".csv", sep = ""))

stopCluster(cl)




