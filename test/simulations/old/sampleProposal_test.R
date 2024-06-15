runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)



run.onceSampleProp <- function(p, n, distro, bs = 200,
                               parent_prob = 1/3, verbose = F,
                               aggType = 3, totalPropCutoff = 2e4, alpha = .1, reCapCut = .1){
  
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
                    highScale = 1, edgeVar = n^(-1/16), dist = distro, uniqueTop = T)
  
  outlingamDirect <- cdcs::directLiNGAM(dat$Y, verbose = verbose, metric = "dhsic")
  
  
  outTime <- system.time(out <- cdcs::bsProp(dat$Y, K = c(2, 3), totalPropCutoff = totalPropCutoff,
                recaptureCutoff = reCapCut, numProposalsEachTime = 100, bs = bs, 
                aggType = aggType, bsMetric = "moment", verbose = verbose))
  
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
           length(trueInProp) == 1, outTime[3], sizes, inc, nrow(out))
  
  names(rec) <- c("pointEst", "trueInProp", "time",
                  "f1Size", "f2Size", "t1Size", "t2Size",
                  "f1Inc", "f2Inc", "t1Inc", "t2Inc", "totalProps")
  
  return(rec)
}


##################
library(cdcs)

sample.size <- 200
rep.runs <- 2

p.list <- c(15)
d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace")
param.grid <- expand.grid(d.list, rep(p.list, each = sample.size / rep.runs))
### Param grid Size 400


distro <- param.grid[runInd, 1]
p <- 12
n <- p^3 * 2



out <- replicate(rep.runs, run.onceSampleProp(p, n, distro, reCapCut = .1))

outTab <- data.frame(p, n, distro, t(out))

colnames(outTab) <- c("p", "n", "distro", "pointEst",
                      "trueInProp", "time",
                      "f1Size", "f2Size", "t1Size", "t2Size",
                      "f1Inc", "f2Inc", "t1Inc", "t2Inc", "totalProps")


write.csv(outTab, paste("../sp/spRes12_",runInd, ".csv", sep = ""))

