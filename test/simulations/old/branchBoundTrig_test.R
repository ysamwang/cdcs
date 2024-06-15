### Simulations shown in appendix comparing different test functions
# Updated: Feb 27
# Set K = 10 for trig functions
# results with trig, moment, comb written to
# bnbResTrig10_*.csv
#
# Updated Mar 23
# results with trig, moment, zoo, comb written to
# bnbRes_comb_10_*.csv

runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)


testFunc <- function(x, type = 1){
  if(type == 1){
    return(scale(sin(x^2)))
    
  } else if (type == 2){
    return(scale(cos(x^2)))
    
  } else if (type == 3){
    return(scale(sin(x) * x))
    
  } else if (type == 4){
    return(scale(cos(x) * x))
    
  } else if (type == 5){
    return(scale(sin(x^2) * x))
    
    
  } else if (type == 6){
    return(scale(cos(x^2) * x))
    
  } else if (type == 7){
    return(scale(tanh(x)))
  } 
  
}



run.onceBnb <- function(p, n, distro, bs = 400, parent_prob = 1/3, verbose = F, cutoff = NULL){
  
  dat <- cdcs::rDAG(p, n, parent_prob = parent_prob, lowScale = .8,
                    highScale = 1, edgeVar = n^(-1/10),
                    dist = distro, uniqueTop = T)
  
    
    Y <- scale(dat$Y)
    
    
    outlingamDirect <- causalXtreme::direct_lingam_search(Y)

    K <- 25
# 
#     gFunc <- function(x, k){
#       if(k == 1){
#         return(scale(x^2))
#       } else if (k ==2 ){
#         return(scale(x^3))
#       } else if (k ==3 ){
#         return(scale(abs(x)^(2.5) * sign(x)))
#       }
#     }
      
    
    
    G1 <- array(0, dim = c(n, K * 2, p) )
    
    G2 <- array(0, dim = c(n, 3, p) )
    G3 <- array(0, dim = c(n, 3, p) )
    
    for(u in 1:p){
      
      for(k in 1:K){
        
        omega <- rnorm(1)
        G1[, 2 * k - 1, u] <- scale(sin(omega * Y[,u]))
        G1[, 2 * k, u] <- scale(cos(Y[,u] * omega))
        
      }
      
      G2[, 1, u] <- scale(Y[, u]^2)
      G2[, 2 , u] <- scale(Y[, u]^3)
      G2[, 3 , u] <- scale(sign(Y[,u])*abs(Y[, u])^(2.5))
      
      for(z in 1:7){
        G3[, z, u] <- testFunc(Y[, u], type =z)
      }
      
    
      
      
      
    }
    
    G4 <- abind::abind(G1, G2, G3, along = 2)
    
    rec <- matrix(0, 4, 6)
    colnames(rec) <- c("testFunc", "size", "cover", "ancest", "time","pointEst")
    
    
    time.rec <- system.time(out <- cdcs::brandAndBound(Y, G1, bs = bs, withinAgg = 2, aggType = 2, alpha = .1,
                                         pValueAgg = "tippet", verbose = verbose))[3]

    rec[1, ] <- c("trig", sum(out$pValue > .1), all(out[1, -1] == 1:p),
                mean(cdcs::getAncest(out[which(out$pValue > .1),-1])[lower.tri(dat$B)] == 1 ), time.rec, all(outlingamDirect == 1:p))

    time.rec <- system.time(out <- cdcs::brandAndBound(Y, G2, bs = bs, withinAgg = 2, aggType = 2, alpha = .1,
                                                       pValueAgg = "tippet", verbose = verbose))[3]
    rec[2, ] <- c("moment", sum(out$pValue > .1), all(out[1, -1] == 1:p),
                  mean(cdcs::getAncest(out[which(out$pValue > .1),-1])[lower.tri(dat$B)] == 1 ), time.rec, all(outlingamDirect == 1:p))

    time.rec <- system.time(out <- cdcs::brandAndBound(Y, G3, bs = bs, withinAgg = 2, aggType = 2, alpha = .1,
                                                       pValueAgg = "tippet", verbose = verbose))[3]
    rec[3, ] <- c("zoo", sum(out$pValue > .1), all(out[1, -1] == 1:p),
                  mean(cdcs::getAncest(out[which(out$pValue > .1),-1])[lower.tri(dat$B)] == 1 ), time.rec, all(outlingamDirect == 1:p))

    
    time.rec <- system.time(out <- cdcs::brandAndBound(Y, G4, bs = bs, withinAgg = 2, aggType = 2, alpha = .1,
                                         pValueAgg = "tippet", verbose = verbose))[3]

    rec[4, ] <- c("comb", sum(out$pValue > .1), all(out[1, -1] == 1:p),
                mean(cdcs::getAncest(out[which(out$pValue > .1),-1])[lower.tri(dat$B)] == 1 ), time.rec, all(outlingamDirect == 1:p))

    
  return(rec)
}

##################
library(cdcs)

sample.size <- 400
rep.runs <- 5
n.list <- c(500, 1000, 2500, 5000)
d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace", "mixed")
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), d.list)
### Param grid Size 480


p <- 10
n <- param.grid[runInd, 1]
distro <- param.grid[runInd, 2]

out <- lapply(1:rep.runs, function(x){run.onceBnb(p, n, distro)})
outTab <- data.frame(p, n, distro, do.call("rbind", out))

write.csv(outTab, paste("../results/bnbTrig10/bnbRes_comb_10_",runInd, ".csv", sep = ""))


#### Analysis ####
# 

sample.size <- 400
rep.runs <- 5
n.list <- c(500, 1000, 2500, 5000)
d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace", "mixed")
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), d.list)



runInd <- 1
outTabTrig <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/bnbCombNew/bnbRes_comb_10_",runInd, ".csv", sep = ""))[,-1]
missing <- c()
for(runInd in 1:nrow(param.grid)){
  if(file.exists(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/bnbCombNew/bnbRes_comb_10_",runInd, ".csv", sep = ""))){
      temp <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/bnbCombNew/bnbRes_comb_10_",runInd, ".csv", sep = ""))[,-1]
      outTabTrig <- rbind(outTabTrig, temp)
  } else {
    missing <- c(missing, runInd)
  }

}





 
resTabTrig <- aggregate(cbind(size, cover) ~ n + distro + distro + testFunc,
                        FUN = mean, data = outTabTrig)
resTabTrig$size <- resTabTrig$size / factorial(10)

out <- reshape(resTabTrig, idvar = c("distro", "n"),
               direction = "wide",
               timevar = "testFunc")



library(tidyverse)
library(kableExtra)

kbl(out[, c(2,1, )], format = "latex", booktabs = T, align = "c",
    linesep = c(rep("", 3), "\\addlinespace"), digits = 3) %>%
  add_header_above(c(" " = 2, "Prop of unrejected" = 3, "Coverage" = 3, " " = 1)) %>%
  collapse_rows(columns = 1)
# 
