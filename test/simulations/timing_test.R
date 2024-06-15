### Timing_test.R ###
# Panel 6 in Figure 2 which considers time vs p
#
# Update 4/25/24: modified to use new max stat
# write files: timingRes_max_*.csv
#

runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)





gFunc <- function(x, type){
  
  if(type == 1){
    return(scale(x^2))
  } else if (type ==2 ){
    return(scale(x^3))
  } else if (type ==3 ){
    return(scale(abs(x)^(2.5) * sign(x)))
    
    
  } else if (type == 4){
    return(scale(tanh(x)))
  } else if (type == 5){
    return(scale(tanh(x) * x))
  } else if(type == 6){
    return(scale(tanh(x) * abs(x)^(1.5)))
    
    
    
  } else if (type == 7){
    return(scale(x^2 / (1 + x^2)))
  } else if (type == 8){
    return(scale(x^3/(1 + abs(x^3))))
    
    
  } else if (type == 9){
    return(scale(sin(x * rnorm(1))))
    
  } else if (type == 10){
    return(scale(cos(x * rnorm(1))))
  }
  
}


run.onceTime <- function(p, n, distro, bs = 400, parent_prob = 1/3, verbose = F, cutoff = NULL){
  
  dat <- cdcs::rDAG(p, n, parent_prob = parent_prob, lowScale = .8,
                    highScale = 1, lowEdge = 0, highEdge = 1,
                    dist = distro, uniqueTop = T)
  
  
  Y <- scale(dat$Y)
  
  
  outlingamDirect <- causalXtreme::direct_lingam_search(Y)
  K <- 10
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
    
    G3[, 1, u] <- gFunc(Y[, u], type =4)
    G3[, 2, u] <- gFunc(Y[, u], type =5)
    G3[, 3, u] <- gFunc(Y[, u], type =6)
    
    
    
    
    
    
  }
  
  G4 <- abind::abind(G1, G2, G3, along = 2)
  
  
  time.rec <- system.time(out <- cdcs::brandAndBound(Y, G4, bs = bs, aggType = 3, alpha = .1,
                                                     pValueAgg = "tippet", verbose = verbose))[3]
  
  rec <- c(sum(out$pValue > .1), all(out[1, -1] == 1:p),
           mean(cdcs::getAncest(out[which(out$pValue > .1),-1])[lower.tri(dat$B)] == 1 ), time.rec, all(outlingamDirect == 1:p))
  
  names(rec) <- c("size", "cover", "ancest", "time", "pointEst")
  
  return(rec)
}


##################
library(cdcs)

sample.size <- 10
rep.runs <- 1
# For 10, 12, 14, 16
# p.list <- seq(10, 20, by = 2)

# For 18, 20 need larger memory
p.list <- c(10, 12, 14, 16, 18, 20)
d.list <- c("gamma")
param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list)
### Param grid Size 60



n <- 10000
p <- param.grid[runInd, 1]
distro <- param.grid[runInd, 2]


out <- run.onceTime(p, n, distro, verbose = T)


outTab <- data.frame(p, n, distro, t(out))

colnames(outTab) <- c("p","n", "distro", c("size", "cover", "ancest", "time", "pointEst"))


write.csv(outTab, paste("../results/timing/timingRes_max_",runInd, ".csv", sep = ""))


