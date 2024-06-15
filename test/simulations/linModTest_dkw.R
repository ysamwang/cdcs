####
# test for single regression 
#
# Update 4/25/24: added max stat
# write files to: powTest_max_*.csv
#
#

runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)


run.once <- function(p, n, distro, bs = 500, parent_prob = 1/2, verbose = F,
                         cutoff = NULL){
  
  
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
      return(scale(tanh(x) * x^2))
      
      
      
    } else if (type == 7){
      return(scale(tanh(x) * sin(x)))
    } else if (type == 8){
      return(scale(tanh(x) * cos(x)))
      
      
    } else if (type == 9){
      return(sin(x * rnorm(1)))
      
    } else if (type == 10){
      return(cos(x * rnorm(1)))
    }
    
  }
  
  
  
  rec <- matrix(0, nrow = 8, ncol = 3)
  colnames(rec) <- paste(c("null_", "alt_", "timeA_"), "M", sep = "")
  
  dat <- cdcs::rDAG(p, n, parent_prob = parent_prob, lowScale = .8,
                    highScale = 1, lowEdge = .1, highEdge = .95,
                    dist = distro, uniqueTop = T)
  
  Y <- scale(dat$Y)
  
  ### Null Hyp
  ind <- p
  child <- Y[, ind, drop = F]
  parents <- Y[, -ind, drop = F]
  
  
  J1 <- 2
  J2 <- 10
  G1 <- array(0, dim = c(n, J1 * 2, ncol(parents)) )
  G2 <- array(0, dim = c(n, J2 * 2, ncol(parents)) )
  G3 <- array(0, dim = c(n, 3, ncol(parents)) )
  G4 <- array(0, dim = c(n, 3, ncol(parents)) )
  G5 <- array(0, dim = c(n, 5, ncol(parents)) )
  

  
  for(u in 1:ncol(parents)){
    
    for(j in 1:J2){
      
      if(j <= J1){
        G1[, 2 * j - 1, u] <- sin(parents[, u]*j)
        G1[, 2 * j, u] <- cos(parents[, u]*j)
      }

      omega <- rnorm(1)
      G2[, 2 * j - 1, u] <- sin(omega * parents[, u])
      G2[, 2 * j, u] <- cos(omega * parents[, u])
      
    }
    
    G3[, 1, u] <- scale(parents[, u]^2, scale = F)
    G3[, 2 , u] <- scale(parents[, u]^3, scale = F)
    G3[, 3 , u] <- scale(sign(parents[,u])*abs(parents[, u])^(2.5), scale = F)
    
    G4[, 1, u] <- scale(parents[, u]^2)
    G4[, 2 , u] <- scale(parents[, u]^3)
    G4[, 3 , u] <- scale(sign(parents[,u])*abs(parents[, u])^(2.5))
    
    
    G5[, 1, u] <- gFunc(parents[, u], 4)
    G5[, 2 , u] <- gFunc(parents[, u], 5)
    G5[, 3 , u] <- gFunc(parents[, u], 6)
    G5[, 4 , u] <- gFunc(parents[, u], 7)
    G5[, 5 , u] <- gFunc(parents[, u], 8)
    
  
  }
  
  G6 <- abind::abind(G1, G4, along = 2)
  G7 <- abind::abind(G2, G4, along = 2)
  G8 <- abind::abind(G5, G4, along = 2)
  
  
  rec[1, 1] <- cdcs::bnbHelperanm(parents, child, G = G1, withinAgg = 3, aggType = 3, bs = bs, intercept = 1)$pVals[1]
  rec[2, 1] <- cdcs::bnbHelperanm(parents, child, G = G2, withinAgg = 3, aggType = 3, bs = bs, intercept = 1)$pVals[1]
  rec[3, 1] <- cdcs::bnbHelperanm(parents, child, G = G3, withinAgg = 3, aggType = 3, bs = bs, intercept = 1)$pVals[1]
  rec[4, 1] <- cdcs::bnbHelperanm(parents, child, G = G4, withinAgg = 3, aggType = 3, bs = bs, intercept = 1)$pVals[1]
  rec[5, 1] <- cdcs::bnbHelperanm(parents, child, G = G5, withinAgg = 3, aggType = 3, bs = bs, intercept = 1)$pVals[1]
  rec[6, 1] <- cdcs::bnbHelperanm(parents, child, G = G6, withinAgg = 3, aggType = 3, bs = bs, intercept = 1)$pVals[1]
  rec[7, 1] <- cdcs::bnbHelperanm(parents, child, G = G7, withinAgg = 3, aggType = 3, bs = bs, intercept = 1)$pVals[1]
  rec[8, 1] <- cdcs::bnbHelperanm(parents, child, G = G8, withinAgg = 3, aggType = 3, bs = bs, intercept = 1)$pVals[1]
  
  
    
    
  
  
  
  ### Alt Hyp
  ind <- 1
  child <- Y[, ind, drop = F]
  parents <- Y[, -ind, drop = F]
  
  J1 <- 2
  J2 <- 10
  G1 <- array(0, dim = c(n, J1 * 2, ncol(parents)) )
  G2 <- array(0, dim = c(n, J2 * 2, ncol(parents)) )
  G3 <- array(0, dim = c(n, 3, ncol(parents)) )
  G4 <- array(0, dim = c(n, 3, ncol(parents)) )
  G5 <- array(0, dim = c(n, 5, ncol(parents)) )
  
  
  
  for(u in 1:ncol(parents)){
    
    for(j in 1:J2){
      
      if(j <= J1){
        G1[, 2 * j - 1, u] <- sin(parents[, u]*j)
        G1[, 2 * j, u] <- cos(parents[, u]*j)
      }
      
      omega <- rnorm(1)
      G2[, 2 * j - 1, u] <- sin(omega * parents[, u])
      G2[, 2 * j, u] <- cos(omega * parents[, u])
      
    }
    
    G3[, 1, u] <- scale(parents[, u]^2, scale = F)
    G3[, 2 , u] <- scale(parents[, u]^3, scale = F)
    G3[, 3 , u] <- scale(sign(parents[,u])*abs(parents[, u])^(2.5), scale = F)
    
    G4[, 1, u] <- scale(parents[, u]^2)
    G4[, 2 , u] <- scale(parents[, u]^3)
    G4[, 3 , u] <- scale(sign(parents[,u])*abs(parents[, u])^(2.5))
    
    
    G5[, 1, u] <- gFunc(parents[, u], 4)
    G5[, 2 , u] <- gFunc(parents[, u], 5)
    G5[, 3 , u] <- gFunc(parents[, u], 6)
    G5[, 4 , u] <- gFunc(parents[, u], 7)
    G5[, 5 , u] <- gFunc(parents[, u], 8)
    
    
  }
  
  G6 <- abind::abind(G1, G4, along = 2)
  G7 <- abind::abind(G2, G4, along = 2)
  G8 <- abind::abind(G5, G4, along = 2)


  rec[1, 3] <- microbenchmark::microbenchmark(rec[1, 2] <- cdcs::bnbHelperanm(parents, child, G = G1, withinAgg = 3, aggType = 3, bs = bs, intercept = 1)$pVals[1], times = 1)$time
  rec[2, 3] <- microbenchmark::microbenchmark(rec[2, 2] <- cdcs::bnbHelperanm(parents, child, G = G2, withinAgg = 3, aggType = 3, bs = bs, intercept = 1)$pVals[1], times = 1)$time
  rec[3, 3] <- microbenchmark::microbenchmark(rec[3, 2] <- cdcs::bnbHelperanm(parents, child, G = G3, withinAgg = 3, aggType = 3, bs = bs, intercept = 1)$pVals[1], times = 1)$time
  rec[4, 3] <- microbenchmark::microbenchmark(rec[4, 2] <- cdcs::bnbHelperanm(parents, child, G = G4, withinAgg = 3, aggType = 3, bs = bs, intercept = 1)$pVals[1], times = 1)$time
  rec[5, 3] <- microbenchmark::microbenchmark(rec[5, 2] <- cdcs::bnbHelperanm(parents, child, G = G5, withinAgg = 3, aggType = 3, bs = bs, intercept = 1)$pVals[1], times = 1)$time
  rec[6, 3] <- microbenchmark::microbenchmark(rec[6, 2] <- cdcs::bnbHelperanm(parents, child, G = G6, withinAgg = 3, aggType = 3, bs = bs, intercept = 1)$pVals[1], times = 1)$time
  rec[7, 3] <- microbenchmark::microbenchmark(rec[7, 2] <- cdcs::bnbHelperanm(parents, child, G = G7, withinAgg = 3, aggType = 3, bs = bs, intercept = 1)$pVals[1], times = 1)$time
  rec[8, 3] <- microbenchmark::microbenchmark(rec[8, 2] <- cdcs::bnbHelperanm(parents, child, G = G8, withinAgg = 3, aggType = 3, bs = bs, intercept = 1)$pVals[1], times = 1)$time
 
  
  
  rec <- cbind(TestFunc = c("trig", "trigRand", "momentCenter", "momentScale", "zoo", "comb_TM", "comb_TRM", "comb_ZM"), rec)
  
  return(rec)
  
}

  

##################
library(cdcs)
sample.size <- 500
rep.runs <- 50
p.list <- c(10, 15, 20, 30, 45)
d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace", "mixed")
pow.list <- c(2, 5/4)
pp.list <- c(1/2)
param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list, pow.list, pp.list)
nrow(param.grid)
## Param Grid 480


p <- param.grid[runInd, 1]
distro <- param.grid[runInd, 2]
pow.i <- param.grid[runInd, 3]
parent_prob <- param.grid[runInd, 4]
n <- round(p^pow.i)


out <- lapply(1:rep.runs, function(x){run.once(p, n, distro, parent_prob = parent_prob, bs = 500)})
outTab <- data.frame(p, n, distro, parent_prob, pow.i, do.call("rbind", out))

write.csv(outTab, paste("../results/lmgdPow/lmgdDAG_dkw_max_", runInd, ".csv", sep = ""), row.names = F)


#####################
# sample.size <- 500
# rep.runs <- 50
# p.list <- c(10, 20, 45)
# d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace", "mixed")
# pow.list <- c(2, 5/4)
# pp.list <- c(1/2)
# param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list, pow.list, pp.list)
# 
# 
# runInd <- 1
# outTab <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/lmgdPow/powTest_max_", runInd, ".csv", sep = ""))
# 
# for(runInd in 2:nrow(param.grid)){
#   temp <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/lmgdPow/powTest_max_", runInd, ".csv", sep = ""))
#   outTab <- rbind(outTab, temp)
# }
# 
# 
# 
# resNull <- aggregate(cbind(null_O1 < .1,
#                            null_O2 < .1,
#                            null_M < .1)~ p + distro + pow.i + TestFunc, data = outTab, FUN = mean)
# colnames(resNull)[5:7] <- c("nullO1", "nullO2", "nullM")
# 
# res <- aggregate(cbind(alt_O1 < .1,
#                        alt_O2 < .1,
#                        alt_M < .1)~ p + distro + pow.i + TestFunc, data = outTab, FUN = mean)
# 
# colnames(res)[5:7] <- c("altO1", "altO2", "altM")
# 
# out <- reshape(res, idvar = c("distro", "p", "pow.i"),
#                direction = "wide",
#                timevar = "TestFunc")
