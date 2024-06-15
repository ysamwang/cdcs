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


run.onceBnb <- function(p, n, distro, bs = 400, parent_prob = 1/3, verbose = F, cutoff = NULL){
  
  dat <- cdcs::rDAG(p, n, parent_prob = parent_prob, lowScale = .8,
                    highScale = 1, edgeVar = n^(-1/10),
                    dist = distro, uniqueTop = T)
  
    
    Y <- scale(dat$Y)
    
    
    outlingamDirect <- causalXtreme::direct_lingam_search(Y)

    J <- 2
    G1 <- array(0, dim = c(n, J * 2, p) )
    G2 <- array(0, dim = c(n, 3, p) )
    
    for(u in 1:p){
      
      for(j in 1:J){
        
        G1[, 2 * j - 1, u] <- sin(j * Y[,u])
        G1[, 2 * j, u] <- cos(j * Y[,u])
        
      }
      
      G2[, 1, u] <- scale(Y[, u]^2)
      G2[, 2 , u] <- scale(Y[, u]^3)
      G2[, 3 , u] <- scale(sign(Y[,u])*abs(Y[, u])^(2.5))
      
    }
    
    G3 <- abind::abind(G1, G2, along = 2)
    
    rec <- matrix(0, 1, 6)
    colnames(rec) <- c("testFunc", "size", "cover", "ancest", "time","pointEst")
    

    time.rec <- system.time(out <- cdcs::brandAndBound(Y, G3, bs = bs, aggType = 3, alpha = .1,
                                                       pValueAgg = "tippet", verbose = verbose))[3]
    rec[1, ] <- c("comb", sum(out$pValue > .1), all(out[1, -1] == 1:p),
                  mean(cdcs::getAncest(out[which(out$pValue > .1),-1])[lower.tri(dat$B)] == 1 ), time.rec, all(outlingamDirect == 1:p))


    
  return(rec)
}

##################
library(cdcs)

sample.size <- 400
rep.runs <- 10
n.list <- c(500, 1000, 2500, 5000)
d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace", "mixed")
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), d.list)


p <- 10
n <- param.grid[runInd, 1]
distro <- param.grid[runInd, 2]


out <- lapply(1:rep.runs, function(x){run.onceBnb(p, n, distro)})


outTab <- data.frame(p, n, distro, do.call("rbind", out))

write.csv(outTab, paste("../results/bnbMax10/bnbRes_Max_10_",runInd, ".csv", sep = ""))
