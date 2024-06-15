####
# Confidence intervals which account for model selection uncertainty
# results shown in Section 6.2
#
# Update 4/25/24: added max stat
# write files to: ciTest_max_*.csv
#
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

run.onceCI <- function(p, n, distro, bs = 400, parent_prob = 1/3, verbose = F, cutoff = NULL){
  
  
  treatment <- 4
  outcome <- 7
  alpha <- .1
  
  
  dat <- cdcs::rDAG(p, n, parent_prob = 1/3, lowScale = .8,
                    highScale = 1, edgeVar = 1/2,
                    dist = distro, uniqueTop = T)

  
 
  trueParam <- solve(diag(rep(1, p)) - dat$B)[outcome, treatment]
  
  Y <- scale(dat$Y)
  Y.centered <- scale(dat$Y, scale = F)
  
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
  
  out <- cdcs::brandAndBound(Y, G3, bs = bs, aggType = 3, alpha = .1,
                             pValueAgg = "tippet", verbose = verbose)
  
  rec <- rep(0, 4)
  
  if(nrow(out) > 0){
    
    ci_out <- cdcs::ci_modSelect(out, treatment, outcome, effectType = "total", alpha = alpha, Y = Y.centered)
    rec[1:2] <- c(intervals::distance_to_nearest(from = trueParam, to = ci_out$ci) == 0, ci_out$length)
    
  } else{
    
    rec[1:2] <- c(0,0)
  }

  
  
  
  if(which(outlingamDirect == treatment) < which(outlingamDirect == outcome)){
    
    adjPoint <- outlingamDirect[1:which(outlingamDirect == treatment)]
    ii <- which(adjPoint == treatment)
    
    mod <- RcppArmadillo::fastLmPure(Y.centered[, adjPoint, drop = F], Y.centered[, outcome, drop = F])
    mult <- -qt(alpha, df = n - length(adjPoint))
    ci_Point <- c(mod$coeff[ii] - mult * mod$stderr[ii], mod$coeff[ii] + mult * mod$stderr[ii])
    
    rec[3:4] <- c(ci_Point[1] < trueParam & ci_Point[2] > trueParam, ci_Point[2] - ci_Point[1])
    
    trueAdjSet <- c(treatment, which(dat$B[treatment, ] != 0) )
    rec[5] <- all(trueAdjSet %in% adjPoint) & all(adjPoint %in% 1:treatment)
    
  } else {

    
    rec[3:5] <- c(0, 0, 0)  
  }

  


    
  names(rec) <- c("cov_ms", "length_ms", "cov_point", "length_point", "pointEst")
  
  return(rec)
}

##################
library(cdcs)
sample.size <- 400
rep.runs <- 5

n.list <- c(250, 500, 1000, 2000)
d.list <- c("laplace", "gamma")
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), d.list)
### Param grid Size 400


p <- 10
n <- param.grid[runInd, 1]
distro <- param.grid[runInd, 2]

out <- lapply(1:rep.runs, function(x){run.onceCI(p, n, distro)})

outTab <- data.frame(p, n, distro, do.call("rbind", out))

write.csv(outTab, paste("../results/ciTest/ciTest_max_",runInd, ".csv", sep = ""))



