runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)

set.seed(100 + runInd)

sim.size <- 10

res100 <-  matrix(0, sim.size, ncol = 6)
res1000 <- matrix(0, sim.size, ncol = 6)

n <- 100
for(i in 1:sim.size){
  Y1 <- rgamma(n, 1, 1) -1
  Y2 <- .5 * Y1 + rgamma(n, 1, 1) - 1  
  
  G <- array(0, dim = c(n, 3, 2))
  G[, , 1] <- cbind(scale(Y1^2), scale(Y1^3), scale(sign(Y1) * abs(Y1)^(2.5))) 
  G[, , 2] <- cbind(scale(Y2^2), scale(Y2^3), scale(sign(Y2) * abs(Y2)^(2.5)))
  
  
  resCor <- RcppArmadillo::fastLm(X = cbind(rep(1, n), Y1), y = Y2)$res
  resIncor <- RcppArmadillo::fastLm(X = cbind(rep(1, n), Y2), y = Y1)$res
  res100[i, 1] <- dHSIC::dhsic.test(resCor, Y1)$p
  res100[i, 2] <- TauStar::tauStarTest(resCor, Y1)$p
  res100[i, 3] <- cdcs::bnbHelperanm(matrix(Y1, ncol = 1), matrix(Y2, ncol = 1)
                                      , G = G[ , , 1, drop = F], withinAgg = 2,
                                      aggType = 2, bs = 500, intercept = 1)$pVals[1]
  
  res100[i, 4] <- dHSIC::dhsic.test(resIncor, Y2)$p
  res100[i, 5] <- TauStar::tauStarTest(resIncor, Y2)$p
  res100[i, 6] <- cdcs::bnbHelperanm(matrix(Y2, ncol = 1), matrix(Y1, ncol = 1)
                                      , G = G[ , , 2, drop = F], withinAgg = 2,
                                      aggType = 2, bs = 500, intercept = 1)$pVals[1]
  
  
  if(i %% 10 == 0){
    cat(i)
    cat(" ")
  }

}

n <- 1000
for(i in 1:sim.size){
  Y1 <- rgamma(n, 1, 1) -1
  Y2 <- .5 * Y1 + rgamma(n, 1, 1) - 1  
  
  G <- array(0, dim = c(n, 3, 2))
  G[, , 1] <- cbind(scale(Y1^2), scale(Y1^3), scale(sign(Y1) * abs(Y1)^(2.5))) 
  G[, , 2] <- cbind(scale(Y2^2), scale(Y2^3), scale(sign(Y2) * abs(Y2)^(2.5)))
  
  
  resCor <- RcppArmadillo::fastLm(X = cbind(rep(1, n), Y1), y = Y2)$res
  resIncor <- RcppArmadillo::fastLm(X = cbind(rep(1, n), Y2), y = Y1)$res
  res1000[i, 1] <- dHSIC::dhsic.test(resCor, Y1)$p
  res1000[i, 2] <- TauStar::tauStarTest(resCor, Y1)$p
  res1000[i, 3] <- cdcs::bnbHelperanm(matrix(Y1, ncol = 1), matrix(Y2, ncol = 1)
                                      , G = G[ , , 1, drop = F], withinAgg = 2,
                                      aggType = 2, bs = 500, intercept = 1)$pVals[1]
  
  res1000[i, 4] <- dHSIC::dhsic.test(resIncor, Y2)$p
  res1000[i, 5] <- TauStar::tauStarTest(resIncor, Y2)$p
  res1000[i, 6] <- cdcs::bnbHelperanm(matrix(Y2, ncol = 1), matrix(Y1, ncol = 1)
                                       , G = G[ , , 2, drop = F], withinAgg = 2,
                                       aggType = 2, bs = 500, intercept = 1)$pVals[1]
  
  
  
  if(i %% 10 == 0){
    cat("i")
    cat(" ")
  }
  
}

colnames(res100) <- colnames(res1000) <- c("dhsicSize", "tauSize", "propSize",
                                           "dhsicPower", "tauPower", "propPower")

write.csv(res100, paste("~/confSetsNew/results/ex1/res100_", runInd, ".csv", sep = ""), row.names = F)
write.csv(res1000, paste("~/confSetsNew/results/ex1/res1000_", runInd, ".csv", sep = ""), row.names = F)
