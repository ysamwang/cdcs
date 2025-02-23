##### Appendix D #####
# 
# Comparison of bootstrap procedure to the limiting Gaussian distribution
#

set.seed(111)
p.list <- c(5, 10, 15, 20, 25)
sim.size <- 1000
rec1 <- rec2 <- rec3 <- rec4 <- matrix(0, sim.size, length(p.list))
rec2 <- matrix(0, sim.size, length(p.list))
rho <- .2
dist <- "ln"
bs <- 500


z <- 1

for(ind in 1:length(p.list)){
  p <- p.list[ind]
  n <- 500
  

  M <- diag(p) + matrix(rho, p, p) - diag(rep(rho, p)) 
  beta <- rep(1, p)
  cat(p)
  cat(": ")

  
  for(i in 1:sim.size){
    
    if(i %% 100 == 0){
      cat(i)
      cat(" ")
    }
    
    X <- scale(exp(mvtnorm::rmvnorm(n, sigma = M)))
    Y <- X %*%  beta + exp(rnorm(n, sd = sqrt(z))) - exp(z/2)
    scale.X2 <- scale(X^2)

    X.int <- cbind(rep(1, n), X)
    H.mat <- t(scale.X2) %*% (diag(n) - X.int %*% solve(t(X.int) %*% X.int, t(X.int)))  
    res <- RcppArmadillo::fastLm(X = X.int, y = Y)$res
    var.est <- sum(res^2/(n-p))
    
    stat <- sum((t(scale.X2) %*% res / sqrt(n))^2)
    # true.sd <- sqrt((exp(1) - 1) * exp(1))
    nullDist1 <- replicate(bs, sum((H.mat %*% rnorm(n, 0, sd = sqrt(var.est)) / sqrt(n))^2))
    nullDist3 <- replicate(500, sum((H.mat %*% (exp(rnorm(n, sd = sqrt(z))) - exp(z/2)) / sqrt(n))^2))
    
  ### Null Hyp
  # ind <- p
  child <- Y
  parents <- X
  
  G <- array(0, dim = c(n , 1, ncol(parents)))
  k <- 1
  for(j in 1:ncol(parents)){
      G[ , k, j] <- scale.X2[, j]  
  }
  
         
    rec1[i, ind] <- (sum(stat < nullDist1) + 1) / (bs + 1) 
    rec2[i, ind] <- (sum(stat < (nullDist1 * n /(n-p))) +1) / (bs + 1)
    rec3[i, ind] <- (sum(stat < nullDist3) +1) / (bs + 1)
    rec4[i, ind] <- cdcs::bnbHelperanm(parents, child, G = G, withinAgg = 3,
                                       aggType = 3, bs = bs, intercept = 1)$pVals[1]
    
    
  }
  cat("\n")
}


getColumnPrint <- function(x, alpha, mult = 1){
  meanx <- colMeans(x <= alpha)
  sdx <- sqrt(meanx * (1 - meanx) / nrow(x))
  
  paste(round(meanx * mult), rep(" (", ncol(x)), 
        round((meanx - 2 * sdx) * mult), rep(", ", ncol(x)),
        round((meanx + 2 * sdx) * mult), rep(")", ncol(x)), sep = "")
}

alpha <- .05
# setEPS()
# postscript("~/Dropbox/Apps/Overleaf/Confidence Sets for Causal Discovery/figures/whyAsympBad_max.eps", width = 6, height = 6)
par(mfrow = c(5,3), mar = c(1, 4, 4, 0))
for(ind in 1:5){
  # hist(rec1[, z], freq = F, main = p, breaks = seq(0, 1, by = .05), )
  # mtext(round(mean(rec1[, z] < .05),3))
  hist(rec2[, ind], freq = F, main = paste("Asymp: p=", p.list[ind], sep = ""), breaks = seq(0, 1, by = .05), xlab = "")
  mtext(round(mean(rec2[, ind] < alpha),3))
  hist(rec3[, ind], freq = F, main = paste("Oracle: p= ", p.list[ind], sep = ""), breaks = seq(0, 1, by = .05), xlab = "")
  mtext(round(mean(rec3[, ind] < alpha),3))
  hist(rec4[, ind], freq = F, main = paste("Proposed: p= ", p.list[ind], sep = ""), breaks = seq(0, 1, by = .05), xlab = "")
  mtext(round(mean(rec4[, ind] < alpha),3))
}
# dev.off()

tab <- rbind(getColumnPrint(rec2, .05, 1000),
      getColumnPrint(rec3, .05, 1000),
      getColumnPrint(rec4, .05, 1000))
rownames(tab) <- c("Asymp", "Oracle", "Proposed")
colnames(tab) <- p.list
xtable::xtable(tab)



##### Appendix C #####
# 
# Naive approach which does not account for estimation of nuissance parameters
#

## setting used in paper ##
# sim.size <- 1000

sim.size <- 50

res100 <-  matrix(0, sim.size, ncol = 10)
res1000 <- matrix(0, sim.size, ncol = 10)

colnames(res100) <- colnames(res1000) <- c("dhsicSize", "tauSize", "dhsicSize_split", "tauSize_split",
                                           "propSize",
                                           "dhsicPower", "tauPower",
                                           "dhsicPower_split", "tauPower_split",
                                           "propPower")

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
  res100[i, 5] <- cdcs::bnbHelperanm(matrix(Y1, ncol = 1), matrix(Y2, ncol = 1)
                                     , G = G[ , , 1, drop = F], withinAgg = 3,
                                     aggType = 3, bs = 500, intercept = 1)$pVals[1]
  
  res100[i, 6] <- dHSIC::dhsic.test(resIncor, Y2)$p
  res100[i, 7] <- TauStar::tauStarTest(resIncor, Y2)$p
  res100[i, 10] <- cdcs::bnbHelperanm(matrix(Y2, ncol = 1), matrix(Y1, ncol = 1)
                                      , G = G[ , , 2, drop = F], withinAgg = 3,
                                      aggType = 3, bs = 500, intercept = 1)$pVals[1]
  
  ## Sample split
  coefCor <- RcppArmadillo::fastLm(X = cbind(rep(1, n/2), Y1[1:(n/2)]), y = Y2[1:(n/2)])$coef
  resCor <-Y2[(n/2 + 1):n] - cbind(rep(1, n/2), Y1[(n/2+1):n]) %*% coefCor
  
  coefIncor <- RcppArmadillo::fastLm(X = cbind(rep(1, n/2), Y2[1:(n/2)]), y = Y1[1:(n/2)])$coef
  resIncor <-Y1[(n/2 + 1):n] - cbind(rep(1, n/2), Y2[(n/2+1):n]) %*% coefIncor
  
  res100[i, 3] <- dHSIC::dhsic.test(resCor, Y1[(n/2+1):n])$p
  res100[i, 4] <- TauStar::tauStarTest(resCor, Y1[(n/2+1):n])$p
  
  res100[i, 8] <- dHSIC::dhsic.test(resIncor, Y2[(n/2+1):n])$p
  res100[i, 9] <- TauStar::tauStarTest(resIncor, Y2[(n/2+1):n])$p

  
  
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
  res1000[i, 5] <- cdcs::bnbHelperanm(matrix(Y1, ncol = 1), matrix(Y2, ncol = 1)
                                     , G = G[ , , 1, drop = F], withinAgg = 3,
                                     aggType = 3, bs = 500, intercept = 1)$pVals[1]
  
  res1000[i, 6] <- dHSIC::dhsic.test(resIncor, Y2)$p
  res1000[i, 7] <- TauStar::tauStarTest(resIncor, Y2)$p
  res1000[i, 10] <- cdcs::bnbHelperanm(matrix(Y2, ncol = 1), matrix(Y1, ncol = 1)
                                      , G = G[ , , 2, drop = F], withinAgg = 3,
                                      aggType = 3, bs = 500, intercept = 1)$pVals[1]
  
  ## Sample split
  coefCor <- RcppArmadillo::fastLm(X = cbind(rep(1, n/2), Y1[1:(n/2)]), y = Y2[1:(n/2)])$coef
  resCor <-Y2[(n/2 + 1):n] - cbind(rep(1, n/2), Y1[(n/2+1):n]) %*% coefCor
  
  coefIncor <- RcppArmadillo::fastLm(X = cbind(rep(1, n/2), Y2[1:(n/2)]), y = Y1[1:(n/2)])$coef
  resIncor <-Y1[(n/2 + 1):n] - cbind(rep(1, n/2), Y2[(n/2+1):n]) %*% coefIncor
  
  res1000[i, 3] <- dHSIC::dhsic.test(resCor, Y1[(n/2+1):n])$p
  res1000[i, 4] <- TauStar::tauStarTest(resCor, Y1[(n/2+1):n])$p
  
  res1000[i, 8] <- dHSIC::dhsic.test(resIncor, Y2[(n/2+1):n])$p
  res1000[i, 9] <- TauStar::tauStarTest(resIncor, Y2[(n/2+1):n])$p
  
  
  if(i %% 10 == 0){
    cat(i)
    cat(" ")
  }
  
}

## Results ##
rbind(colMeans(res100 < .1),
      colMeans(res1000 < .1))