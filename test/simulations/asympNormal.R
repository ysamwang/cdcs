##################
# Supplement D
#
# Comparison of bootstrap procedure to the limiting Gaussian distribution
#

set.seed(1111)
p.list <- c(5, 10, 15, 20, 25)
sim.size <- 2000
rec1 <- rec2 <- rec3 <- rec4 <- matrix(0, sim.size, length(p.list))
rec2 <- matrix(0, sim.size, length(p.list))
rho <- .2
dist <- "ln"
bs <- 500


z <- log(1 + sqrt(5)) - log(2)

# par(mfrow = c(5,4))
for(z in 1:length(p.list)){
  p <- p.list[z]
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

    # X <- lcmix::rmvgamma(n, corr = M)
    # Y <- X %*%  beta + rgamma(n, 2, 1) - 2
    # scale.X2 <- scale(X^2)
    
    # X <- scale(LaplacesDemon::rmvl(n, mu = rep(0, p), Sigma = M))
    # Y <- X %*%  beta + rlaplace(n, 0, 1)    
    # scale.X2 <- scale(X^3)
    
    X.int <- cbind(rep(1, n), X)
    H.mat <- t(scale.X2) %*% (diag(n) - X.int %*% solve(t(X.int) %*% X.int, t(X.int)))  
    res <- RcppArmadillo::fastLm(X = X.int, y = Y)$res
    sd.est <- sum(res^2/(n-p))
    
    stat <- sum((t(scale.X2) %*% res / sqrt(n))^2)
    # true.sd <- sqrt((exp(1) - 1) * exp(1))
    nullDist1 <- replicate(bs, sum((H.mat %*% rnorm(n, 0, sd = sqrt(sd.est)) / sqrt(n))^2))
    # nullDist2 <- replicate(500, sum((H.mat %*% rnorm(n, 0, sd = sqrt(sd.est)) / sqrt(n-p))^2))
    nullDist3 <- replicate(500, sum((H.mat %*% (exp(rnorm(n, sd = sqrt(z))) - exp(z/2)) / sqrt(n))^2))
    # nullDist3 <- replicate(500, sum((H.mat %*% rlaplace(n, 0, 1) / sqrt(n))^2))
    
      ### Null Hyp
  # ind <- p
  child <- Y
  parents <- X
  
  G <- array(0, dim = c(n , 1, ncol(parents)))
  k <- 1
  for(j in 1:ncol(parents)){
      G[ , k, j] <- scale.X2[, j]  
  }
  
         
    rec1[i, z] <- (sum(stat < nullDist1) + 1) / (bs + 1) 
    rec2[i, z] <- (sum(stat < (nullDist1 * n /(n-p))) +1) / (bs + 1)
    rec3[i, z] <- (sum(stat < nullDist3) +1) / (bs + 1)
    rec4[i, z] <- cdcs::bnbHelperanm(parents, child, G = G, withinAgg = 3, aggType = 3, bs = bs, intercept = 1)$pVals[1]
    
    
  }
  # cat("\n")
  # hist(rec1[, z], freq = F, main = p, breaks = seq(0, 1, by = .05))
  # mtext(round(mean(rec1[, z] < .05),3))
  # hist(rec2[, z], freq = F, main = p, breaks = seq(0, 1, by = .05))
  # mtext(round(mean(rec2[, z] < .05),3))
  # hist(rec3[, z], freq = F, main = p, breaks = seq(0, 1, by = .05))
  # mtext(round(mean(rec3[, z] < .05),3))
  # hist(rec4[, z], freq = F, main = p, breaks = seq(0, 1, by = .05))
  # mtext(round(mean(rec4[, z] < .05),3))
  cat("\n")
}


# colMeans(rec1 < .1)
# colMeans(rec2 < .1)
# colMeans(rec3 < .1)
# colMeans(rec4 < .1)

getColumnPrint <- function(x, alpha, mult = 1){
  meanx <- colMeans(x <= alpha)
  sdx <- sqrt(meanx * (1 - meanx) / nrow(x))
  
  paste(round(meanx * mult), rep(" (", ncol(x)), 
        round((meanx - 2 * sdx) * mult), rep(", ", ncol(x)),
        round((meanx + 2 * sdx) * mult), rep(")", ncol(x)), sep = "")
}


setEPS()
postscript("~/Dropbox/Apps/Overleaf/Confidence Sets for Causal Discovery/figures/whyAsympBad_max.eps", width = 6, height = 6)
par(mfrow = c(5,3), mar = c(1, 4, 4, 0))
for(z in 1:5){
  # hist(rec1[, z], freq = F, main = p, breaks = seq(0, 1, by = .05), )
  # mtext(round(mean(rec1[, z] < .05),3))
  hist(rec2[, z], freq = F, main = paste("Asymp: p=", p.list[z], sep = ""), breaks = seq(0, 1, by = .05), xlab = "")
  mtext(round(mean(rec2[, z] < .05),3))
  hist(rec3[, z], freq = F, main = paste("Oracle: p= ", p.list[z], sep = ""), breaks = seq(0, 1, by = .05), xlab = "")
  mtext(round(mean(rec3[, z] < .05),3))
  hist(rec4[, z], freq = F, main = paste("Proposed: p= ", p.list[z], sep = ""), breaks = seq(0, 1, by = .05), xlab = "")
  mtext(round(mean(rec4[, z] < .05),3))
}
dev.off()

tab <- rbind(getColumnPrint(rec2, .05, 1000),
      getColumnPrint(rec3, .05, 1000),
      getColumnPrint(rec4, .05, 1000))
rownames(tab) <- c("Asymp", "Oracle", "Proposed")
colnames(tab) <- p.list
xtable::xtable(tab)
