sim.size <- 5000
singleTest <- function(X, Y, K, bs = 200){
  n <- length(X)
  
  X.aug <- cbind(rep(1, length(X)),X)
  res <- RcppArmadillo::fastLm(X = X.aug, y = Y)$res
  

  Z <- t(X^K) %*% (diag(n) - X.aug %*% solve(t(X.aug) %*% X.aug) %*% t(X.aug)) / sqrt(n)
  testStat <- c(abs(t(X^K) %*% res) / sqrt(n))

  .helper <- function(){
    abs(Z %*% sample(res, replace = T))  
  }
  
  nullDist <- replicate(bs, .helper())
  
  return(mean(nullDist > testStat))    
}


run.once <- function(n, beta = .5){
  rec <- matrix(0, 1, 10)
  X <- rgamma(n, 1, 1) - 1
  Y <- beta * X + rgamma(n, 1, 1) - 1
  res1 <- RcppArmadillo::fastLm(X = cbind(rep(1, n), X), y = Y)$res
  res2 <- RcppArmadillo::fastLm(X = cbind(rep(1, n), Y), y = X)$res
  
  Xb <- rgamma(n, 1, 1) - 1
  Yb <- beta * Xb + rgamma(n, 1, 1) - 1
  res1b <- RcppArmadillo::fastLm(X = cbind(rep(1, n), Xb), y = Yb)$res
  res2b <- RcppArmadillo::fastLm(X = cbind(rep(1, n), Yb), y = Xb)$res
  
  
  rec[, 1:5] <- c(dHSIC::dhsic.test(X, res1)$p,
                TauStar::tauStarTest(X, res1, mode = "permutation")$p,
                t.test(X^2 * res1)$p.v,
                emplik::el.test(X^2 * res1, 0)$P,
                singleTest(X, Y, 2))
  
  rec[, 6:10] <- c(dHSIC::dhsic.test(Y, res2)$p,
                   TauStar::tauStarTest(Y, res2, mode = "permutation")$p,
                   t.test(Y^2 * res2)$p.v,
                   emplik::el.test(Y^2 * res2, 0)$P,
                   singleTest(Y, X, 2))
  
  colnames(rec) <- c("dhisc0", "tau0", "t_test0", "el0", "dkw0",
                     "dhisc1", "tau1", "t_test1", "el1", "dkw1")
  
  return(rec)  
}

cl <- parallel::makeCluster(15)
parallel::clusterExport(cl, varlist = ls())
rec <- parallel::parSapply(cl, 1:sim.size, function(i){run.once(100)})
rec <- t(rec)
colnames(rec) <- c("dhisc0", "tau0", "t_test0", "el0", "dkw0",
                   "dhisc1", "tau1", "t_test1", "el1", "dkw1")

write.csv(rec, "univariateEx_n100.csv")

rec <- parallel::parSapply(cl, 1:sim.size, function(i){run.once(1000)})
rec <- t(rec)
colnames(rec) <- c("dhisc0", "tau0", "t_test0", "el0", "dkw0",
                   "dhisc1", "tau1", "t_test1", "el1", "dkw1")
write.csv(rec, "univariateEx_n1000.csv")


run.onceSplit <- function(n, beta = .5){
  rec <- matrix(0, 1, 10)
  X <- rgamma(n, 1, 1) - 1
  Y <- beta * X + rgamma(n, 1, 1) - 1
  
  
  ## Used to estimate beta ##
  Xb <- rgamma(n, 1, 1) - 1
  Yb <- beta * Xb + rgamma(n, 1, 1) - 1
  betaHat1 <- RcppArmadillo::fastLm(X = cbind(rep(1, n), Xb), y = Yb)$coef
  betaHat2 <- RcppArmadillo::fastLm(X = cbind(rep(1, n), Yb), y = Xb)$coef
  
  
  
  res1 <- Y - X * betaHat1[2] - betaHat1[1]
  res2 <- X - Y * betaHat2[2] - betaHat2[1]
  
  
  rec[, 1:5] <- c(dHSIC::dhsic.test(X, res1)$p,
                  TauStar::tauStarTest(X, res1, mode = "permutation")$p,
                  t.test(X^2 * res1)$p.v,
                  emplik::el.test(X^2 * res1, 0)$P,
                  singleTest(X, Y, 2))
  
  rec[, 6:10] <- c(dHSIC::dhsic.test(Y, res2)$p,
                   TauStar::tauStarTest(Y, res2, mode = "permutation")$p,
                   t.test(Y^2 * res2)$p.v,
                   emplik::el.test(Y^2 * res2, 0)$P,
                   singleTest(Y, X, 2))
  
  colnames(rec) <- c("dhisc0", "tau0", "t_test0", "el0", "dkw0",
                     "dhisc1", "tau1", "t_test1", "el1", "dkw1")
  
  return(rec)  
}

cl <- parallel::makeCluster(15)
parallel::clusterExport(cl, varlist = ls())
rec <- parallel::parSapply(cl, 1:sim.size, function(i){run.onceSplit(100)})
rec <- t(rec)
colnames(rec) <- c("dhisc0", "tau0", "t_test0", "el0", "dkw0",
                   "dhisc1", "tau1", "t_test1", "el1", "dkw1")

write.csv(rec, "univariateEx_split_n100.csv")

rec <- parallel::parSapply(cl, 1:sim.size, function(i){run.onceSplit(1000)})
rec <- t(rec)
colnames(rec) <- c("dhisc0", "tau0", "t_test0", "el0", "dkw0",
                   "dhisc1", "tau1", "t_test1", "el1", "dkw1")
write.csv(rec, "univariateEx_split_n1000.csv")

















#######################################
# c1 <- rgb(100, 216, 255, max = 255, alpha = 60, names = "lt.blue")
# c2 <- rgb(255,192,203, max = 255, alpha = 60, names = "lt.pink")
# 
# 
# 
# # pdf("test/uni_n100.pdf", width = 5.5, height = 4)
# dat <- read.csv("test/results/univariateEx_n100.csv")[, -1]
# par(mfrow = c(2,3), mar = c(4, 4, 2, 2), oma = c(0, 0, 2, 0))
# 
# plotNames <- c("dHSIC", "Tau", "T Test", "Emp. Like.", "Resid. Boot.")
# 
# for(i in 1:5){
#   plot(-1, -1, xlim = c(0, 1), ylim = c(0, dim(dat)[1]), ylab = "frequency", xlab = "p-values")
#   hist(dat[,i], breaks = seq(0, 1, by = .05), add = T, col = c1)
#   hist(dat[,i + 5], breaks = seq(0, 1, by = .05), add = T, col = c2)
#   mtext(plotNames[i])
#   abline(h = dim(dat)[1] * .05, col = "black", lwd = 1)
# }
# 
# mtext("n = 100", outer = TRUE, cex = 1.5)
# # dev.off()
# # 
# # 
# # 
# # pdf("test/uni_n1000.pdf", width = 5.5, height = 4)
# dat1 <- read.csv("test/results/univariateEx_n1000.csv")[, -1]
# par(mfrow = c(2,3), mar = c(4, 4, 2, 2), oma = c(0, 0, 2, 0))
# 
# plotNames <- c("dHSIC", "Tau", "T Test", "Emp. Like.", "Resid. Boot.")
# 
# for(i in 1:5){
#   plot(-1, -1, xlim = c(0, 1), ylim = c(0, dim(dat)[1]), ylab = "frequency", xlab = "p-values")
#   hist(dat1[,i], breaks = seq(0, 1, by = .05), add = T, col = c1)
#   hist(dat1[,i + 5], breaks = seq(0, 1, by = .05), add = T, col = c2)
#   mtext(plotNames[i])
#   abline(h = dim(dat)[1] * .05, col = "black", lwd = 1)
# }
# 
# mtext("n = 1000", outer = TRUE, cex = 1.5)
# # dev.off()
# # 
# # 
# out <- cbind(t(rbind(colMeans(dat < .05)
# ,colMeans(dat < .1))),
# t(rbind(colMeans(dat1 < .05)
#         ,colMeans(dat1 < .1))))
# 
# out.final <- cbind(out[1:5, ], out[6:10,])
# colnames(out.final) <- 
# 
# xtable::xtable(out.final, digits = 3)


######################################
c1 <- rgb(100, 216, 255, max = 255, alpha = 60, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 60, names = "lt.pink")



# pdf("test/uni_n100.pdf", width = 5.5, height = 4)
dat <- read.csv("test/results/univariateEx_split_n100.csv")[, -1]
par(mfrow = c(2,3), mar = c(4, 4, 2, 2), oma = c(0, 0, 2, 0))

plotNames <- c("dHSIC", "Tau", "T Test", "Emp. Like.", "Resid. Boot.")

for(i in 1:5){
  plot(-1, -1, xlim = c(0, 1), ylim = c(0, dim(dat)[1]), ylab = "frequency", xlab = "p-values")
  hist(dat[,i], breaks = seq(0, 1, by = .05), add = T, col = c1)
  hist(dat[,i + 5], breaks = seq(0, 1, by = .05), add = T, col = c2)
  mtext(plotNames[i])
  abline(h = dim(dat)[1] * .05, col = "black", lwd = 1)
}

mtext("n = 100", outer = TRUE, cex = 1.5)
# dev.off()
#
#
#
# pdf("test/uni_n1000.pdf", width = 5.5, height = 4)
dat1 <- read.csv("test/results/univariateEx_split_n1000.csv")[, -1]
par(mfrow = c(2,3), mar = c(4, 4, 2, 2), oma = c(0, 0, 2, 0))

plotNames <- c("dHSIC", "Tau", "T Test", "Emp. Like.", "Resid. Boot.")

for(i in 1:5){
  plot(-1, -1, xlim = c(0, 1), ylim = c(0, dim(dat)[1]), ylab = "frequency", xlab = "p-values")
  hist(dat1[,i], breaks = seq(0, 1, by = .05), add = T, col = c1)
  hist(dat1[,i + 5], breaks = seq(0, 1, by = .05), add = T, col = c2)
  mtext(plotNames[i])
  abline(h = dim(dat)[1] * .05, col = "black", lwd = 1)
}

mtext("n = 1000", outer = TRUE, cex = 1.5)
# dev.off()
#
#
out <- cbind(t(rbind(colMeans(dat < .05)
,colMeans(dat < .1))),
t(rbind(colMeans(dat1 < .05)
        ,colMeans(dat1 < .1))))

out.final <- cbind(out[1:5, ], out[6:10,])


xtable::xtable(out.final, digits = 3)
