pg <- expand.grid(n = c(100, 200, 500, 1000), p = c(10, 20, 50, 80))
totalRec <- matrix(0, dim(pg)[1], 4)
ss <- 1000
u <- 9
v <- 8


library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cl <- makeCluster(6) #not to overload your computer
registerDoParallel(cl)

for(j in 1:dim(pg)[1]){
  cat(j)
  cat(": ")
  n <- pg[j, 1]
  p <- pg[j, 2]
  
  rec <- foreach(i=1:ss, .combine = rbind) %dopar% {
    M <- toeplitz(.8^(0:(p-1)))
    M.sub <- M[-v, -v]
    X <- lcmix::rmvgamma(n, corr = M)
    S <- t(X[,-v]) %*% X[,-v]
    Z <- X[, u]^2 %*% X[,-v] %*% solve(S) %*% t(X[,-v]) %*% X[, v] / sqrt(n)
    Z1 <- X[, u]^2 %*% X[,-v] %*% solve(M.sub * n) %*% t(X[,-v]) %*% X[, v] / sqrt(n)
    
    c(Z, Z1)
  }
  totalRec[j, ] <- c(var(rec[,1]), mean(rec[,1]), mean(rec[,2]), mean((rec[,1] - mean(rec[,2]))^2))

}

total <- cbind(pg, totalRec, totalRec[, 1] * pg[,1], totalRec[, 4] * pg[,1]/log(pg[,2]) )



