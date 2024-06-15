n.list <- c(100, 200, 500, 800, 1000)
p.list <- c(10)
pg <- expand.grid(n = n.list, p = p.list)
totalRec <- matrix(0, dim(pg)[1], 21)
ss <- 10000
u <- 9
v <- 8


library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cl <- makeCluster(2) #not to overload your computer
registerDoParallel(cl)

for(j in 1:dim(pg)[1]){
  cat(j)
  cat(": ")
  n <- pg[j, 1]
  p <- pg[j, 2]
  
  rec <- foreach(i=1:ss, .combine = rbind) %dopar% {
    
    M <- toeplitz(.8^(0:(p-1)))
    M.sub <- M[-v, -v]
    X <- lcmix::rmvgamma(n, corr = M) - 1
    S <- t(X[,-v]) %*% X[,-v]
    s.half <- expm::sqrtm(solve(S))
    
     
    Z <- X[, u]^2 %*% X[,-v] %*% solve(S) %*% t(X[,-v]) %*% X[, v] / sqrt(n)
    Z1 <- s.half %*% t(X[,-v]) %*% X[, v]
    Z2 <- t(X[, u]^2) %*% X[,-v] %*% s.half / sqrt(n)
    
    c(Z, Z1, Z2)
  }
  
  
  totalRec[j, ] <- c(var(rec[, 1]), mean(rec[, 1]), sum(colMeans(rec[, 2:10]) * colMeans(rec[, 11:19])),
                     apply(rec[, -1], MAR = 2, var))

  
}


(totalRec[,2] - totalRec[, 3]) * sqrt(n.list) 




# total <- cbind(pg, totalRec, totalRec[, 1] * pg[,1], totalRec[, 4] * pg[,1]/log(pg[,2]) )



