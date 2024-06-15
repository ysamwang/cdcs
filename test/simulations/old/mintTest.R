p.list <- 10
a.list <- c(1)
b.list <- c(1)
n.list <- c(25, 50, 100, 200, 500)
param.grid <- expand.grid(a = a.list, b = b.list, n = n.list, p = p.list) #, a.list)
ss <- 1000

total.rec1 <- total.rec5 <- matrix(0, nrow(param.grid), 6)

library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(23) #not to overload your computer
registerDoParallel(cl)



for(j in 1:nrow(param.grid)){
  
  a <- param.grid[j, 1]
  b <- param.grid[j, 2]
  n <- param.grid[j, 3]
  p <- param.grid[j, 4]
  
  
  cat(paste(paste(c("a","b", "n", "p"), param.grid[j, ], sep = " = "), collapse = "; "))
  cat("\n")

  rec1 <- foreach(i=1:ss, .combine=rbind) %dopar%  {
    rec <- rep(0, 6)
    beta <- runif(p + 1, -2, 2)
    X <- cbind(rep(1, n), matrix(runif(n * p, -.5, .5), n, p))
    Y <- X %*% beta + (rnorm(n))
    # scale.X <- cbind(rep(1, n), X[,-1] * a)
    
    regOut <- RcppArmadillo::fastLm(X = X, y = Y)
    

    # proj <- diag(n) - X %*% solve(t(X) %*% X) %*% t(X)
    # test.stat <- sqrt(sum((t(X^2) %*% errs / sqrt(n))^2))
    # 
    # nullDist <- replicate(200, sqrt(sum((t(X^2) %*% proj %*% sample(errs, replace = T) / sqrt(n-p))^2)))
    

    rec[1] <- IndepTest::MINTregression(X,Y, 10, 20, w=FALSE, sample(regOut$res, n * 500, replace = T))
    rec[2] <- IndepTest::MINTregression(X,Y, 5, 10,w=FALSE, sample(regOut$res, n * 500, replace = T))
    rec[3] <- IndepTest::MINTregression(X,Y, 2, 5, w=FALSE, sample(regOut$res, n * 500, replace = T))
    # rec[4] <- IndepTest::MINTregression(X,Y, 5, 10,w=FALSE, sample(regOut$res, n * 500, replace = T))
    rec[4] <- IndepTest::MINTregression(X,Y, 10, 20, w=FALSE, rnorm( n * 500))
    rec[5] <- IndepTest::MINTregression(X,Y, 5, 10, w=FALSE, rnorm( n * 500))
    rec[6] <- IndepTest::MINTregression(X,Y, 2, 5, w=FALSE, rnorm( n * 500))
    
    
    # rec[4] <- IndepTest::MINTregression(scale.X,Y, 10, 20, w=FALSE, sample(regOut$res, n * 500, replace = T))
    # rec[5] <- IndepTest::MINTregression(scale.X,Y, 5, 10,w=FALSE, sample(regOut$res, n * 500, replace = T))
    # rec[6] <- IndepTest::MINTregression(scale.X,Y, 5, 10,w=FALSE, rnorm( n * 500))
    
    rec
  }
  total.rec1[j, ] <- colMeans(rec1 < .1)
  total.rec5[j, ] <- colMeans(rec1 < .05)
  cat(".1 :")
  cat(round(total.rec1[j, ], 3))
  cat("\n")
  cat(".05: ")
  cat(round(total.rec5[j, ], 3))
  
  cat("\n")
}

stopCluster(cl)

cbind(param.grid, total.rec1)

cbind(param.grid, total.rec5)
