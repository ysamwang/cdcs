p <- 10
n.list <- c(30, 50, 100, 200, 500)
a.list <- c(1, 2, 4, 8)
param.grid <- expand.grid(n.list, a.list)
ss <- 1000

total.rec <- matrix(0, nrow(param.grid), 4)

for(j in 1:nrow(param.grid)){

  n <- param.grid[j, 1]
  a <- param.grid[j, 2]
  
  cat(n)
  cat(", ")
  cat(a)
  cat(" : ")
  rec <- matrix(0, ss, 4)
  
  for(i in 1:ss){
    beta <- rnorm(p+1)
    X <- cbind(rep(1, n), matrix(rgamma(n *p, 1, a) - 1/a, n, p))
    Y <- X %*% beta + (rgamma(n, 1, a) - 1/a)
    
    beta.hat <- RcppArmadillo::fastLm(X = X, y = Y)$coef
    
    errs <- Y - X %*% beta.hat
    errs <- (errs - mean(errs)) 
    
    proj <- diag(n) - X %*% solve(t(X) %*% X) %*% t(X)
    test.stat <- sqrt(sum((t(X^2) %*% errs / sqrt(n))^2))
    
    nullDist <- replicate(250, sqrt(sum((t(X^2) %*% proj %*% sample(errs, replace = T) / sqrt(n-p))^2)))
    
    rec[i,1] <- mean(test.stat <= nullDist)
    rec[i,2] <- IndepTest::MINTregression(X,Y,5,10,w=FALSE, rnorm(10000))
    rec[i,3] <- IndepTest::MINTregression(X,Y,5,10,w=FALSE, rnorm(10000, sd = sd(errs)))
    rec[i,4] <- IndepTest::MINTregression(X,Y,5,10,w=FALSE, sample(errs, 10000, replace = T))
    
    if(i %%1e2 == 0){
      cat(i)
      cat(" ")
    }
  }
  cat("\n")
  
  total.rec[j, ] <- colMeans(rec < .1)
  cat(round(total.rec[j, ], 3))
  cat("\n")
}



