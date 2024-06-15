p <- 8
n <- 1e5
dat <- cdcs::rDAG(p, n, parent_prob = 1, dist = "gamma", intercept = F)
B <- dat$B
mu <- dat$mu
Y.mu <- solve(diag(p) - B, mu)

v <- 3
U <- setdiff(c(1:6), 3)

Y.mod <- cbind(rep(1,n), dat$Y)


sigma <- rbind(c(1, Y.mu),
               cbind(Y.mu, solve(diag(p) - B) %*% (mu %*% t(mu) + diag(p)) %*%  t(solve(diag(p) - B))))

sig.hat <- t(Y.mod) %*% Y.mod / n
sigma - sig.hat


delta <- round(solve(sigma[c(1, U + 1), c(1, U + 1)]) %*% sigma[c(1, U + 1),v+1], 8)



n.list <- c(50, 100, 250, 500, 1000)
ss <- 5e4
bs.draws <- 400

rec <- rec1 <- rec2 <- rep(0, length(n.list))
for(j in 1:length(n.list)){
  n <- n.list[j]
  cat("N = ")
  cat(n)
  cat("\n")
  
  tempRec <- tempRec1 <- tempRec2 <- rep(0, ss)
  for(i in 1:ss){
    dat <- cdcs::rDAG(p, n, parent_prob = 1/3, BInput = B, intInput = mu, dist = "gamma", intercept = F)
    dat$Y <- dat$Y
    X <- cbind(rep(1, n), dat$Y[, U, drop = F])
    y <- dat$Y[, v, drop=F]
  
    res <- RcppArmadillo::fastLm(X = X, y = y)$res

    resTrue <- y - X %*% delta
    
    tempRec[i] <- max(abs(res - resTrue))
    
    Z <- t(X^2)
    
    tempRec1[i] <- sqrt(sum((Z %*% res / sqrt(n))^2)) / length(U)
    

    tempRec2[i] <- mean(replicate(bs.draws, sqrt(sum((Z %*% sample(res, replace = T) / sqrt(n))^2 )) / length(U)))
    # tempRec1[i] <- mean((Z %*% sample(res, replace = T) / sqrt(n))^2)
    if(i %%1e4 == 0){
      cat(i)
      cat(" ")
    }
  }
  rec[j] <- mean(tempRec)
  rec1[j] <- var(tempRec1)
  rec2[j] <- mean(tempRec2)
  cat("\n")
  cat("Res: ")
  cat(round(rec2[j],3))
  cat("\n")
}

round(rec1 / log(n.list), 3)
round(rec2 / log(n.list), 3)
