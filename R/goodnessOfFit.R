#' Confidence sets for Causal Discovery
#'
#'
#' Goodness of fit test from Sen and Sen (2014)
#'
#' @param Y an n vector with the observed outcome variables
#' @param X an n by p matrix with the covariates
#' @param bs the number of bootstrap resamples for the null distribution
#' @param statType 
#' #' \itemize{
#' \item "hsic": uses dHSIC
#' \item "moment": uses moment based statistic with degree K
#' \item "dcor": uses distance correlation 
#' }
#' @param K degree of moment for moment based statistic
#' @sampleSplit Split the sample using the second half to estimate beta and first half to run independence test
#' @return
#' a p-value testing the hypothesis 
senSen2014 <- function(X, Y, bs = 200, statType = "hsic", K = c(2,3), cutoff = NULL,
                       sampleSplit = F, kernel = "gaussian"){

  if(!(statType %in% c("hsic", "moment", "dcor"))){
    stop("Invalid statType")
  }
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  if(is.null(cutoff)){cutoff <- log(n) * 3}

  
  if(sampleSplit){
    ### If sample splitting ###
    n <- floor(n/2)
    betaHat <- RcppArmadillo::fastLm(X = X[(n+1):nrow(X), ], y = Y)
    X <- X[1:n, ]
    Y <- Y[1:n, ]
    errs.uncentered <- Y - X %*% betaHat
    errs <- errs.uncentered - mean(errs.uncentered)
    
  } else{
    
    regOutput <- RcppArmadillo::fastLm(X = X, y = Y)
    # center the errors
    betaHat <- regOutput$coeff
    errs.uncentered <- regOutput$res
    errs <- errs.uncentered - mean(errs.uncentered)
    
  } 


  
  ## calculate the test statistic
   if(statType == "hsic"){
     
     testStat <- n * dHSIC::dhsic(X, errs, kernel = kernel)$dHSIC
     
    } else if (statType == "dcor"){
       testStat <- n * energy::dcor(X, errs)
    
    } else {
     
      if(is.null(cutoff)){cutoff <- log(n) * 3}
      
      G <- array(0, dim = c(n, length(K), p))
      for(k in 1:length(K)){
        G[, k, ] <- sign(X^K[k]) * pmin(abs(X^K[k]), cutoff)
      }

      testStat <- max(sapply(1:p,
                         function(pp){ max(abs(t(G[, , pp]) %*% errs / sqrt(n)))}))
   }
  
    

                
  

  
  .drawBS <- function(){
    
    # draw predictors with replacement
    X.bs <- X[sample(n, replace = T), , drop = F]
    # form pseudo data
    yCheck <- X.bs %*% betaHat + sample(errs, replace = T)
    # regress pseudo data
    regOutputBS <- RcppArmadillo::fastLm(X = X.bs, y = yCheck)
    
    # calculate test stat
    if(statType == "hsic"){
      
      testStatBS <- n * dHSIC::dhsic(X.bs, regOutputBS$res, kernel = kernel)$dHSIC
      
      } else if (statType == "dcor"){
      
        testStatBS <- n * energy::dcor(X, regOutputBS$res)
    
      } else {
      
        G <- array(0, dim = c(n, length(K), p))
        for(k in 1:length(K)){
          G[, k, ] <- sign(X.bs^K[k]) * pmin(abs(X.bs^K[k]), cutoff)
        }
        
        testStatBS <- max(sapply(1:p, function(pp){
                     max(abs(t(G[, , pp]) %*% regOutputBS$res / sqrt(n - p))) }))
    }

    return(testStatBS)
  }  
  
  nullDist <- replicate(bs, .drawBS())  
  

  return((sum(nullDist >= testStat) + 1) / (bs + 1))
  
}


#' Confidence sets for Causal Discovery
#'
#'
#' Goodness of fit test
#'
#' @param Y an n vector with the observed outcome variables
#' @param X an n by p matrix with the covariates
#' @param bs the number of bootstrap resamples for the null distribution
#' @param statType 
#' #' \itemize{
#' \item "hsic": uses dHSIC
#' \item "moment": uses moment based statistic with degree K
#' \item "dcor": uses distance correlation 
#' }
#' @param K degree of moment for moment based statistic
#' @sampleSplit Split the sample using the second half to estimate beta and first half to run independence test
#' @return
#' a p-value testing the hypothesis 
dkw2020 <- function(X, Y, bs = 200, statType = "moment", K = c(2, 3),
                    withinAgg = 3, sampleSplit = F, kernel = "gaussian", cutoff = NULL){
  
  

  n <- dim(X)[1]
  #if(sd(X[,1]) >1e-8){X <- cbind(rep(1, n), X) }
  p <- dim(X)[2]
  if(is.null(cutoff)){cutoff <- log(n) * 3}

  
  if(statType == "moment"){
    
    G <- array(0, dim = c(n, length(K), p))
    for(k in 1:length(K)){
      G[, k, ] <- sign(X^K[k]) * pmin(abs(X^K[k]), cutoff)
    }
    
    out <- gofTest(covariates = X, Y = Y,
            G = G, bs = bs, withinAgg = withinAgg, intercept = 0,
            sampleSplit = sampleSplit)
    
    return(c(
      (sum(out$nullDistOne > out$testStatOne) + 1)/ (bs + 1),
      (sum(out$nullDistTwo > out$testStatTwo) + 1)/ (bs + 1),
             (sum(out$nullDistTwo > out$testStatTwo) +1 ) / (bs + 1)
             )
    )
    
  } else {
    
    
    
    if(sampleSplit){
      ### If sample splitting ###
      n <- floor(n/2)
      betaHat <- RcppArmadillo::fastLm(X = X[(n+1):nrow(X), ], y = Y)
      X <- X[1:n, ]
      Y <- Y[1:n, ]
      hatMat <- diag(n) - X %*% solve(t(X) %*% X, t(X))
      errs <- Y - X %*% betaHat
      
      
    } else{

      hatMat <- diag(n) - X %*% solve(t(X) %*% X, t(X))
      errs <- hatMat %*% Y
      
    } 
    
    if(statType == "hsic"){
      testStat <- n * dHSIC::dhsic(X, errs, kernel = kernel)$dHSIC
    } else if (statType == "dcor"){
      testStat <- n * energy::dcor(X, errs)
    }

    
    
    errs <- errs - mean(errs)
    
    .drawBS <- function(){
      err.BS <- hatMat %*% sample(errs, replace = T)
      
      if(statType == "hsic"){
        testStatBS <- n * dHSIC::dhsic(X, err.BS, kernel = kernel)$dHSIC
      } else if (statType == "dcor"){
        testStatBS <- n * energy::dcor(X, err.BS)
      }
      

      
      return(testStatBS)
    }  
    
    nullDist <- replicate(bs, .drawBS())  
    return((sum(nullDist >= testStat) + 1) / (bs + 1)) 
  }
  
  
}

  
  
  

senSenWrapper <- function(x, X, y,Boots = 200){
  
  out <- RcppArmadillo::fastLm(X = X, y = y)
  
  
  bs <- TestIndBoots(x,X,out$res,b = out$coef,hx =1,hy = 1,Boots = Boots)
  
  return(bs[2])
}




TestIndBoots <- function(x,X,e,b,hx,hy,Boots = 100){
  
  T_hat <- HSIC(x,e,2*hx*hx,2*hy*hy);   ## Compute the test-statistic
  # plot(T_hat,0,cex = 1.0, col = "dark red");
  # Implementing the bootstrap procedure
  
  T_hat_B <- matrix(0,Boots,1);
  e_0 <- e - mean(e);				## centered residuals
  
  for(j in 1:Boots){
    idx <- sample(n,n,replace=T);	## with replacement samples from the errors
    e_B <- e_0[idx];
    
    idx2 <- sample(n,n,replace=T);	## with replacement samples from the predictors
    x_B <- x[idx2,];
    X_B <- X[idx2,];
    
    
    yhat_B <- X_B%*%b + e_B;			## Create the bootstrap response values
    bhat_B <- solve(t(X_B)%*%X_B)%*%t(X_B)%*%yhat_B;
    e_hat_B <- yhat_B - X_B%*%bhat_B;		## Bootstrap residuals
    
    
    hx_B <- hx; hy_B <- hy;
    T_hat_B[j] <- HSIC(x_B,e_hat_B,2*hx_B*hx_B,2*hy_B*hy_B);
    
  }
  pval <- mean(T_hat_B >= T_hat);
  # points(T_hat_B,matrix(0,Boots,1),cex = 0.3);
  return(cbind(T_hat,pval));
  
}


HSIC <- function(X,Y,hx,hy){
  
  n <- length(Y); 
  
  K <- matrix(0,n,n);
  L <- matrix(0,n,n);
  H <- diag(n) - matrix(1,n,n)/n;
  
  J <- matrix(1,n,1); 
  # K = exp(-(X%*%t(J) - J%*%t(X))^2/hx);
  L = exp(-(Y%*%t(J) - J%*%t(Y))^2/hy);
  
  K = exp(-as.matrix(dist(X, method = "euclidean", diag = TRUE, upper = TRUE))^2/hx);
  
  # for (i in 1:n){
  #    for (j in 1:n){
  #        K[i,j] <- exp(-norm(as.matrix(X[i,] - X[j,]))^2/hx);
  #        # L[i,j] <- exp(-norm(as.matrix(Y[i,] - Y[j,]))^2/hy);
  #    }
  #}
  
  testStat = sum(diag(K%*%H%*%L%*%H))/n;
  return(testStat);
}


