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
#' }
#' @param K degree of moment for moment based statistic
#' @return
#' a p-value testing the hypothesis 
senSen2014 <- function(X, Y, bs = 200, statType = "hsic", K = 3){

  if(!(statType %in% c("hsic", "moment"))){
    stop("Invalid statType")
  }
  
  n <- dim(X)[1]
  p <- dim(X)[2]

  
  
  ## Regress Y onto X
  regOutput <- RcppArmadillo::fastLm(X = X, y = Y)
  # center the errors

  errs <- regOutput$res - mean(regOutput$res)
  
  ## calculate the test statistic
   if(statType == "hsic"){
     testStat <- n * dHSIC::dhsic(X, regOutput$res)$dHSIC
     
   } else {
     
     testStat <- c(mean(abs(t(X^K[1]) %*% regOutput$res / sqrt(n))), mean((t(X^K[1]) %*% regOutput$res / sqrt(n))^2 ))
     
     if(length(K) > 1){
       
       for(k in 2:length(K)){
         
         testStat <- testStat + c(mean(abs(t(X^K[k]) %*% regOutput$res / sqrt(n))),
                                  mean((t(X^K[k]) %*% regOutput$res / sqrt(n))^2 ))
       }
       
     }
   }
  
    

                
  

  
  .drawBS <- function(){
    
    # draw predictors with replacement
    X.bs <- X[sample(n, replace = T), , drop = F]
    # form pseudo data
    yCheck <- X.bs %*% regOutput$coeff + sample(errs, replace = T)
    # regress pseudo data
    regOutputBS <- RcppArmadillo::fastLm(X = X.bs, y = yCheck)
    
    # calculate test stat
    if(statType == "hsic"){
      
      testStatBS <- n * dHSIC::dhsic(X.bs, regOutputBS$res)$dHSIC
      
    } else {
      
      #calculate tau's and scale by n-p
      intermed <- abs(t(X.bs^K[1]) %*% regOutputBS$res / sqrt(n-p))
      # calc test stat
      testStatBS <- c(mean(intermed), mean(intermed^2))
      
      # if more than 1 value for K is being used
      if(length(K) > 1){
        for(k in 2:length(K)){

          intermed <- abs(t(X.bs^K[k]) %*% regOutputBS$res / sqrt(n-p))
          testStatBS <- testStatBS + c(mean(intermed), mean(intermed^2))
          
        }
      }
    }

    return(testStatBS)
  }  
  
  nullDist <- replicate(bs, .drawBS())  
  
  if(statType == "hsic"){
    return(mean(nullDist >= testStat))
    
  } else {
    
    return(c(mean(nullDist[1, ] >= testStat[1]),
           mean(nullDist[2, ] >= testStat[2])))
    
  }
  
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
#' }
#' @param K degree of moment for moment based statistic
#' @return
#' a p-value testing the hypothesis 
dkw2020 <- function(X, Y, bs = 200, statType = "moment", K = c(2, 3), aggType = 3){
  
  n <- dim(X)[1]
  p <- dim(X)[2]

  
  if(statType == "moment"){
  
    out <- singleTestcppMulti(dat = X, Y = Y, K = K, aggType = 3, bs = bs, intercept = 0)
    
    return(c(mean(out$nullDistOne >= out$testStatOne),
             mean(out$nullDistTwo >= out$testStatTwo)))
    
  } else {
    
    
    hatMat <- diag(n) - X %*% solve(t(X) %*% X, t(X))
    errs <- hatMat %*% Y
    

    testStat <- n * dHSIC::dhsic(X, errs)$dHSIC
    
    
    errs <- errs - mean(errs)
    
    .drawBS <- function(){
      err.BS <- hatMat %*% sample(errs, replace = T)
      
      testStatBS <- n * dHSIC::dhsic(X, err.BS)$dHSIC
      
      return(testStatBS)
    }  
    
    nullDist <- replicate(bs, .drawBS())  
    return(mean(nullDist >= testStat)) 
  }
  

}





### From Sen and Sen (2014)
### Code written by Sen and Sen
### http://www.stat.columbia.edu/~bodhi/research/HSIC-IndependenceTest-Predictor&Error.zip
# HSIC_senSen <- function(X,Y,hx = 1,hy = 1){
#   
#   n <- length(Y); 
#   
#   K <- matrix(0,n,n);
#   L <- matrix(0,n,n);
#   H <- diag(n) - matrix(1,n,n)/n;
#   
#   J <- matrix(1,n,1); 
#   # K = exp(-(X%*%t(J) - J%*%t(X))^2/hx);
#   L = exp(-(Y%*%t(J) - J%*%t(Y))^2/hy);
#   
#   K = exp(-as.matrix(dist(X, method = "euclidean", diag = TRUE, upper = TRUE))^2/hx);
#   
#   # for (i in 1:n){
#   #    for (j in 1:n){
#   #        K[i,j] <- exp(-norm(as.matrix(X[i,] - X[j,]))^2/hx);
#   #        # L[i,j] <- exp(-norm(as.matrix(Y[i,] - Y[j,]))^2/hy);
#   #    }
#   #}
#   
#   testStat = sum(diag(K%*%H%*%L%*%H))/n;
#   return(testStat);
# }
# 
