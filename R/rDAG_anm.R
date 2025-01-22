#' Confidence sets for Causal Discovery
#'
#'
#' Sample data from an additive noise model 
#'
#' @param p number of variables
#' @param n number of observations
#' @param parent_prob the probability of an edge between any two nodes
#' @param lowScale lower bound on sd of error terms
#' @param highScale upper bound on sd of error terms
#' @param dist the distribution of the error terms. Choices are: "gauss", "unif", "lognorm", "gamma", "weibull", "laplace", "mixed"
#' @param uniqueTop Whether to enforce a unique topological ordering so that u -> u+1 for all u
#' @param AdjMat an adjacency matrix which can be passed in to specify the edges instead of randomly selecting a graph
#' @param funcType Either "GP", "poly" "fourier", "cam", or "tanh"
#' \itemize{
#' \item GP: sample from a gaussian process with parameters sigma and h
#' \item poly: sample from a random polynomial where the scores for the hermite polynomials have sd 1/k! 
#' \item fourier: sample from random fourier function where scores are 1/M^4
#' \item cam: function from Buhlmann et al 2013
#' \item tanh: tanh function
#' }
#' @param h bandwidth parameter if generating functions from Gaussian process "GP"
#' @param sigma parameter if generating functions from Gaussian process "GP"
#' @param M paramter
#' @param noParentMult multiplier for sd of errors if a node has no parent (i.e., is a root note)
#' @return
#' \itemize{
#' \item B the p x p adjacency matrix
#' \item Y the n x p data
#' \item errs the error realizations
#' \item scales the variance of the errors
#' } 
#' @export
#' 
rDAG_anm <- function(p, n, parent_prob = 1/3, lowScale = 1/5, highScale = sqrt(2)/5,
                     dist = "gauss", uniqueTop = T, AdjMat = NULL, funcType = "GP", 
                     h = 1, sigma = 1, M = 50, noParentMult = 5) {
  
  
  if(is.null(AdjMat)){  
    ## Setup mat for directed edges
    AdjMat <- matrix(0, nrow = p, ncol = p)
    
    ## draw which edges will be active
    possibleEdges <- which(lower.tri(AdjMat, diag = F))
    actualEdges <- possibleEdges[which(rbinom(length(possibleEdges), size = 1, prob = parent_prob) == 1)]
    
    if(uniqueTop){
      actualEdges <- union(actualEdges, which(row(AdjMat) == (col(AdjMat) + 1))) 
    }
    
    
    AdjMat[actualEdges] <- 1
    
  }
  
  noParents <- apply(AdjMat, MARGIN = 1, sum) == 0
  potentialScales <- runif(p, lowScale, highScale)
  scales <- ifelse(noParents, potentialScales * noParentMult, potentialScales)
  
  
  if(dist == "gauss") {
    
    errs <- matrix(rnorm(n * p) * scales, nrow = n, ncol = p, byrow = T)
    
  } else if(dist == "unif"){
    
    errs <- matrix(runif(n * p, -sqrt(3), sqrt(3))  * scales, nrow = n, ncol = p, byrow = T)
    
  } else if (dist == "lognorm") {
    
    errs <- matrix((exp(rnorm(n * p)) - exp(1/2)) / sqrt(exp(1) * (exp(1) - 1)) * scales, nrow = n, ncol = p , byrow = T)
    
  } else if (dist == "gamma"){
    
    errs <- matrix((rgamma(n * p, 1, 1) - 1) * scales, nrow = n, ncol = p, byrow = T)
    
  } else if (dist == "t"){
    
    errs <- matrix(rt(n * p, df = 10)  * scales, nrow = n, ncol = p, byrow = T)
    
  } else if (dist == "weibull"){
    a <- 3/4
    
    errs <- matrix( ((rweibull(n * p, shape = a, 1) - gamma(1 + 1/a)) / sqrt(gamma(1 + 2/a) - gamma(1 + 1/a)^2)) * scales,
                    nrow = n, ncol = p, byrow = T)
    
  } else if (dist == "laplace"){
    
    errs <- matrix(rmutil::rlaplace(n * p, 0, 1/sqrt(2)) * scales, nrow = n, ncol = p, byrow = T)
    
  } else if (dist == "mixed"){
    
    errs <- matrix(0, nrow = n, ncol = p)
    for(j in 1:p){
      
      distTemp <- sample(c("unif", "lognorm", "gamma", "weibull", "laplace"), size = 1)
      
      if(distTemp == "unif"){
        
        errs[,j] <- runif(n, -sqrt(3), sqrt(3))  * scales[j]
        
      } else if (distTemp == "lognorm") {
        
        errs[,j] <- (exp(rnorm(n)) - exp(1/2)) / sqrt(exp(1) * (exp(1) - 1)) * scales[j]
        
      } else if (distTemp == "gamma"){
        
        errs[,j] <- (rgamma(n , 1, 1) - 1) * scales[j]
        
      } else if (distTemp == "weibull"){
        a <- 3/4
        
        errs[,j] <- ((rweibull(n , shape = a, 1) - gamma(1 + 1/a)) / sqrt(gamma(1 + 2/a) - gamma(1 + 1/a)^2)) * scales[j]
      } else if (distTemp == "laplace"){
        
        errs[,j] <- rmutil::rlaplace(n , 0, 1/sqrt(2)) * scales[j]
      }
    }
  }
  
  Y <- matrix(0, n, p)
  
  Y[, 1] <- errs[, 1] 
  
  
    
    for(s in 2:(p)){
      for(t in 1:(s-1)){
        if(AdjMat[s, t] ==1){
          if(funcType == "GP"){
            
            f_st <- cdcs::sampleGP(Y[, t], h)
            
          } else if (funcType == "poly"){
            f_st <- cdcs::samplePoly(Y[, t], M)$f / sum(AdjMat[s,])
            
          } else if (funcType == "fourier"){
            f_st <- cdcs::sampleFourier(Y[, t], M)$f
            
          } else if (funcType == "cam"){
            f_st <- cdcs::sampleCF(Y[, t])$f
            
          } else if (funcType == "tanh"){
            f_st <- cdcs::sampleTanh(Y[, t])$f
            
          }
          
          
          
          Y[, s] <- Y[, s] + f_st
          
        } 
      }
      
      Y[, s] <- Y[, s] + errs[, s]
      
    }
    
  
  
  
  return(list(B = AdjMat, Y = Y, errs = errs, scales = scales))
}



sampleGP <- function(X, h = 1, sigma = 1){
  cov_fun <- function(X, h = h){
    outer(X, X, function(a, b){sigma^2 * exp(-(a - b)^2/(2 * h^2))})
  }
  
  xEval <- c(mvtnorm::rmvnorm(1, sigma = cov_fun(X, h = h)))  
  return(xEval)
}


# Evaluate the first M hermite polynomials on X
# Normalized to have variance 1 under standard normal
getHermite <- function(X, M){
  mat <- matrix(0, length(X), M + 1)
  
  mat[, 1] <- 1
  mat[, 2] <- X
  for(i in 2:M){
    mat[,i + 1] <- X * mat[, i] - (i-1) * mat[, i-1]
  }
  
  mat <- mat[,-1]
  mat <- scale(mat)
  return(mat)
}

# Sample a polynomial from the hermite polynomials
# where sd of scores 1/k! 
samplePoly <- function(X, M){
  mat <- cdcs::getHermite(X, M)
  scores <- rnorm(M, mean = 0, sd = 1) / factorial(1:M)
  # scores <- rnorm(M, mean = 0, sd = 1)
  funcEval <- mat %*% scores
  return(list(f = c(funcEval), scores = scores))
}

getFourier <- function(X, M){
  mat <- sapply(1:M,FUN = function(s){ifelse(rep(rbinom(1, size = 1, prob = 1/2), length(X)),
                                             sin(X * s/2 ), cos(X * s/2 ))} )
  return(mat) 
}

sampleFourier <- function(X, M){
  mat <- cdcs::getFourier(X, M)
  scores <- rnorm(M, mean = 0, sd = 1/(1:M)^2)
  funcEval <- mat %*% scores
  return(list(f = c(funcEval), scores = scores))
}



sampleCF <- function(X){
  
  a <- rexp(1, 4) + 1
  b <- sample(c(-1, 1), 1) * runif(1, .5, 2)
  c <- runif(1, -2, 2)
  funcEval <- a * b*(X + c) / (1 + abs(b*(X + c)))
  
  return(list(f = c(funcEval), param = c(a, b, c)))
}

sampleTanh <- function(X){
  
  a <- rexp(1, 4) + 1
  b <- sample(c(-1, 1), 1) * runif(1, .5, 2)
  c <- runif(1, -2, 2)
  funcEval <- a * tanh(X * b + c)
  
  return(list(f = c(funcEval), param = c(a, b, c)))
}

