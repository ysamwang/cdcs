#' Sample data from an additive noise model 
#'
#' @param p number of variables
#' @param n number of observations
#' @param parent_prob the probability of an edge between any two nodes
#' @param lowScale lower bound on variance of error terms
#' @param highScale upper bound on variance of error terms
#' @param lowEdge lower bound on edgeWeights
#' @param highEdge upper bound on edgeWeights 
#' @param edgeVar if null, defaults to uniform draws. If specified, this is the variance of the edgeweights which are drawn from gamma distribution
#' @param intercept draw random intercept or leave mean 0
#' @param dist the distribution of the error terms. Choices are: "gauss", "unif", "lognorm", "gamma", "weibull", "laplace", "mixed"
#' @param uniqueTop Whether to enforce a unique topological ordering so that u -> u+1 for all u
#' @param AdjMat an adjacency matrix which can be passed in to specify the edges instead of randomly selecting a graph
#' @param BInput a matrix of edges which can be passed in instead of randomly drawing edges
#' @param scalesInput a vector of error variances which can be passed in instead of randomly drawing
#' @param intInput a vector of intercepts which can be passed in instead of randomly drawing
#' @param posAndNeg whether edges should be both positive and negative or only positive
#' @return
#' #' \itemize{
#' \item B the matrix of edge weights
#' \item Y the n x p data
#' \item errs the error realizations
#' \item scales the variance of the errors
#' \item mu the intercept terms
#' }
#' 
rDAG <- function(p, n, parent_prob = 1/3, 
                 lowScale = 1, highScale = 1,
                 lowEdge = .3, highEdge = 1, edgeVar = NULL,
                 intercept = F, dist = "gauss", uniqueTop = T,
                 BInput = NULL, scalesInput = NULL, intInput = NULL, posAndNeg = T) {

  if(posAndNeg){
    signsAllowed <- c(-1, 1)
  } else {
    signsAllowed <- c(1)
  }
  
  if(is.null(BInput)){  
    ## Setup mat for directed edges
    B <- matrix(0, nrow = p, ncol = p)
    
    ## draw which edges will be active
    possibleEdges <- which(lower.tri(B, diag = F))
    actualEdges <- possibleEdges[which(rbinom(length(possibleEdges), size = 1, prob = parent_prob) == 1)]
    
    if(uniqueTop){
      actualEdges <- union(actualEdges, which(row(B) == (col(B) + 1))) 
    }
    
    ## fill in edge weights 
    if(is.null(edgeVar)){
      
      B[actualEdges] <- runif(length(actualEdges), lowEdge, highEdge)  * sample(signsAllowed, size = length(actualEdges), replace = T)
      
    } else {
      
      B[actualEdges] <-  sample(signsAllowed, size = length(actualEdges), replace = T) * rgamma(length(actualEdges), edgeVar, 1)
    }
    
  } else {
    
    B <- BInput
  }
  
  if(is.null(scalesInput)){
    
    scales <- runif(p, lowScale, highScale)
    
  } else {
    
    scales <- scalesInput
  }

  

  
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

    
  if(intercept){
    if(is.null(intInput)){
      mu <-  rnorm(p, 0, 1)  
    } else {
      mu <- intInput
    }
    
  } else {
    mu <-  rep(0, p)
  }
  
  
  Y <- mu +  solve(diag(rep(1, p)) - B, t(errs))
  Y <- t(Y)   

  
  return(list(B = B, Y = Y, errs = errs, scales = scales, mu = mu))
}
