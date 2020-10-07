rDAG <- function(p, n, parent_prob, lowScale = 1, highScale = 1, lowEdge = .3, highEdge = 1, edgeVar = NULL,
                 dist = "gauss", uniqueTop = T, BInput = NULL, scalesInput = NULL) {

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
      
      B[actualEdges] <- runif(length(actualEdges), lowEdge, highEdge)  * sample(c(-1, 1), size = length(actualEdges), replace = T)
      
    } else {
      
      B[actualEdges] <-  sample(c(-1, 1), size = length(actualEdges), replace = T) * rgamma(length(actualEdges), edgeVar, 1)
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
    
  } else if (dist == "weibull"){
    a <- 3/4
    
    errs <- matrix( ((rweibull(n * p, shape = a, 1) - gamma(1 + 1/a)) / sqrt(gamma(1 + 2/a) - gamma(1 + 1/a)^2)) * scales,
                    nrow = n, ncol = p, byrow = T)
  } else if (dist == "laplace"){
    
    errs <- matrix(rmutil::rlaplace(n * p, 0, 1/sqrt(2)) * scales, nrow = n, ncol = p, byrow = T)
    
  }
  
  Y <- solve(diag(rep(1, p)) - B, t(errs))
  Y <- t(Y)
  
  return(list(B = B, Y = Y, errs = errs, scales = scales))
}
