directLiNGAM <- function(Y, verbose = T, metric = "dsic"){
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  
  
  nodesLeft <- 1:p
  ordering <- c()
  
  
  .getKernelMeasure <- function(j){
    nodesLeftMinusJ <- setdiff(1:p, c(ordering, j))
    
    residMat <- sapply(nodesLeftMinusJ, function(i){RcppArmadillo::fastLm(X = cbind(rep(1, n), Y[, j]), y = Y[, i])$res})
    
    if(metric == "dhsic"){
      
      ret <- sum(apply(residMat, MAR = 2, function(y){dhsicMod(Y[, j], y,
                                                               resample = (n > 500),
                                                               sampleSize = 500,
                                                               numResamples =  10)}))
      
    } else {
      
      ret <- sum(apply(residMat, MAR = 2, function(y){TauStar::tStar(x = Y[, j], y,                                                                     resample = (n > 1000),
                                                                     sampleSize = 1000,
                                                                     numResamples =  10)}))
    }
      
    return(ret)
  }
  
  while(length(ordering) < p-1){
    ## Get measures of independence of residual vs regressor
    measures <- sapply(nodesLeft, .getKernelMeasure)
    names(measures) <- nodesLeft
    
    ## Select a root
    root <- nodesLeft[which.min(measures)]
    
    
    ## Print details
    if(verbose){
      cat(paste("=== Step: ", length(ordering) + 1 ," ===", sep = ""))
      cat("\n")
      cat("Order so far: ")
      cat(ordering)
      cat("\n")
      print(measures)
      cat("Selected Root: ")
      cat(root)
      cat("\n\n")
    }
    
    ## Update nodesLeft and ordering
    nodesLeft <- setdiff(nodesLeft, root)
    ordering <- c(ordering, root)
    
    # regress all nodes left onto the root
    Y[, nodesLeft] <- sapply(nodesLeft, function(i){RcppArmadillo::fastLm(X = cbind(rep(1, n), Y[, root,  drop = F]), y = Y[, i,  drop = F])$res})
    
    
  }
  ordering <- c(ordering, nodesLeft)
  
  return(ordering)
}


dhsicMod <- function(x, y, resample = T, sampleSize = 500, numResamples = 10){
  
  .resampleOnce <- function(){
    ind <- sample(length(x), size = sampleSize)
    dHSIC::dhsic(x[ind], y[ind])$dHSIC
  }
    
  if(resample){
      ret <- mean(replicate(numResamples, .resampleOnce()))
  } else {
    ret <- dHSIC::dhsic(x,y)$dHSIC
  }
  
  return(ret)
  
}






