#' Confidence sets for Causal Discovery
#'
#'
#' Branch and bound for an exhaustive search
#'
#' @param Y an n by p matrix with the observations
#' @param G a 3d array containing the evaluated test functions. G[i, j, u] is \eqn{h_j(Y_{u,i})} 
#' @param bs the number of bootstrap resamples for the null distribution
#' @param aggType the aggregation used for the test statistic
#' #' \itemize{
#' \item "1" \eqn{\sum_{u \prec pr(v), j} \vert h_j(Y_u)^T \eta_v \vert}
#' \item "2" \eqn{\sum_{u \prec pr(v), j} \vert h_j(Y_u)^T \eta_v \vert^2}
#' \item "3" or "inf": \eqn{\max_{u \prec pr(v), j} \vert h_j(Y_u)^T \eta_v \vert} 
#' }
#' @param pValueAgg the procedure for aggregating p-values 
#' #' \itemize{
#' \item "tippet": \eqn{T = \min p_{\theta(v)})} which follows a Beta(1, p-1) distribution
#' \item "fisher": \eqn{T = -2 \sum_{\theta(v) > 1} \log( p_{\theta(v)})} which follows a chi-squared with p-1 degrees of freedom
#' }
#' @param intercept 1 indicates that an intercept should be used; 0 indicates no intercept is used
#' @param verbose T or F. Indicates whether to print updates to console
#' @param alpha the size of the test
#' @return
#' a vector of p-values corresponding to test 
#' @export
brandAndBound <- function(Y, G, bs = 200,
                          aggType = 3, alpha = .05,
                          pValueAgg = "tippett", intercept = 1, verbose = T){
  
  
  ## Dimensions of Y
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  
  
  ## Check if aggType is correct ##
  if(!(aggType %in% c(1, 2, 3, "inf"))){
    stop("Incorrect aggType: must be either '1', '2', '3', 'inf'")
  }
  
  ## Check if intercept is correct ##
  if(!(intercept %in% c(T, F))){
    stop("intercept must be either TRUE or FALSE")
  }
  
  
  ### Change aggType and intercept to numeric for passing onto helpers   
  if(aggType == "inf"){
    aggType <- 3
  }
  
  # make intercept 0/1
  intercept <- as.numeric(intercept)
  

  
  ### Helper function to update the counter for the p-values 
  .updatePvals <- function(nativeInd){
    
    ### Get the set and trackers from the larger function environment 
    matchInd <- hashInd[nativeInd]
    currentSeq_i <- unlist(currentSeq[nativeInd, -1])
    currentTrack <- currentSeq[nativeInd, 1]
    
    ## get P-values computed from function environment ##
    freshPVals <- unlist(uniqueRes[[matchInd]]["pVals"])
    possibleChildren <- unlist(uniqueRes[[matchInd]]["possibleChildren"])
    
    
    # update orderings using either fisher or tippet
    if(pValueAgg == "fisher"){
      
      newTrack <- currentTrack - 2 * log(freshPVals)
      continueOn <- newTrack < cutoff
      
    } else {
      
      newTrack <- pmin(currentTrack, freshPVals)
      continueOn <- newTrack > cutoff  
      
    }
    
    
    
    # Check if there are any orderings which haven't passed the cut-off 
    # i.e., still viable orderings
    if(sum(continueOn) > 0){
      
      # data frame to return
      # Col 1: updated tracker
      # Remaining columns: Orderings which have not passed the cut-off
      newSeq <- data.frame(newTrack[which(continueOn)],
                           matrix(currentSeq_i, nrow = sum(continueOn), ncol = length(currentSeq_i), byrow = T),
                           possibleChildren[which(continueOn)])
      
    } else {
      
      # if no orderings have passed the cut-off, return an empty data frame    
      newSeq <- data.frame(matrix(0, nrow = 0, ncol = length(currentSeq_i) + 2))
      
    }
    
    # names for each column
    names(newSeq) <- c("currentTrack", paste("V", 1:(length(currentSeq_i) + 1), sep = ""))
    rownames(newSeq) <- NULL
    
    return(newSeq)
  }
  
  
  
  ### takes in a set of ancestors and tests whether any node not included in the set 
  ### could be a descendant of the ancestors
  .testAncest <- function(ancest){
    
    ## will come in as a list because gets pulled out of data.frame
    ancest <- unlist(ancest)
    
    ## set of any nodes not in the set which could be potential children
    possibleChildren <- setdiff(1:p, ancest)
    
    ## subtract 1 for cpp indexing
    pVals <- cdcs::bnbHelperGof(ancest - 1, possibleChildren - 1,
                                Y, G, withinAgg = aggType, aggType = aggType, 
                                bs = bs, intercept = intercept)
    
    ## data frame with:
    ## Column 1: p-values (converted to uniform)
    ## Column 2: possible children
    ret <- data.frame(pVals =  cdcs::convToCont(pVals$pVals, bs), possibleChildren = possibleChildren)
    return(ret)
    
  }
  
  
  ### Calculates a cut-off based on whether the aggregation type is a max or sum
  if(pValueAgg == "fisher"){
    
    ### Cut off using fishers aggrgation method ###
    cutoff <- qchisq(1 - alpha, df = 2 * (p-1))
    currentTrack <- rep(0, p)
    
  } else {
    
    ### Calculate the quantile of the minimum of uniforms ###
    cutoff <- qbeta(alpha, 1, p-1)
    currentTrack <- rep(1, p)
    
  }
  
  
  # Initial values
  # Start with each node and initial currentTrack value
  currentSeq <- data.frame(currentTrack, matrix(1:p, ncol = 1))
  names(currentSeq) <- c("currentTrack", "V1")
  
  
  while(dim(currentSeq)[1] > 0 & dim(currentSeq)[2] <= p){
    
    # hash involves all sequences which are the same set
    # computation only depends on v and set an(v), so the "ordering" of an(v)  doesn't matter
    hash <- apply(currentSeq[, -1, drop = F], MARGIN = 1, function(x){paste(sort(unlist(x)), collapse = ".")})
    
    # uniqueHash is the list of unique ancestral sets
    uniqueHash <- unique(hash)
    
    # each ongoing ordering maps to a specific set
    hashInd <- match(hash, uniqueHash)
    
    
    if(verbose){
      print("==========================================================")
      print(paste("Current Order: ", dim(currentSeq)[2] -1, sep = ""))
      print(paste("Number of Perm: ", dim(currentSeq)[1] , "( ",
                  round(dim(currentSeq)[1] / prod(c(p:(p-(dim(currentSeq)[2] -2)))), 3) ," )", sep = ""))
      print(paste("Number of Comb: ", length(uniqueHash), sep = "" ))
      cat("\n")
    }  
    
    
    # For each value in uniqueHash, get a representative set and run .testAncest on that set
    uniqueRes <- apply(currentSeq[match(uniqueHash, hash), -1, drop = F], MARGIN = 1, .testAncest)
    
    ### Update Sequences ###
    updatedSeq <- lapply(1:length(hash), .updatePvals)
    currentSeq <- as.data.frame(data.table::rbindlist(updatedSeq))

  }
  
  ## Check if there are any orderings to return
  if(dim(currentSeq)[1] == 0){
    
    # return empty data frame
    currentSeq <- data.frame(matrix(0, nrow = 0, ncol = p + 1))
  
    } else {
  
      ## if there are orderings to be returned
      ## Col 1: final p-value of aggregated p-values 
      ## then form data frame with orderings
      if(pValueAgg == "fisher"){
        
        pVals <- pchisq(currentSeq[, 1], df = 2 * (p-1), lower.tail = F)
      
        } else {
        
        pVals <- pbeta(currentSeq[, 1], 1, p-1, lower.tail = T)
        
      }
      
      currentSeq[, 1] <- pVals
  }
  names(currentSeq) <- c("pValue", paste("V", 1:p, sep = ""))
  
  return(currentSeq)
}


