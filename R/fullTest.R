#' Confidence sets for Causal Discovery
#'
#'
#' Test whether a specific ordering is consistent with the data
#'
#' @param Y an n by p matrix with the observations
#' @param ord a total ordering of 1, ..., p which will be tested
#' @param K a vector containing integer values greater than 2 which will be used for the test statistic
#' @param bs the number of bootstrap resamples for the null distribution
#' @param aggType the aggregation used for the test statistic
#' #' \itemize{
#' \item "inf": L-inf 
#' \item "1" L-1
#' \item "2" L-2 
#' \item "3" Computes all of the above
#' }
#' @param intercept T indicates that an intercept should be used; 0 indicates no intercept is used
#' @return
#' a vector of p-values corresponding to test 
testOrdering <- function(Y, ord, K = c(2, 3), bs = 500, aggType = 1, intercept = T){
  n <- dim(Y)[1]
  p <- dim(Y)[2]
    
  ## Check if aggType is correct ##
  if(!(aggType %in% c("inf", 1, 2, 3))){
    stop("Incorrect aggType: must be either 'inf', '1', '2', '3'")
  }
  
  ## Check if intercept is correct ##
  if(!(intercept %in% c(T, F))){
    stop("intercept must be either TRUE or FALSE")
  }
  

  ## Check if K is correct ##
  if(any(K < 2) | any(K != round(K))  ){
    stop("K must contain integers greater than or equal to 2")
  }
  
  
  ### Change aggType and intercept to numeric for passing onto helpers   
  if(aggType == "inf"){
    aggType <- 0
  }
  # make intercept 0/1
  intercept <- intercept + 0
  
  
  ## Initialize final p-values to NA
  fisherInfPval <-  tipettInfPval <-
    fisherOnePval <- tipettOnePval <-
    fisherTwoPval <- tipettTwoPval <- NA
  
  
  fisherInfStat <-  tipettInfStat <-
    fisherOneStat <- tipettOneStat <-
    fisherTwoStat <- tipettTwoStat <- NA
  
  
  ## all sets of "ancestors"
  ## minus 1 for c++ indexing
  
  if(length(unique(ord)) < length(ord)){
    print(ord)
  }
  
  ordSets <- sapply(2:p, function(x){ord[1:x] - 1})
  


  
  out <- lapply(ordSets, cdcs::exhaustiveHelperMulti, dat = Y, K = K, aggType = aggType,
                bs = bs, intercept = intercept)
  
  
  if(aggType %in% c(0, 3)){
    pValTotalInf <- sapply(out, function(res){mean(res$nullDistInf >= res$testStatInf)})
    pValTotalInf <- cdcs::convToCont(pValTotalInf, bs)
    
    fisherInfStat <- -2 * sum(log(pValTotalInf))
    tipettInfStat <- min(pValTotalInf)
    
    fisherInfPval <- pchisq(fisherInfStat, df = 2 * (p - 1) )
    tipettInfPval <- pbeta(tipettInfStat, 1, p-1)
    
  } 
  
  if (aggType %in% c(1, 3)) {
  
    pValTotalOne <- sapply(out, function(res){mean(res$nullDistOne >= res$testStatOne)})
    
    pValTotalOne <- cdcs::convToCont(pValTotalOne, bs)
    
    
    fisherOneStat <- -2 * sum(log(pValTotalOne))
    tipettOneStat <- min(pValTotalOne)
    fisherOnePval <- pchisq(fisherOneStat, df =2 *(p - 1))
    tipettOnePval <- pbeta(tipettOneStat, 1, p-1)
    
  } 
  
  if (aggType %in% c(2, 3)) {
    

    pValTotalTwo <- sapply(out, function(res){mean(res$nullDistTwo >= res$testStatTwo)})
    
    pValTotalTwo <- cdcs::convToCont(pValTotalTwo, bs)
    
    fisherTwoStat <- -2 * sum(log(pValTotalTwo))
    tipettTwoStat <- min(pValTotalTwo)
    fisherTwoPval <- pchisq(fisherTwoStat, df = 2 * (p - 1))
    tipettTwoPval <- pbeta(tipettTwoStat, 1, p-1)
    
  }
  
  
ret <- c(fisherInfStat = fisherInfStat,
         fisherOneStat = fisherOneStat,
         fisherTwoStat = fisherTwoStat,
         tipettInfStat = tipettInfStat,
         tipettOneStat = tipettOneStat,
         tipettTwoStat = tipettTwoStat,
         
         fisherInfPval = fisherInfPval,
         fisherOnePval = fisherOnePval,
         fisherTwoPval = fisherTwoPval,
         tipettInfPval = tipettInfPval,
         tipettOnePval = tipettOnePval,
         tipettTwoPval = tipettTwoPval)


return(ret)

  
}







#' Confidence sets for Causal Discovery
#'
#'
#' Tests all total orderings to see if they are consistent with the data
#'
#' @param Y an n by p matrix with the observations
#' @param K a vector containing integer values greater than 2 which will be used for the test statistic
#' @param bs the number of bootstrap resamples for the null distribution
#' @param aggType the aggregation used for the test statistic
#' #' \itemize{
#' \item "inf": \eqn{\max_{j \prec v} \vert (Y_j^K)^T \eta_v \vert} 
#' \item "1" \eqn{\sum_{j \prec v} \vert (Y_j^K)^T \eta_v \vert}
#' \item "2" \eqn{\sum_{j \prec v} \vert (Y_j^K)^T \eta_v \vert}
#' \item "3" \eqn{\sum_{j \prec v} \vert (Y_j^K)^T \eta_v \vert}
#' }
#' @param intercept 1 indicates that an intercept should be used; 0 indicates no intercept is used
#' @return
#' a matrix containing p-values corresponding to various tests for all permutations of 1,..., p 
exhaustiveTest <- function(Y, K = c(2,3), aggType = "all", bs = 200, intercept = 1){
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  
  
  ## Check if aggType is correct ##
  if(!(aggType %in% c("inf", 1, 2, 3))){
    stop("Incorrect aggType: must be either 'inf', '1', '2', '3'")
  }
  
  ## Check if intercept is correct ##
  if(!(intercept %in% c(T, F))){
    stop("intercept must be either TRUE or FALSE")
  }
  
  
  ## Check if K is correct ##
  if(any(K < 2) | any(K != round(K))  ){
    stop("K must contain integers greater than or equal to 2")
  }
  
  
  ### Change aggType and intercept to numeric for passing onto helpers   
  if(aggType == "inf"){
    aggType <- 0
  }
  # make intercept 0/1
  intercept <- intercept + 0
  
  
  
  # Same as testOrdering, but doesn't require passing in other data
  .exhaustiveInner <- function(ordering){
    
    ## Initialize final p-values to NA
    fisherInfPval <-  tipettInfPval <-
      fisherOnePval <- tipettOnePval <-
      fisherTwoPval <- tipettTwoPval <- NA
    
    
    fisherInfStat <-  tipettInfStat <-
      fisherOneStat <- tipettOneStat <-
      fisherTwoStat <- tipettTwoStat <- NA
    
    
    ## all sets of "ancestors"
    ## minus 1 for c++ indexing
    ordSets <- sapply(2:p, function(x){ord[1:x] - 1})
    
    
    
    
    out <- lapply(ordSets, cdcs::exhaustiveHelperMulti, dat = Y, K = K, aggType = aggType, bs = bs, intercept = intercept)
    
    
    if(aggType %in% c(0, 3)){
      pValTotalInf <- sapply(out, function(res){mean(res$nullDistInf >= res$testStatInf)})
      pValTotalInf <- cdcs::convToCont(pValTotalInf, bs)
      
      fisherInfStat <- -2 * sum(log(pValTotalInf))
      tipettInfStat <- min(pValTotalInf)
      
      fisherInfPval <- pchisq(fisherInfStat, df = 2 * (p - 1) )
      tipettInfPval <- pbeta(tipettInfStat, 1, p-1)
      
    } 
    
    if (aggType %in% c(1, 3)) {
      
      pValTotalOne <- sapply(out, function(res){mean(res$nullDistOne >= res$testStatOne)})
      
      pValTotalOne <- cdcs::convToCont(pValTotalOne, bs)
      
      
      fisherOneStat <- -2 * sum(log(pValTotalOne))
      tipettOneStat <- min(pValTotalOne)
      fisherOnePval <- pchisq(fisherOneStat, df =2 *(p - 1))
      tipettOnePval <- pbeta(tipettOneStat, 1, p-1)
      
    } 
    
    if (aggType %in% c(2, 3)) {
      
      
      pValTotalTwo <- sapply(out, function(res){mean(res$nullDistTwo >= res$testStatTwo)})
      
      pValTotalTwo <- cdcs::convToCont(pValTotalTwo, bs)
      
      fisherTwoStat <- -2 * sum(log(pValTotalTwo))
      tipettTwoStat <- min(pValTotalTwo)
      fisherTwoPval <- pchisq(fisherTwoStat, df = 2 * (p - 1))
      tipettTwoPval <- pbeta(tipettTwoStat, 1, p-1)
      
    }
    
    
    ret <- c(fisherInfStat = fisherInfStat,
             fisherOneStat = fisherOneStat,
             fisherTwoStat = fisherTwoStat,
             tipettInfStat = tipettInfStat,
             tipettOneStat = tipettOneStat,
             tipettTwoStat = tipettTwoStat,
             
             fisherInfPval = fisherInfPval,
             fisherOnePval = fisherOnePval,
             fisherTwoPval = fisherTwoPval,
             tipettInfPval = tipettInfPval,
             tipettOnePval = tipettInfPval,
             tipettTwoPval = tipettInfPval)
    
    
    return(ret)
  }
  

  
  allPerms <- combinat::permn(1:p)
  allPermPvals <- t(sapply(allPerms, FUN = .exhaustiveInner))
  ret <- cbind(allPermPvals, do.call("rbind",allPerms))
  colnames(ret) <- c("fisherInfStat", "fisherOneStat", "fisherTwoStat",
                     "tipettInfStat", "tipettOneStat", "tipettTwoStat",
                     "fisherInfPval", "fisherOnePval", "fisherTwoPval",
                     "tipettInfPval", "tipettOnePval", "tipettTwoPval",
                     paste("V", 1:p, sep = "")
                     )
  
  return(ret)
}