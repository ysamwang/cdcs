#' Confidence sets for Causal Discovery
#'
#'
#' Test whether a specific ordering is consistent with the data
#'
#' @param Y an n by p matrix with the observations
#' @param K a vector containing integer values greater than 2 which will be used for the test statistic
#' @param numProposals The number of proposals to test for the confidence set
#' @param proposalBs the size of the bootstrap sample used to generate the proposals 
#' @param bs the number of bootstrap resamples for the null distribution when testing any propposal
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
bsProp <- function(Y, K = c(2, 3), propBsCutoff = 1000, numProposals = 100, recaptureCutoff = .1, proposalBs,
                   bs = 200, aggType = 1, intercept = T){
  
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
  
  
  
  #### Generate proposals ###
  .bootstrapProp <- function(){
    cdcs::directLiNGAM(Y[sample(n, size = proposalBs, replace = T), , drop  = F], verbose = F)  
  }
  
  
  recapture <- 1
  initProposal <- t(replicate(numProposals, .bootstrapProp()))
  proposals <- mgcv::uniquecombs(initProposal)
  
  ## as long as a significant portion of the newly generated proposals
  ## have not yet been seen, or the total number of proposals is not too large
  while(recapture > recaptureCutoff & dim(proposals)[1] < propBsCutoff){
   
    l <- dim(proposals)[1]
    initProposal <- t(replicate(numProposals, .bootstrapProp()))
    proposals <- mgcv::uniquecombs(rbind(proposals, initProposal))
    recapture <- (dim(proposals)[1] - l) / numProposals

  }
  
  
  ### calculate p-values for all proposals ###
  pvals <- apply(proposals, MAR = 1,
               function(ord){cdcs::testOrdering(dat$Y, ord, K = K,
                                                bs = bs, aggType = aggType)[7:12] } )
  

  
  ### combine ###
  res <- cbind(t(pvals), proposals)
  colnames(res)[7:ncol(res)] <- paste("V", 1:p, sep = "")

  return(res)
  
}
