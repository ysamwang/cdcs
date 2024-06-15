#' Confidence sets for Causal Discovery
#'
#'
#' Take a set of orderings and see which ancestral relations hold
#'
#' @param tab a matrix with p columns (corresponding to variables) and each row is an ordering
#' @return
#' a p x p matrix where A[i,j] indicates the proportion of orderings where j precedes i
getAncest <- function(tab){
  
  # if the confidence set is empty, then return 0 
  if(nrow(tab) == 0){
    p <- dim(tab)[2]
    return(matrix(0, p, p))
  }
  
  ## Checks if first column of tab is p-values or not
  # if so, remove first column
  if(!all(tab[,1] == round(tab[,1]))){
    tab <- tab[, -1]
  }
  
    # takes in a list of pairs
    # returns proportion of time the first node in the pair preceedes the 2nd node 
  .estOrder <- function(ind){
    check <- apply(tab, MAR = 1, function(x){which(x == ind[1]) < which(x == ind[2])})
    
    return(mean(check))
  }
  
  p <- dim(tab)[2]
  
  pairs <- combn(1:p, 2)
  checkRes <- apply(pairs, MAR = 2, .estOrder)
  
  A <- matrix(0, p, p)
  A[cbind(pairs[2,], pairs[1,])] <- checkRes
  
  return(A)
  
}