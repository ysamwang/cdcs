## convert discrete pvalues to continuous uniform
# by sampling uniformly within the discrete buckets
convToCont <- function(pval, bs){

  cutPoints <- 1 / (bs + 1) * 0:(bs+1)
  
  ind <- pval * bs + 1

  return(runif(length(ind), cutPoints[ind], cutPoints[ind + 1]))
}
