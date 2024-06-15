#' Confidence sets for Causal Discovery
#'
#'
#' Take a set of orderings and get the confidence interval for a specific
#' causal effect which accounts for model uncertainty
#'
#' @param tab a matrix with p columns (corresponding to variables) and each row is an ordering
#' @param treatment index of treatment of interest
#' @param outcome index of outcome of interest
#' @param effectType either 'total' or 'direct'
#' @param alpha the confidence level (conditional on the selected graph) 
#' 
#' @return
#' ci is the interval
#' length is the length
ci_modSelect <- function(tab, treat, outc, effectType = "total", alpha = .05, Y){
  
  ## Checks if first column of tab is p-values or not
  # if so, remove first column
  if(any(tab[,1] != round(tab[,1]))){
    
    tab <- tab[, -1]
    
  }
  
  n <- nrow(Y)
  n_orderings <- nrow(tab)
  
  if(nrow(tab) > 0){
    
    # rows in tab where treatment preceeds the outcome
    ind <- apply(tab, MAR = 1, function(x){which(x == treat) < which(x == outc)})
    
    if(sum(ind) > 0){
      
      if(effectType == "total"){
        
        ## Extra processing to force the object into a list 
        adj_Set <- unique(lapply(
          apply(tab[which(ind), , drop = F], MAR = 1, function(x){list(unname(sort(x[1:which(x == treat)])))})
                                 , "[[", 1))
        
      } else if(effectType == "direct"){
        
        ## Extra processing to force the object into a list
        adj_Set <- unique(lapply(
          apply(tab[which(ind), , drop = F], MAR = 1, function(x){list(unname(sort(x[1:(which(x == outc) - 1)])))})
          ,"[[", 1))
        
        
      } 
      
      
      helper1 <- function(adj){
        
        ii <- which(adj == treat)
        
        mod <- RcppArmadillo::fastLmPure(Y[, adj, drop = F], Y[, outc, drop = F])
        mult <- -qt(alpha/2, df = n - length(adj))
        
        
        if(kappa(Y[, adj, drop = F]) > 1e8){
          print("Err! (Cond Number)")
          print(adj)
          print(outc)
          print(tab)
        }
        
        if(outc %in% adj){
          print("Err! (Outcome in Adj)")
          print(adj)
          print(outc)
          print(tab)
        }
        
        if( length(unique(adj)) != length(adj) ){
          print("Err! (Dupe in Adj)")
          print(adj)
          print(outc)
          print(tab)
        }
        return(matrix(c(mod$coeff[ii] - mult * mod$stderr[ii], mod$coeff[ii] + mult * mod$stderr[ii]), 1, 2))
      }
      
      
      
      ## get CI for each adjustment set
      ci <- t(sapply(adj_Set, helper1))
      
      
      # If there are orderings where outc preceeds treatment, also include 0
      if(sum(ind) < n_orderings){
        ci <- rbind(ci, c(0, 0))
      }
      

      ci_final <- intervals::interval_union(intervals::Intervals(ci))
      
      
      return(list(ci = ci_final, length = sum(intervals::size(ci_final)), numAdjSets = length(adj_Set)))
      
      
    } else {
      
      # no non-rejected orderings where treatment precedes outcome
      return(list(ci = intervals::Intervals(c(0,0)), length = 0, numAdjSets = 1))
      
    }
  } else {
    
    # no not rejected orderings at all
    return(list(ci = intervals::Intervals(c(0,0)), length = 0, numAdjSets = 0))
  }
}
