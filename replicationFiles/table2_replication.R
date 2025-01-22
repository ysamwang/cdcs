### Replication file for Table 1 in main manuscript ###
#
#  Confidence intervals which account for model selection uncertainty

library(cdcs)
library(intervals)
library(tidyverse)



### Function used to run one replication of the simulation ###
# p: the number of variables
# n: the number of samples
# distro: the distribution of the errors (see documentation of rDAG for options)
# bs: the number of bootstrap samples to take
# parent_prob: the probability of adding an edge in the graph
# verbose: print details at each step 
run.onceCI <- function(p, n, distro, bs = 400, parent_prob = 1/3, verbose = F){
  
  # index of treatment
  treatment <- 4
  # index of outcome
  outcome <- 7
  alpha <- .1
  
  
  # Generate data from linear SEM
  # see ?rDAG for parameters
  # Graph has a unique causal ordering of 1:p
  dat <- cdcs::rDAG(p, n, parent_prob = 1/3, lowScale = .8,
                    highScale = 1, edgeVar = 1/2,
                    dist = distro, uniqueTop = T)

  
  # Caluclate the true total effect of 4 onto 7
  trueParam <- solve(diag(rep(1, p)) - dat$B)[outcome, treatment]
  
  # scale data
  Y <- scale(dat$Y)
  
  # Get point estimate of causal ordering
  outlingamDirect <- causalXtreme::direct_lingam_search(Y)

  
  # setup test functions  
  J <- 2
  G1 <- array(0, dim = c(n, J * 2, p) )
  G2 <- array(0, dim = c(n, 3, p) )
  
  for(u in 1:p){
    
    for(j in 1:J){
      
      G1[, 2 * j - 1, u] <- sin(j * Y[,u])
      G1[, 2 * j, u] <- cos(j * Y[,u])
      
    }
    
    G2[, 1, u] <- scale(Y[, u]^2)
    G2[, 2 , u] <- scale(Y[, u]^3)
    G2[, 3 , u] <- scale(sign(Y[,u])*abs(Y[, u])^(2.5))
    
  }
  
  G3 <- abind::abind(G1, G2, along = 2)
  
  # calculate confidence set of causal orderings
  out <- cdcs::brandAndBound(Y, G3, bs = bs, aggType = 3, alpha = .1,
                             pValueAgg = "tippet", verbose = verbose)
  
  # contains 4 columns
  # col 1: does honest CI contain the truth
  # col 2: length of honest CI
  # col 3: does naive CI contain the truth
  # col 4: length of naive CI
  # col 5: is the adjustment set from the naive method valid?
  rec <- rep(0, 5)
  names(rec) <- c("cov_honest", "length_honest", "cov_naive", "length_naive", "pointEst")
  
  ## if the confidence set is non-empty, compute the honest CI
  if(nrow(out) > 0){
    
    # use centered, but not scaled data, since scaling will change the estimated regressino coefficients
    Y.centered <- scale(dat$Y, scale = F)
    # get confidence intervals for each causal ordering in the confidence set
    ci_out <- cdcs::ci_modSelect(out, treatment, outcome, effectType = "total", alpha = alpha, Y = Y.centered)
    
    # Check if true parameter is in the interval
    rec[1:2] <- c(intervals::distance_to_nearest(from = trueParam, to = ci_out$ci) == 0, ci_out$length)
    
  } else{
    
    ## if confidence set is empty, then does not cover and has length 0
    rec[1:2] <- c(0,0)
  }

  
  
  ## if the treatment preceeds the outcome in the point estimate of the
  ## causal ordering, then calculate a naive CI 
  if(which(outlingamDirect == treatment) < which(outlingamDirect == outcome)){
    
    ## adjustment set implied by the point estimate
    adjPoint <- outlingamDirect[1:which(outlingamDirect == treatment)]
    ii <- which(adjPoint == treatment)
    
    # Regress outcome onto the adjustment set and calculate the confidence interval 
    mod <- RcppArmadillo::fastLmPure(Y.centered[, adjPoint, drop = F], Y.centered[, outcome, drop = F])
    mult <- -qt(alpha, df = n - length(adjPoint))
    ci_Point <- c(mod$coeff[ii] - mult * mod$stderr[ii], mod$coeff[ii] + mult * mod$stderr[ii])
    
    # check if the interval contains the truth and calulate the length
    rec[3:4] <- c(ci_Point[1] < trueParam & ci_Point[2] > trueParam, ci_Point[2] - ci_Point[1])
    #
    trueAdjSet <- c(treatment, which(dat$B[treatment, ] != 0) )
    rec[5] <- all(trueAdjSet %in% adjPoint) & all(adjPoint %in% 1:treatment)
    
  } else {

    ## if the treatment does not preceed the outcome in the point estimate
    ## then the CI is empty
    rec[3:5] <- c(0, 0, 0)  
  }
    

  
  return(rec)
}

##################

# sample.size <- 400
# rep.runs <- 5
# p <- 10
# n.list <- c(250, 500, 1000, 2000)
# d.list <- c("laplace", "gamma")
# param.grid <- expand.grid(rep(n.list, sample.size), d.list)

sample.size <- 2
p <- 8
n.list <- c(1000)
d.list <- c("laplace", "gamma")
param.grid <- expand.grid(rep(n.list, sample.size), d.list)


colnames(param.grid) <- c("n", "distr")

## output will be stored in this data frame
out <- data.frame(matrix(0, nrow(param.grid), 5))
colnames(out) <- c("cov_honest",
                   "length_honest", "cov_naive", "length_naive", "adjSet")

for(runInd in 1:nrow(param.grid)){
  
  n <- param.grid[runInd, 1]
  distro <- param.grid[runInd, 2]
  
  cat("=== n = ")
  cat(n)
  cat("; dist = ")
  cat(toString(distro))
  cat(" ===\n")
  
  out[runInd, ] <- run.onceCI(p, n, distro, verbose = T)
  
}


### Summary of all runs ###
# cov_honest: coverage probability of honest CI
# length_honest: average length of honest CI
# cov_naive: coverage probability of naive CI
# length_naive: average length of naive CI
# adjSet: the proportion of times the point estimate yields a valid adjustment set
results <- data.frame(param.grid, out) %>% group_by(n, distr) %>%
  summarize(cov_honest = mean(cov_honest),
            length_honest = mean(length_honest),
            cov_naive = mean(cov_naive),
            length_naive = mean(length_naive),
            adjSet = mean(adjSet))

results






