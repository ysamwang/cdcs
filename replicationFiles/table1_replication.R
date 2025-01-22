### Replication file for Table 1 in main manuscript ###
#
# Performs goodness of fit test for single regression
# Tests empirical level when the null hypothesis is true
# and tests power when null is false
#


source("comparison/HOLS_procedure.R") # source code for HOLS
library(hdi)
library(cdcs)
library(RPtests)
library(IndepTest)
library(microbenchmark)
library(tidyverse)


### Function used to run one replication of the simulation ###
# p: the number of variables
# n: the number of samples
# distro: the distribution of the errors (see documentation of rDAG for options)
# parent_prob: the probability of adding an edge in the graph
# verbose: print details at each step 
run.once <- function(p, n, distro, bs = 500, parent_prob = 1/2, verbose = F){
  
  # matrix to record results
  rec <- matrix(0, nrow = 7, ncol = 3)
  colnames(rec) <- c("null", "alt", "time")
  
  
  # Generate data from linear SEM
  # see ?rDAG for parameters
  # Graph has a unique causal ordering of 1:p
  dat <- cdcs::rDAG(p, n, parent_prob = parent_prob,
                    lowScale = .8, highScale = 1,
                    lowEdge = .1, highEdge = .95,
                    dist = distro, uniqueTop = T)
  
  # scale and center data
  Y <- scale(dat$Y)
  
  ### Test when the Null Hypothesis is true ###
  # H0: pa(p) \subseteq 1:(p-1) \subseteq nd(p)
  # i.e., testing whether 1:(p-1) are the parents of p 
  ind <- p
  child <- Y[, ind, drop = F]
  parents <- Y[, -ind, drop = F]
  
  
  # Generate test functions for proposed procedure
  # used for bnbHelperanm
  J1 <- 2
  G1 <- array(0, dim = c(n, J1 * 2, ncol(parents)) )
  G2 <- array(0, dim = c(n, 3, ncol(parents)) )

  for(u in 1:ncol(parents)){
    
    for(j in 1:J1){
      
      G1[, 2 * j - 1, u] <- sin(parents[, u]*j)
      G1[, 2 * j, u] <- cos(parents[, u]*j)
      
    }
    
    G2[, 1, u] <- scale(parents[, u]^2)
    G2[, 2 , u] <- scale(parents[, u]^3)
    G2[, 3 , u] <- scale(sign(parents[,u])*abs(parents[, u])^(2.5))
  
  }
  
  G3 <- abind::abind(G1, G2, along = 2)

  
  
  ### Record p-value for all tests ###
  rec[1, 1] <- cdcs::bnbHelperanm(parents, child, G = G3, withinAgg = 3, aggType = 3, bs = bs, intercept = 1)$pVals[1]
  rec[2, 1] <- cdcs::senSen2014(parents, child, bs = bs)
  rec[3, 1] <- RPtests::RPtest(parents, child, B = bs, resid_type = "OLS")
  rec[4, 1] <- RPtests::RPtest(parents, child, B = bs, resid_type = "Lasso")
  rec[5, 1] <- IndepTest::MINTregression(parents, child, max(round(n / 20), 3), max(round(n / 10), 5), w=FALSE, rnorm(n * bs))
  rec[6, 1] <- HOLS.check(parents, child, simulated.pval = F, nsim = bs)$pval.glob
  rec[7, 1] <- HOLS.check(parents, child, simulated.pval = T, nsim = bs)$pval.glob.sim
  
  
  ### Test when the Alternative Hypothesis is true ###
  # H0: pa(1) \subseteq 2:p \subseteq nd(1)
  # i.e., testing whether 2:p are the parents of 1 
  ind <- 1
  child <- Y[, ind, drop = F] # child variable
  parents <- Y[, -ind, drop = F] # testing parents
  
  
  # Generate test functions for proposed procedure
  J1 <- 2
  G1 <- array(0, dim = c(n, J1 * 2, ncol(parents)) )
  G2 <- array(0, dim = c(n, 3, ncol(parents)) )
  
  
  for(u in 1:ncol(parents)){
    
    for(j in 1:J1){
      G1[, 2 * j - 1, u] <- sin(parents[, u]*j)
      G1[, 2 * j, u] <- cos(parents[, u]*j)
    }
    
    G2[, 1, u] <- scale(parents[, u]^2)
    G2[, 2 , u] <- scale(parents[, u]^3)
    G2[, 3 , u] <- scale(sign(parents[,u])*abs(parents[, u])^(2.5))
    
  }
  
  G3 <- abind::abind(G1, G2, along = 2)
  
  
  ### Record p-value for all tests ###
  rec[1, 3] <- microbenchmark::microbenchmark(rec[1, 2] <- cdcs::bnbHelperanm(parents, child, G = G3, withinAgg = 3, aggType = 3, bs = bs, intercept = 1)$pVals[1], times = 1)$time
  rec[2, 3] <- microbenchmark::microbenchmark(rec[2, 2] <- cdcs::senSen2014(parents, child, bs = bs), times = 1)$time
  rec[3, 3] <- microbenchmark::microbenchmark(rec[3, 2] <- RPtests::RPtest(parents, child, B = bs, resid_type = "OLS"), times = 1)$time
  rec[4, 3] <- microbenchmark::microbenchmark(rec[4, 2] <- RPtests::RPtest(parents, child, B = bs, resid_type = "Lasso"), times = 1)$time
  rec[5, 3] <- microbenchmark::microbenchmark(rec[5, 2] <- IndepTest::MINTregression(parents, child, max(round(n / 20), 3), max(round(n / 10), 5), w=FALSE, rnorm(n * bs)), times = 1)$time
  rec[6, 3] <- microbenchmark::microbenchmark(rec[6, 2] <- HOLS.check(parents, child, simulated.pval = F, nsim = bs)$pval.glob, times = 1)$time
  rec[7, 3] <- microbenchmark::microbenchmark(rec[7, 2] <- HOLS.check(parents, child, simulated.pval = T, nsim = bs)$pval.glob.sim, times = 1)$time
  
  
  rec <- data.frame(p = p, n = n, distro = toString(distro), parent_prob = parent_prob,
                    method = c("W", "SS", "RPols", "RPlasso", "Mint", "HOLS", "HOLSsim"), rec)
  
  return(rec)
  
}

  

##################

### Settings from the paper ###
# sample.size <- 500
# p.list <- c(10, 15, 20, 30, 45)
# d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace", "mixed")
# pow.list <- c(2, 5/4)
# pp.list <- c(1/2)


### Faster settings ###
sample.size <- 5 # Number of replicates
# vector which contains settings for p, the number of variables
p.list <- c(5) 
# vector which contains the distribution of the errors
d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace", "mixed")
# number of samples based on p; i.e., n = round(p^pow.list)
pow.list <- c(2)
# probability of adding an edge between ancestor and descendant 
pp.list <- c(1/2)

param.grid <- expand.grid(rep(p.list, sample.size), d.list, pow.list, pp.list)
colnames(param.grid) <- c("p", "distr", "deg", "parent_prob")

## output will be stored in this data frame
out <- data.frame(matrix(0, nrow(param.grid) * 7, 8))
colnames(out) <- c("p", "n", "distro", "parent_prob", "method", "null_pval", "alt_pval", "time")
for(runInd in 1:nrow(param.grid)){
  p <- param.grid[runInd, 1]
  distro <- param.grid[runInd, 2]
  pow.i <- param.grid[runInd, 3]
  parent_prob <- param.grid[runInd, 4]
  n <- round(p^pow.i)
  
  cat("=== n = ")
  cat(n)
  cat("; dist = ")
  cat(toString(distro))
  cat(" ===\n")
  
  out[(1+(runInd-1)*7):(runInd *7), ] <- run.once(p, n, distro, verbose = T)
  
}


### Summary of all runs ###
# method: W is Wang, Kolar Drton; SS is Sen and Sen 2014; RPols and RPlasso
#   are the low-dimensional and high-dimensional variants from 
#   Shah and Buhlmann 2017; HOLS and HOLSsim are from 
#   Schultheiss, Buhlmann, and Yuan 2023
# size: the empirical size of an alpha = .1 hypothesis test
# power: the empirical power of an alpha = .1 hypothesis test
# time: time for computation in seconds
results <- out %>% group_by(p, n, parent_prob, distro, method) %>%
  summarize(size = mean(null_pval < .1),
            power = mean(alt_pval < .1),
            time = mean(time /1e9))

results




