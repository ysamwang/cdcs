### Replication file for Figure 3 in main manuscript ###
# Simulations examining performance of the confidence set procedures for 
# non-linear SEMs


## run.onceBnb is a helper function which takes in
# p: the number of variables
# n: the sample size
# distro: the distribution of the errors
# parent_prob: the probability of an edge beween an ancestor and a descendant
# verbose: whether to print details at every iteration
# J: number of random fourier features to use
# funcType: type of functions used in SEM. Choices are
#     GP: sample from a gaussian process with parameters sigma and h
#     poly: sample from a random polynomial where the scores for the hermite polynomials have sd 1/k! 
#     fourier: sample from random fourier function where scores are 1/M^4
#     cam: function from Buhlmann et al 2013
#     tanh: tanh function
# h: bandwidth if using a Gaussian process
# Mtrue: degree of polynomial if using poly
# k.list: list of parameters either the df for bsplines or the number of polynomial terms 
# basis: the basis used to model f_v. Choices are "bspline" or "poly"
run.onceBnb <- function(p = 7, n, distro, bs = 400, parent_prob = 1/4, J = c(10),
                        funcType = "GP", h = 1, Mtrue= 7, k.list =c(1, 2, 3, 4, 5),
                        basis = "bspline", verbose = F) {
  
  
  # helper functions for various test functions
  testFunc <- function(x, type = 1){
    if(type == 1){
      return(scale(sin(x^2)))
      
    } else if (type == 2){
      return(scale(cos(x^2)))
      
    } else if (type == 3){
      return(scale(sin(x) * x))
      
    } else if (type == 4){
      return(scale(cos(x) * x))
      
    } else if (type == 5){
      return(scale(sin(x^2) * x))
      
      
    } else if (type == 6){
      return(scale(cos(x^2) * x))
      
    } else if (type == 7){
      return(scale(tanh(x)))
    } 
    
  }
  

  rec <- matrix(0, 0, 5)
  colnames(rec) <- c("K", "size", "cover", "ancest", "time")

  
  ## Generate data
  out <- cdcs::rDAG_anm(p = p, n = n, parent_prob = parent_prob,
                           funcType = funcType, M = Mtrue, h = h, dist = distro,
                           lowScale = 1/5, highScale = sqrt(2)/5, noParentMult = 5)
  Y <- scale(out$Y)
    

    G1 <- array(0, dim = c(n, J * 2, p) )
    G2 <- array(0, dim = c(n, 3, p) )
    G3 <- array(0, dim = c(n, 7, p) )
    
    
  
    # Create test functions
    for(u in 1:p){
      
        for(j in 1:J){
          
          omega <- rnorm(1)
          G1[, 2 * j - 1, u] <- scale(sin(omega * Y[,u]))
          G1[, 2 * j, u] <- scale(cos(Y[,u] * omega))
          
        }
      
      G2[, 1, u] <- scale(Y[, u]^2)
      G2[, 2 , u] <- scale(Y[, u]^3)
      G2[, 3 , u] <- scale(sign(Y[,u])*abs(Y[, u])^(2.5))
      
      for(z in 1:7){
        G3[, z, u] <- testFunc(Y[, u], type =z)
      }

    }
    
    G4 <- abind::abind(G1, G2, G3, along = 2)
    
    
    
    
    # test for each value of k (the size of the basis for the sieve estimator)
    for(k in 1:length(k.list)){
      
      time.rec <- system.time(out <- cdcs::brandAndBound_anm(Y, G4, 
                                                             bs =bs, aggType = 3, alpha = .1,
                                                             pValueAgg = "tippet", intercept = 1,
                                                             basis = basis,
                                                             K = k.list[k], verbose = verbose))[3]
      
      rec <- rbind(rec,
                   c(k.list[k], sum(out$pValue > .1),
                     all(out[1, -1] == 1:p),
                     mean(cdcs::getAncest(out[which(out$pValue > .1),-1])[lower.tri(matrix(0, p,p))] == 1 ),
                     time.rec))
    }
      
  
  return(rec)

}


###########################
library(cdcs)

### Settings used in the paper for polynomial ###
# sample.size <- 400
# n.list <- c(2500, 5000, 7500, 10000)
# p.list <- c(7)
# d.list <- c("gamma", "laplace")
# func.list <- c("poly"")
# basis.list <- c("poly")
# prob.list <- c(1/3)
# k.list <- c(2, 3, 4, 5)
# J <- 10

### Settings used in the paper for sigmoid type ###
# sample.size <- 400
# n.list <- c(2500, 5000, 7500, 10000)
# d.list <- c("gamma", "laplace")
# p.list <- c(7)
# func.list <- c("cam")
# basis.list <- c("bspline")
# prob.list <- c(1/3)
# p.list <- c(7)
# k.list <- c(20, 40, 60)
# J <- 25

sample.size <- 2 # number of replicates
n.list <- c(5000) # list of sample sizes
p.list <- c(7) # list of graph size
d.list <- c("gamma", "laplace") # list of error distributions
func.list <- c("poly") # list of different functions types for SEM
basis.list <- c("poly") # list of different basis types for sieve estimator
prob.list <- c(1/3) # probability of adding an edge between ancestor and descendant 
k.list <- c(2, 5) # size of basis for sieve estimator
J <- 10 # number of random fourier features to use for test functions

param.grid <- expand.grid(rep(n.list, sample.size), p.list, d.list, func.list, basis.list,
                          prob.list)

colnames(param.grid) <- c("n","p", "distr", "funcType", "basisType", "parent_prob")

out <- matrix(0, length(k.list) * nrow(param.grid), 5)
colnames(out) <- c("K", "size", "cover", "ancest", "time")

for(runInd in 1:nrow(param.grid)){
  
  n <- param.grid[runInd, 1]
  p <- param.grid[runInd, 2]
  distro <- param.grid[runInd, 3]
  funcType <- param.grid[runInd, 4]
  basis <- param.grid[runInd, 5]
  parentProb <- param.grid[runInd, 6]
  
  cat("=== n = ")
  cat(n)
  cat("; dist = ")
  cat(toString(distro))
  cat(" ===\n")
  
  out[(1+(runInd-1)*length(k.list)):(runInd *length(k.list)), ] <- run.onceBnb(p = p, n =n, distro = distro,
                               parent_prob = parentProb, funcType = funcType,
                               basis = basis, Mtrue = 5,
                               J = J, k.list = k.list, verbose = T)
  
}


### Summary of all runs ###
# p: the size of the graph
# n: the sample size
# distr: the distribution of the errors
# funcType: the type of function used in the SEM for generating data
# basisType: the type of basis used for the sieve estimator
# K: the size of the basis used for the sieve estimator
# size: the number of orderings in the confidence set
# cover: whether the confidence set contains the true ordering
# ances: the proportion of all ancestal relations which are contained in the lower envelope
# time: time in seconds required to compute the confidence set
results <- data.frame(param.grid, out) %>% group_by(p, n, distr, funcType, basisType, K) %>%
  summarize(size = mean(size),
            cover = mean(cover),
            ancest = mean(ancest),
            time = mean(time))

results %>% print(n = 10)




