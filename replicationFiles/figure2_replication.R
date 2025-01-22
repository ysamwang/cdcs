### Simulations shown in Figure 1 examining performance of the confidence set procedures


## run.onceBnb is a helper function which takes in
# p: the number of variables
# n: the sample size
# distro: the distribution of the errors
# parent_prob: the probability of an edge beween an ancestor and a descendant
# verbose: whether to print details at every iteration

run.onceBnb <- function(p, n, distro, bs = 400, parent_prob = 1/3, verbose = T){
  
  
  # Sample the data from a linear SEM
  dat <- cdcs::rDAG(p, n, parent_prob = parent_prob, lowScale = .8,
                    highScale = 1, edgeVar = n^(-1/10),
                    dist = distro, uniqueTop = T)
  
  # scale and center the data
  Y <- scale(dat$Y)
  
  # Use direct lingam to calculate a point estimate
  outlingamDirect <- causalXtreme::direct_lingam_search(Y)

  # Form the test functions
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
    
    rec <- matrix(0, 1, 5)
    colnames(rec) <- c("size", "cover", "ancest", "time","pointEst")
    
    # compute the confidence interval
    time.rec <- system.time(out <- cdcs::brandAndBound(Y, G3, bs = bs, aggType = 3, alpha = .1,
                                                       pValueAgg = "tippet", verbose = verbose))[3]
    
    # record the:
    # size of the confidence set
    # whether it covers the true order (1, 2, ..., p)
    # proportion of ancestral relations in the lower confidence set
    # the time (in seconds)
    # whether the point estimate is correct
    rec[1, ] <- c(sum(out$pValue > .1), all(out[1, -1] == 1:p),
                  mean(cdcs::getAncest(out[which(out$pValue > .1),-1])[lower.tri(dat$B)] == 1 ),
                  time.rec, all(outlingamDirect == 1:p))

  return(rec)
}

##################

library(cdcs)

## Settings used to create Figure 1
# p <- 10
# n.list <- c(500, 1000, 2500, 5000)
# sample.size <- 400

## Alternative settings which are faster
p <- 8
n.list <- c(1000)
sample.size <- 2

d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace", "mixed")
param.grid <- expand.grid(rep(n.list, sample.size), d.list)
colnames(param.grid) <- c("n", "distr")

out <- matrix(0, nrow(param.grid), 5)
colnames(out) <- c("size", "cover", "ancest", "time","pointEst")
for(runInd in 1:nrow(param.grid)){
  n <- param.grid[runInd, 1]
  distro <- param.grid[runInd, 2]
  cat("=== n = ")
  cat(n)
  cat("; dist = ")
  cat(toString(distro))
  cat(" ===\n")
  
  out[runInd, ] <- run.onceBnb(p, n, distro, verbose = T)

}



### Summary of all runs ###
# n: the sample size
# distr: the distribution of the errors
# size: the number of orderings in the confidence set
# cover: whether the confidence set contains the true ordering
# ances: the proportion of all ancestal relations which are contained in the lower envelope
# time: time in seconds required to compute the confidence set
# pointEst: 1 if the point estimate from DirectLiNGAM is correct; 0 otherwise
results <- data.frame(param.grid, out) %>% group_by( n, distr) %>%
  summarize(size = mean(size),
            cover = mean(cover),
            ancest = mean(ancest),
            time = mean(time))








