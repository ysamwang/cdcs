### Updated ###
# 4/25/24
# Max stat
# write to nl_exp_hard_max_*.csv
#

runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)



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

run.onceBnb <- function(p.list = c(7), n, distro, bs = 400, parent_prob = 1/4, j.list = c(5, 25, 50),
                        funcType = "GP", h = 1, Mtrue= 7, k.list =c(1, 2, 3, 4, 5),
                        basis = "bspline") {
  
  
  rec <- matrix(0, 0, 8)
  colnames(rec) <- c("p","K", "J", "testFunc", "size", "cover", "ancest", "time")
  
  for(pInd in 1:length(p.list)){
    p <- p.list[pInd]
    
    out <- cdcs::rDAG_anm(p = p, n = n, parent_prob = parent_prob,
                          funcType = funcType, M = Mtrue, h = h, dist = distro,
                          lowScale = 1/5, highScale = sqrt(2)/5, noParentMult = 5)
    Y <- scale(out$Y)
    
    
    G1 <- array(0, dim = c(n, max(j.list) * 2, p) )
    G2 <- array(0, dim = c(n, 3, p) )
    G3 <- array(0, dim = c(n, 7, p) )
    
    
    
    # Create test functions
    for(u in 1:p){
      
      for(j in 1:max(j.list)){
        
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
    
    
    ## go through all j.list length
    for(jInd in 1:length(j.list)){
      
      for(k in 1:length(k.list)){
        time.rec <- system.time(out <- cdcs::brandAndBound_anm(Y, G1[ , 1:(2*j.list[jInd]), , drop = F], 
                                                               bs =bs, aggType = 3, alpha = .1,
                                                               pValueAgg = "tippet", intercept = 1,
                                                               verbose = F, basis = basis,
                                                               K = k.list[k]))[3]
        
        rec <- rbind(rec,
                     c(p, k.list[k], j.list[jInd], "trig", sum(out$pValue > .1), all(out[1, -1] == 1:p),
                       mean(cdcs::getAncest(out[which(out$pValue > .1),-1])[lower.tri(matrix(0, p,p))] == 1 ),
                       time.rec))
      }
    }
    
    ## moments
    for(k in 1:length(k.list)){
      
      time.rec <- system.time(out <- cdcs::brandAndBound_anm(Y, G2, 
                                                             bs =bs, aggType = 3, alpha = .1,
                                                             pValueAgg = "tippet", intercept = 1,
                                                             verbose = F, basis = basis,
                                                             K = k.list[k]))[3]
      
      rec <- rbind(rec,
                   c(p, k.list[k], 3, "moment", sum(out$pValue > .1), all(out[1, -1] == 1:p),
                     mean(cdcs::getAncest(out[which(out$pValue > .1),-1])[lower.tri(matrix(0, p,p))] == 1 ),
                     time.rec))
    }
    
    # Zoo
    for(k in 1:length(k.list)){
      
      time.rec <- system.time(out <- cdcs::brandAndBound_anm(Y, G3, 
                                                             bs =bs, aggType = 3, alpha = .1,
                                                             pValueAgg = "tippet", intercept = 1,
                                                             verbose = F, basis = basis,
                                                             K = k.list[k]))[3]
      
      rec <- rbind(rec,
                   c(p, k.list[k], 8, "zoo", sum(out$pValue > .1), all(out[1, -1] == 1:p),
                     mean(cdcs::getAncest(out[which(out$pValue > .1),-1])[lower.tri(matrix(0, p,p))] == 1 ),
                     time.rec))
    }
    
    # Comb
    for(k in 1:length(k.list)){
      
      time.rec <- system.time(out <- cdcs::brandAndBound_anm(Y, G4, 
                                                             bs =bs, aggType = 3, alpha = .1,
                                                             pValueAgg = "tippet", intercept = 1,
                                                             verbose = F, basis = basis,
                                                             K = k.list[k]))[3]
      
      rec <- rbind(rec,
                   c(p, k.list[k], j.list[jInd] + 3 + 8, "comb", sum(out$pValue > .1), all(out[1, -1] == 1:p),
                     mean(cdcs::getAncest(out[which(out$pValue > .1),-1])[lower.tri(matrix(0, p,p))] == 1 ),
                     time.rec))
    }
    
    
    
  }
  
  return(rec)
  
}


# system.time(out <- run.onceBnb(p.list = c(5), n =2500, distro = "gamma",
#             h = 1, parent_prob = 3/4,
#             funcType = "cam", basis = "bspline",
#             m.list = c(30, 40), k.list = c(5, 10)))


###########################
library(cdcs)

#runInd <- 800 + runInd

sample.size <- 400
rep.runs <- 2
n.list <- c(2500, 5000, 7500, 10000)
d.list <- c("gamma", "laplace")
func.list <- c("cam")
basis.list <- c("bspline")
prob.list <- c(1/3)
p.list <- c(7)

param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), d.list, func.list, basis.list,
                          prob.list)



n <- param.grid[runInd, 1]
distro <- param.grid[runInd, 2]
funcType <- param.grid[runInd, 3]
basis <- param.grid[runInd, 4]
parentProb <- param.grid[runInd, 5]


if(basis == "bspline"){
  k.list <- c(20, 40, 60)
  j.list <- c(25)
} else {
  k.list <- c(2, 3, 4, 5)
  j.list <- c(10)
}


out <- lapply(1:rep.runs, function(x){run.onceBnb(p.list = p.list, n =n, distro = distro,
                                                  parent_prob = parentProb, funcType = funcType,
                                                  basis = basis, Mtrue = 5,
                                                  j.list = j.list, k.list = k.list)})



tab <- do.call("rbind", out)
outTab <- data.frame(n, distro, basis, funcType, parentProb, tab)
write.csv(outTab, paste("../results/nl/nl_exp_hard_",runInd, ".csv", sep = ""))

#################
# 
# sample.size <- 400
# rep.runs <- 2
# n.list <- c(2500, 5000, 7500, 10000)
# d.list <- c("gamma", "laplace")
# func.list <- c("cam")
# basis.list <- c("bspline")
# prob.list <- c(1/3)
# p.list <- c(7)
# 
# param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), d.list, func.list, basis.list,
#                           prob.list)
# 
# 
# 
# missing <- c()
# for(runInd in 1:nrow(param.grid)){
#   if(file.exists(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/nl_new/nl_exp_hard_",runInd, ".csv", sep = ""))){
#     temp <- data.frame(read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/nl_new/nl_exp_hard_",runInd, ".csv", sep = ""))[,-1])
#     outTab <- rbind(outTab, temp)
#   } else {
#     missing <- c(missing, runInd)
#   }
# 
# }
# 
# missing
# 
# outTab$cover <- ifelse(is.na(outTab$cover), FALSE, outTab$cover) 
# resTab <- aggregate(cbind(size, cover) ~ n + m + parentProb + k + distro +p, dat = outTab, FUN = mean)
# resTab$size <- resTab$size / factorial(resTab$p)
# resTab$size <- round(resTab$size, 2)
# out <- reshape(resTab[, -which(names(resTab) %in% c("p"))], idvar = c("m", "distro", "parentProb", "n"),
#                direction = "wide",
#                timevar = "k")
# 
# multip <- -qnorm(.025 / 12)
# 
# ratioTab <- cbind(out[, c(4,3,2,1)],
#                   momentCover = ifelse(out$cover.3 > .9 - multip * sqrt(.9 * .1 / 400), round(out$cover.3,3), paste("\\textbf{", round(out$cover.3,3),"}", sep = "")),
#                   Trig10Cover = ifelse(out$cover.10 > .9 - multip * sqrt(.9 * .1 / 400), round(out$cover.10,3), paste("\\textbf{", round(out$cover.10,3),"}", sep = "")),
#                   Trig100Cover = ifelse(out$cover.100 > .9 - multip * sqrt(.9 * .1 / 400), round(out$cover.100,3), paste("\\textbf{", round(out$cover.100,3),"}", sep = "")),
#                   CombCover = ifelse(out$cover.103 > .9 - multip * sqrt(.9 * .1 / 400), round(out$cover.103,3), paste("\\textbf{", round(out$cover.103,3),"}", sep = "")),
# 
#                   momentSize = ifelse(out$cover.3 > .9 - multip * sqrt(.9 * .1 / 400), round(out$size.3,3), ""),
#                   Trig10Size = ifelse(out$cover.10 > .9 - multip * sqrt(.9 * .1 / 400), round(out$size.10,3), ""),
#                   Trig100Size = ifelse(out$cover.100 > .9 - multip * sqrt(.9 * .1 / 400), round(out$size.100,3), ""),
#                   CombSize = ifelse(out$cover.103 > .9 - multip * sqrt(.9 * .1 / 400), round(out$size.103,3), ""))
# 
# library(tidyverse)
# library(kableExtra)
# 
# printTab <- ratioTab
# 
# printTab$parentProb <- ifelse(printTab$parentProb == 1/4, "$\\pi = 1/4$", "$\\pi = 3/4$")
# kbl(printTab[which(printTab$distro == "gamma"), -1], format = "latex", booktabs = T, align = "c",
#     linesep = c(rep("", 3), "\\addlinespace"), digits = 3, escape = F) %>%
#   add_header_above(c(" " = 3, "Coverage" = 4, "Prop of Unrejected" = 4)) %>%
#   collapse_rows(columns = c(1,2), valign = "middle")
# 
# kbl(printTab[which(printTab$distro == "laplace"), -1], format = "latex", booktabs = T, align = "c",
#     linesep = c(rep("", 3), "\\addlinespace"), digits = 3, escape = F, row.names = F) %>%
#   add_header_above(c(" " = 3, "Coverage" = 4, "Prop of Unrejected" = 4)) %>%
#   collapse_rows(columns = c(1,2), valign = "middle")
# 

