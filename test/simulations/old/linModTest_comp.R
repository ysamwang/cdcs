### Running comparison for testing single linear regression
# 
#
#
runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)




run.once <- function(p, n, distro, bs = 500, parent_prob = 1/2, verbose = F,
                         cutoff = NULL){
  source("comparison/HOLS_procedure.R")
  k.list <- c(1:6)
  
  gFunc <- function(x, type){
    if(type == 1){
      return(scale(x^2))
    } else if (type ==2 ){
      return(scale(x^3))
    } else if (type ==3 ){
      return(scale(abs(x)^(2.5) * sign(x)))
    }  else if (type == 4){
      return(scale(sin(x) * x))
      
    } else if (type == 5){
      return(scale(cos(x) * x))
      
    } else if (type == 6){
      return(scale(tanh(x)))
    } else if(type == 7){
      return(scale(sin(x^2)))
      
    } else if (type == 8){
      return(scale(cos(x^2)))
      
    } else if (type == 9){
      return(scale(sin(x^2) * x))
      
    } else if (type == 10){
      return(scale(cos(x^2) * x))
      
    } 
    
  }
  
  rec <- matrix(0, nrow = 1, ncol = 30)
  colnames(rec) <- paste(rep(c("null_", "alt_", "timeA_"), each = 10), c("dkwPow1", "dkwPow2","dkwPow1a", "dkwPow2a", "senSen",
                                                                        "rptest_ols","rptest_lasso",  "mint" , "hols", "holsSim"), sep = "")
  
  dat <- cdcs::rDAG(p, n, parent_prob = parent_prob, lowScale = .8,
                    highScale = 1, lowEdge = .1, highEdge = .95,
                    dist = distro, uniqueTop = T)
  
  Y <- scale(dat$Y)
  
  ### Null Hyp
  ind <- p
  child <- Y[, ind, drop = F]
  parents <- Y[, -ind, drop = F]
  
  G <- array(0, dim = c(n , length(k.list), ncol(parents)))

  for(j in 1:ncol(parents)){
    for(k in 1:length(k.list)){
      G[ , k, j] <- gFunc(parents[, j], k.list[k])  
    }
  }
  
         
  rec[1] <- cdcs::bnbHelperanm(parents, child, G = G, withinAgg = 2, aggType = 1, bs = bs, intercept = 1)$pVals[1]
  rec[2] <- cdcs::bnbHelperanm(parents, child, G = G, withinAgg = 2, aggType = 2, bs = bs, intercept = 1)$pVals[1]
  rec[3] <- cdcs::bnbHelperanm_comb(parents, child, G = G, withinAgg = 2, aggType = 1, bs = bs, intercept = 1)$pVals[1]
  rec[4] <- cdcs::bnbHelperanm_comb(parents, child, G = G, withinAgg = 2, aggType = 2, bs = bs, intercept = 1)$pVals[1]
  rec[5] <- cdcs::bnbHelperanm_comb1(parents, child, G = G, withinAgg = 2, aggType = 1, bs = bs, intercept = 1)$pVals[1]
  rec[6] <- cdcs::bnbHelperanm_comb1(parents, child, G = G, withinAgg = 2, aggType = 2, bs = bs, intercept = 1)$pVals[1]
  rec[7] <- RPtests::RPtest(parents, child, B = bs, resid_type = "Lasso")
  rec[8] <- IndepTest::MINTregression(parents, child, max(round(n / 20), 3), max(round(n / 10), 5), w=FALSE, rnorm(n * bs))
  rec[9] <- HOLS.check(parents, child, simulated.pval = F, nsim = bs)$pval.glob
  rec[10] <- HOLS.check(parents, child, simulated.pval = T, nsim = bs)$pval.glob.sim
  
  
  
  ### Alt Hyp
  ind <- 1
  child <- Y[, ind, drop = F]
  parents <- Y[, -ind, drop = F]
  
  G <- array(0, dim = c(n , length(k.list), ncol(parents)))
  
  for(j in 1:ncol(parents)){
    for(k in 1:length(k.list)){
      G[ , k, j] <- gFunc(parents[, j], k.list[k])  
    }
  }
  
  offset <- 10
  
  rec[2 * offset + 1] <- microbenchmark::microbenchmark(rec[offset + 1] <- cdcs::bnbHelperanm(parents, child, G = G, withinAgg = 2,
                                                                           aggType = 1, bs = bs, intercept = 1)$pVals[1], times = 1)$time
  rec[2 * offset + 2] <- microbenchmark::microbenchmark(rec[offset + 2] <- cdcs::bnbHelperanm(parents, child, G = G, withinAgg = 2,
                                                                           aggType = 2, bs = bs, intercept = 1)$pVals[1], times = 1)$time
  rec[2 * offset + 3] <- microbenchmark::microbenchmark(rec[offset + 3] <- cdcs::bnbHelperanm_comb(parents, child, G = G, withinAgg = 2,
                                                                           aggType = 1, bs = bs, intercept = 1)$pVals[1], times = 1)$time
  rec[2 * offset + 4] <- microbenchmark::microbenchmark(rec[offset + 4] <- cdcs::bnbHelperanm_comb(parents, child, G = G, withinAgg = 2,
                                                                           aggType = 2, bs = bs, intercept = 1)$pVals[1], times = 1)$time
  rec[2 * offset + 5] <- microbenchmark::microbenchmark(rec[offset + 5] <- cdcs::bnbHelperanm_comb1(parents, child, G = G, withinAgg = 2,
                                                                                                   aggType = 1, bs = bs, intercept = 1)$pVals[1], times = 1)$time
  rec[2 * offset + 6] <- microbenchmark::microbenchmark(rec[offset + 6] <- cdcs::bnbHelperanm_comb1(parents, child, G = G, withinAgg = 2,
                                                                                                   aggType = 2, bs = bs, intercept = 1)$pVals[1], times = 1)$time
  rec[2 * offset + 5] <- microbenchmark::microbenchmark( rec[offset + 5] <- cdcs::senSen2014(parents, child, bs = bs), times = 1)$time
  rec[2 * offset + 6] <- microbenchmark::microbenchmark(rec[offset + 6] <- RPtests::RPtest(parents, child, B = bs, resid_type = "OLS"), times = 1)$time
  rec[2 * offset + 7] <- microbenchmark::microbenchmark(rec[offset + 7] <- RPtests::RPtest(parents, child, B = bs, resid_type = "Lasso"), times = 1)$time
  rec[2 * offset + 8] <- microbenchmark::microbenchmark( rec[offset + 8] <- IndepTest::MINTregression(parents, child, round(max(n / 20, 3)),
                                                                                                      round(max(n / 10, 5)), w=FALSE, rnorm(n * bs)), times = 1)$time
  rec[2 * offset + 9] <- microbenchmark::microbenchmark( rec[offset + 9] <- HOLS.check(parents, child, simulated.pval = F, nsim = bs)$pval.glob, times = 1)$time
  rec[2 * offset + 10] <- microbenchmark::microbenchmark( rec[offset + 10] <- HOLS.check(parents, child,  simulated.pval = T, nsim = bs)$pval.glob.sim, times = 1)$time


  return(rec)
  
}

  

##################
library(cdcs)
sample.size <- 500
rep.runs <- 50
p.list <- c(10, 15, 20, 30, 45)
d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace", "mixed")
pow.list <- c(2, 5/4)
pp.list <- c(1/2)
param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list, pow.list, pp.list)
nrow(param.grid)
## Param Grid 480


p <- param.grid[runInd, 1]
distro <- param.grid[runInd, 2]
pow.i <- param.grid[runInd, 3]
parent_prob <- param.grid[runInd, 4]
n <- round(p^pow.i)


out <- lapply(1:rep.runs, function(x){run.once(p, n, distro, parent_prob = parent_prob, bs = 500)})
outTab <- data.frame(p, n, distro, parent_prob, pow.i, do.call("rbind", out))

colnames(outTab) <- c("p", "n", "distro", "parent_prob", "pow", 
                      paste(rep(c("null_", "alt_", "timeA_"), each = 10), c("dkwPow1",
                                                                           "dkwPow2","dkwPow1a", "dkwPow2a", "dkwPow1b", "dkwPow2b", 
                                                                           "rptest_lasso",  "mint" , "hols", "holsSim"), sep = ""))

write.csv(outTab, paste("../results/lmgdDAGNew/lmgdDAG_comp_", runInd, ".csv", sep = ""), row.names = F)




# ##########
# sample.size <- 500
# rep.runs <- 100
# p.list <- c(10, 15, 20, 30, 45)
# d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace", "mixed")
# pow.list <- c(2, 5/4)
# pp.list <- c(1/2)
# param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list, pow.list, pp.list)
# nrow(param.grid)
# #
# runInd <- 1
# outTab <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/lmgdDAGNew/lmgdDAG_alt_", runInd, ".csv", sep = ""))
# 
# colnames(outTab) <- c("p", "n", "distro", "parent_prob", "pow",
#                       paste(rep(c("null_", "alt_", "timeA_"), each = 10), c("dkwPow1",
#                                                                             "dkwPow2","dkwPow1a", "dkwPow2a", "senSen",
#                                                                             "rptest_ols","rptest_lasso",  "mint" , "hols", "holsSim"), sep = ""))
# 
# for(runInd in 2:nrow(param.grid)){
#   temp <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/lmgdDAGNew/lmgdDAG_alt_", runInd, ".csv", sep = ""))
#   colnames(temp) <- c("p", "n", "distro", "parent_prob", "pow",
#                         paste(rep(c("null_", "alt_", "timeA_"), each = 10), c("dkwPow1",
#                                                                               "dkwPow2","dkwPow1a", "dkwPow2a", "senSen",
#                                                                               "rptest_ols","rptest_lasso",  "mint" , "hols", "holsSim"), sep = ""))
# 
#   outTab <- rbind(outTab, temp)
# }
# 
# out_Orig <- aggregate(cbind(null_dkwPow1 < .1,
#                             null_dkwPow2 < .1,
#                             alt_dkwPow1 < .1,
#                             alt_dkwPow2 < .1)~ p + distro + parent_prob + pow, data = outTab, FUN = mean)
# 
# sample.size <- 500
# rep.runs <- 50
# p.list <- c(10, 15, 20, 30, 45)
# d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace", "mixed")
# pow.list <- c(2, 5/4)
# pp.list <- c(1/2)
# param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list, pow.list, pp.list)
# nrow(param.grid)
# 
# runInd <- 1
# outTab <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/lmgdDAGNew/lmgdDAG_alt_comb1_", runInd, ".csv", sep = ""))
# 
# colnames(outTab) <- c("p", "n", "distro", "parent_prob", "pow", 
#                       paste(rep(c("null_", "alt_", "timeA_"), each = 10), c("dkwPow1",
#                                                                             "dkwPow2","dkwPow1a", "dkwPow2a", "dkwPow1b", "dkwPow2b", 
#                                                                             "rptest_lasso",  "mint" , "hols", "holsSim"), sep = ""))
# 
# for(runInd in 2:nrow(param.grid)){
#   temp <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/lmgdDAGNew/lmgdDAG_alt_comb1_", runInd, ".csv", sep = ""))
#   colnames(temp) <- c("p", "n", "distro", "parent_prob", "pow", 
#                         paste(rep(c("null_", "alt_", "timeA_"), each = 10), c("dkwPow1",
#                                                                               "dkwPow2","dkwPow1a", "dkwPow2a", "dkwPow1b", "dkwPow2b", 
#                                                                               "rptest_lasso",  "mint" , "hols", "holsSim"), sep = ""))
# 
#   outTab <- rbind(outTab, temp)
# }
# 
# out_comb <- aggregate(cbind(null_dkwPow2 < .1,
#                             null_dkwPow2a < .1,
#                             null_dkwPow2b < .1,
#                             alt_dkwPow2 < .1,
#                             alt_dkwPow2a < .1,
#                             alt_dkwPow2b < .1)~ p + distro + parent_prob + pow, data = outTab, FUN = mean)
# 
# 
# colnames(out_comb) <- c("p", "distro", "parent_prob", "pow", "sizeOrig2",
#                         "sizeComb2", "sizeMax2", "powerOrig2","powerComb2", "powerMax2")
# 
# out_comb
# 
# out_comb$ratio <- out_comb$powerOrig2 / out_comb$powerComb2
# out_comb[which(out_comb$p %in% c(10,20, 45)), ]

# # 
# out_null <- aggregate(cbind(null_dkwPow1 < .1,
#                             null_dkwPow2 < .1,
#                             null_dkwPow1a < .1,
#                             null_dkwPow2a < .1,
#                             null_senSen < .1,
#                             null_rptest_ols < .1,
#                             null_rptest_lasso < .1,
#                             null_mint < .1,
#                             null_hols < .1,
#                             null_holsSim < .1)~ p + distro + parent_prob + pow, data = outTab, FUN = mean)
# names(out_null) <- c("p", "dist", "parent_prob", "pow", "dkwPow1", "dkwPow2", "dkwPow1a", "dkwPow2a", "senSen", "rptest_ols","rptest_lasso", "mint", "hols", "holsSim")
# 
# out_alt <- aggregate(cbind(alt_dkwPow1 < .1,
#                            alt_dkwPow2 < .1,
#                            alt_dkwPow1a < .1,
#                            alt_dkwPow2a < .1,
#                            alt_senSen < .1,
#                            alt_rptest_ols < .1,
#                            alt_rptest_lasso < .1,
#                            alt_mint < .1,
#                            alt_hols < .1,
#                            alt_holsSim < .1)~ p + distro + parent_prob + pow, data = outTab, FUN = mean)
# names(out_alt) <- c("p", "dist", "parent_prob", "pow", "dkwPow1", "dkwPow2", "dkwPow1a", "dkwPow2a", "senSen", "rptest_ols","rptest_lasso", "mint", "hols", "holsSim")
# 
# 
# out_time <- aggregate(cbind(timeA_dkwPow2, timeA_senSen, timeA_rptest_ols,
#                             timeA_rptest_lasso, timeA_mint, timeA_hols) / 1e9 ~ p + parent_prob + pow, data = outTab, FUN = mean)
# names(out_time) <- c("p", "parent_prob", "pow", "dkw", "senSen", "rptest_ols","rptest_lasso", "mint", "hols")
# 
# 
# 
# 
# 
# nullTab2 <- out_null[which(out_null$pow == 2 & abs(out_null$parent_prob - 1/2) < 1e-5),c(2,1,6, 9:13)]
# altTab2 <- out_alt[which(out_alt$pow == 2 & abs(out_alt$parent_prob - 1/2) < 1e-5),c(6, 9:13)]
# tempTab <- cbind(nullTab2, altTab2)
# 
# for(j in 1:nrow(tempTab)){
#   z <- as.numeric(unlist(tempTab[j, 3:8]))
#   invalidLevel <- z > (.1 + 1.96 * sqrt(.1 * .9 /  500))
#   tempTab[j, 3:8] <- ifelse(invalidLevel ,paste0("\\bftab ", z *1000) , z *1000)
# 
#   z <- as.numeric(unlist(tempTab[j, 9:14]))
#   mz <- max(z * (1-invalidLevel))
#   topPower <- (mz - z * (1-invalidLevel)) <= (1.96 * sqrt(mz * (1-mz) / 500 + z * (1-z) / 500))
#   tempTab[j, 9:14] <- ifelse(!invalidLevel, ifelse(topPower, paste0("\\bftab ", z *1000) , z *1000), "")
# 
# }
# first <- !duplicated(tempTab[[1]])
# tempTab[[1]][!first] <- ""
# tempTab[[1]][first] <- paste0("\\hline\\multirow{", 5, "}{*}{\\STAB{\\rotatebox[origin=c]{90}{\\bftab ", tempTab[[1]][first], "}}}")
# print(xtable::xtable(tempTab, digits = 0),
#       sanitize.text.function = force, include.rownames = F,
#       add.to.row = list(pos = list(-1),
#                         command = c("\\multicolumn{2}{|c|}{} & \\multicolumn{7}{c|}{Size} & \\multicolumn{7}{c|}{Power} \\\\\n")
#                         )
# )
# 
# 
# 
# nullTab5 <- out_null[which(out_null$pow == 5/4 & abs(out_null$parent_prob - 1/2) < 1e-5),c(2,1,6, 9:13)]
# altTab5 <- out_alt[which(out_alt$pow == 5/4 & abs(out_alt$parent_prob - 1/2) < 1e-5),c(6, 9:13)]
# tempTab <- cbind(nullTab5, altTab5)
# 
# for(j in 1:nrow(tempTab)){
#   z <- as.numeric(unlist(tempTab[j, 3:8]))
#   invalidLevel <- z > (.1 + 1.96 * sqrt(.1 * .9 /  500))
#   tempTab[j, 3:8] <- ifelse(invalidLevel ,paste0("\\bftab ", z *1000) , z *1000)
#   
#   z <- as.numeric(unlist(tempTab[j, 9:14])) 
#   mz <- max(z * (1-invalidLevel))
#   topPower <- (mz - z * (1-invalidLevel)) <= (1.96 * sqrt(mz * (1-mz) / 500 + z * (1-z) / 500))
#   tempTab[j, 9:14] <- ifelse(!invalidLevel, ifelse(topPower, paste0("\\bftab ", z *1000) , z *1000), "")
#   
# }
# first <- !duplicated(tempTab[[1]])
# tempTab[[1]][!first] <- ""
# tempTab[[1]][first] <- paste0("\\hline\\multirow{", 5, "}{*}{\\STAB{\\rotatebox[origin=c]{90}{\\bftab ", tempTab[[1]][first], "}}}")
# print(xtable::xtable(tempTab, digits = 0),
#       sanitize.text.function = force, include.rownames = F,
#       add.to.row = list(pos = list(-1),
#                         command = c("\\multicolumn{2}{|c|}{} & \\multicolumn{7}{c|}{Size} & \\multicolumn{7}{c|}{Power} \\\\\n")
#       )
# )
# 
# 
# ############# Main Paper Plots #########
# 
# 
# 
# nullTab5 <- out_null[which(out_null$pow == 5/4 & abs(out_null$parent_prob - 1/2) < 1e-5 &
#                              out_null$p %in% c(10, 20, 45)),c(4, 2,1,6, 9:13)]
# altTab5 <- out_alt[which(out_alt$pow == 5/4 & abs(out_alt$parent_prob - 1/2) < 1e-5 &
#                            out_null$p %in% c(10, 20, 45)), c(6, 9:13)]
# 
# tempTab5 <- cbind(nullTab5, altTab5)
# 
# 
# 
# for(j in 1:nrow(tempTab5)){
#   z <- as.numeric(unlist(tempTab5[j, 4:9]))
#   invalidLevel <- z > (.1 + 1.96 * sqrt(.1 * .9 /  500))
#   tempTab5[j,  4:9] <- ifelse(invalidLevel ,paste0("\\bftab ", round(z *100)) , round(z *100))
#   
#   z <- as.numeric(unlist(tempTab5[j, 10:15]))
#   mz <- max(z * (1-invalidLevel))
#   topPower <- (mz - z * (1-invalidLevel)) <= (1.96 * sqrt(mz * (1-mz) / 500 + z * (1-z) / 500))
#   tempTab5[j, 10:15] <- ifelse(!invalidLevel, ifelse(topPower, paste0("\\bftab ", round(z *100)) , round(z *100)), "")
#   
# }
# 
# 
# first <- !duplicated(tempTab5[[2]])
# tempTab5[[2]][!first] <- ""
# tempTab5[[2]][first] <- paste0("\\multirow{", 3, "}{*}{\\rotatebox[origin=c]{0}{\\bftab ", c("gamma", "laplace", "lognormal", "mixed", "uniform", "weibull"), "}}")
# tempTab5[[1]][-1] <- ""
# tempTab5[[1]][which(first)] <- "\\cline{2-15}"
# tempTab5[[1]][1] <- paste0("\\hline\\multirow{", 18, "}{*}{\\rotatebox[origin=c]{90}{$n \\approx p^{5/4}$}}")
# 
# print(xtable::xtable(tempTab5, digits = 0),
#       sanitize.text.function = force, include.rownames = F,
#       add.to.row = list(pos = list(-1),
#                         command = c("\\multicolumn{3}{|c|}{} & \\multicolumn{6}{c|}{Size} & \\multicolumn{6}{c|}{Power} \\\\\n")
#       )
# )
# 
# nullTab2 <- out_null[which(out_null$pow == 2 & abs(out_null$parent_prob - 1/2) < 1e-5 &
#                              out_null$p %in% c(10, 20, 45)),c(4,2,1,6, 9:13)]
# altTab2 <- out_alt[which(out_alt$pow == 2 & abs(out_alt$parent_prob - 1/2) < 1e-5 &
#                            out_null$p %in% c(10, 20, 45)), c(6, 9:13)]
# 
# tempTab2 <- cbind(nullTab2, altTab2)
#   
# 
# 
# for(j in 1:nrow(tempTab2)){
#   
#   z <- as.numeric(unlist(tempTab2[j, 4:9]))
#   invalidLevel <- z > (.1 + 1.96 * sqrt(.1 * .9 /  500))
#   tempTab2[j,  4:9] <- ifelse(invalidLevel ,paste0("\\bftab ", round(z *100)) , round(z *100))
#   
#   z <- as.numeric(unlist(tempTab2[j, 10:15]))
#   mz <- max(z * (1-invalidLevel))
#   topPower <- (mz - z * (1-invalidLevel)) <= (1.96 * sqrt(mz * (1-mz) / 500 + z * (1-z) / 500))
#   tempTab2[j, 10:15] <- ifelse(!invalidLevel, ifelse(topPower, paste0("\\bftab ", round(z *100)) , round(z *100)), "")
#   
# }
# 
# first <- !duplicated(tempTab2[[2]])
# tempTab2[[2]][!first] <- ""
# tempTab2[[2]][first] <- paste0("\\multirow{", 3, "}{*}{\\rotatebox[origin=c]{0}{\\bftab ", c("gamma", "laplace", "lognormal", "mixed", "uniform", "weibull"), "}}")
# tempTab2[[1]][-1] <- ""
# tempTab2[[1]][which(first)] <- "\\cline{2-15}"
# tempTab2[[1]][1] <- paste0("\\hline\\multirow{", 18, "}{*}{\\rotatebox[origin=c]{90}{$n = p^{2}$}}")
# 
# print(xtable::xtable(tempTab2, digits = 0),
#       sanitize.text.function = force, include.rownames = F,
#       add.to.row = list(pos = list(-1),
#                         command = c("\\multicolumn{3}{|c|}{} & \\multicolumn{6}{c|}{Size} & \\multicolumn{6}{c|}{Power} \\\\\n")
#       )
# )
# 
# 
# 
# timeTab <- cbind(out_time[which(out_time$pow == 5/4 & abs(out_time$parent_prob - 1/2) < 1e-5 &
#                             out_time$p %in% c(10, 20, 45)), c(1, 4:9)],
#                  out_time[which(out_time$pow == 2 & abs(out_time$parent_prob - 1/2) < 1e-5 &
#                                   out_time$p %in% c(10, 20, 45)), c(4:9)])
# numformat <- function(val) { sub("^(-?)0.", "\\1.", sprintf("%.3f", val)) }
# 
# 
# 
# 
# timeTab[, -1] <- ifelse(c(unlist(timeTab[, -1])) < 1, c(unlist(sapply(timeTab[,-1], numformat))), round(c(unlist(timeTab[, -1])), 0))
# print(xtable::xtable(timeTab), include.rownames = F )
# 
