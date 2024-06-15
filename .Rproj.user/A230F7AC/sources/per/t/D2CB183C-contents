#### Paper Plots #####
#

#### Table 1 ####
# Updated for max stats

sample.size <- 500
rep.runs <- 100
p.list <- c(10, 15, 20, 30, 45)
d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace", "mixed")
pow.list <- c(2, 5/4)
pp.list <- c(1/2)
param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list, pow.list, pp.list)
nrow(param.grid)


runInd <- 1
outTabComp <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/lmgdDAGNew/lmgdDAG_comp_", runInd, ".csv", sep = ""))



for(runInd in 2:nrow(param.grid)){
  temp <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/lmgdDAGNew/lmgdDAG_comp_", runInd, ".csv", sep = ""))
  outTabComp <- rbind(outTabComp, temp)
  
}

out_comp <- aggregate(cbind(null_senSen, null_rptest_ols, null_rptest_lasso,
                            null_mint, null_hols, alt_senSen, alt_rptest_ols, alt_rptest_lasso,
                            alt_mint, alt_hols) < .1 ~ p + distro + pow, data = outTabComp, FUN = mean)

time_comp <- aggregate(cbind(timeA_senSen, timeA_rptest_ols, timeA_rptest_lasso,
                            timeA_mint, timeA_hols) ~ p + pow, data = outTabComp, FUN = mean)


sample.size <- 500
rep.runs <- 50
p.list <- c(10, 15, 20, 30, 45)
d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace", "mixed")
pow.list <- c(2, 5/4)
pp.list <- c(1/2)
param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list, pow.list, pp.list)
nrow(param.grid)

runInd <- 1
outTabDKW <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/lmgdPow/lmgdDAG_dkw_max_", runInd, ".csv", sep = ""))
for(runInd in 2:nrow(param.grid)){
  temp <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/lmgdPow/lmgdDAG_dkw_max_", runInd, ".csv", sep = ""))
  outTabDKW <- rbind(outTabDKW, temp)
  
}

out_DKW <- aggregate(cbind(null_M, alt_M) < .1 ~ TestFunc + p  + distro + pow.i, data = outTabDKW, FUN = mean)
time_DKW <- aggregate(cbind(timeA_M) ~ TestFunc + p  + pow.i, data = outTabDKW, FUN = mean)



### Tab1(a) for timing ###

tFunc <- "comb_TM"
timeTab <- cbind(p = c(10, 20, 45), 
                 cbind(time_DKW$timeA_M[which(time_DKW$pow.i == 1.25 &time_DKW$TestFunc == tFunc & time_DKW$p %in% c(10, 20, 45))],
                 time_comp[which(time_comp$pow == 1.25 & time_comp$p %in% c(10, 20, 45)), -c(1,2)],
                 time_DKW$timeA_M[which(time_DKW$pow.i == 2 &time_DKW$TestFunc == tFunc & time_DKW$p %in% c(10, 20, 45))],
                 time_comp[which(time_comp$pow == 2 & time_comp$p %in% c(10, 20, 45)), -c(1,2)]
                 ) / 1e9)
colnames(timeTab) <- c("p", "W", "S", "RO", "RL", "M", "H", "W", "S", "RO", "RL", "M", "H") 

timeTab[,2:7]/ timeTab[,2]
timeTab[,8:13]/ timeTab[,8]

numformat <- function(val) { sub("^(-?)0.", "\\1.", sprintf("%.3f", val)) }
timeTab[, -1] <- ifelse(c(unlist(timeTab[, -1])) < 1, c(unlist(sapply(timeTab[,-1], numformat))), round(c(unlist(timeTab[, -1])), 0))
print(xtable::xtable(timeTab), include.rownames = F )





### Table 1b for comparison ###
tFunc <- "comb_TM"



perfTab <- rbind(cbind(out_comp[which(out_comp$pow == 1.25), c(3, 2, 1)],
                 W=out_DKW[which(out_DKW$pow.i == 1.25 &out_DKW$TestFunc == tFunc), "null_M"],
                 out_comp[which(out_comp$pow == 1.25),(4:8)],
                 Wa= out_DKW[which(out_DKW$pow.i == 1.25 &out_DKW$TestFunc == tFunc), "alt_M"],
                 out_comp[which(out_comp$pow == 1.25),(9:13)]),
                 
                 cbind(out_comp[which(out_comp$pow == 2), c(3, 2, 1)],
                 W= out_DKW[which(out_DKW$pow.i == 2 &out_DKW$TestFunc == tFunc), "null_M"],
                 out_comp[which(out_comp$pow == 2),(4:8)],
                 Wa=out_DKW[which(out_DKW$pow.i == 2 &out_DKW$TestFunc == tFunc), "alt_M"],
                 out_comp[which(out_comp$pow == 2),(9:13)]))

colnames(perfTab) <- c("pow", "dist", "p", "W", "S", "RO", "RL", "M", "H", "W", "S", "RO", "RL", "M", "H") 


## n = p^(5/4)
temp <- perfTab[which(perfTab$pow == 1.25 & perfTab$p %in% c(10, 20, 45)), ]

multTestAdj <- 18
for(j in 1:nrow(temp)){
  z <- as.numeric(unlist(temp[j, 4:9]))
  invalidLevel <- z > (.1 + abs(qnorm(.05/multTestAdj)) * sqrt(.1 * .9 /  500))
  temp[j,  4:9] <- ifelse(invalidLevel ,paste0("\\bftab ", round(z *100)) , round(z *100))

  z <- as.numeric(unlist(temp[j, 10:15]))
  mz <- max(z * (1-invalidLevel))
  topPower <- (mz - z * (1-invalidLevel)) <= (1.65 * sqrt(mz * (1-mz) / 500 + z * (1-z) / 500))
  temp[j, 10:15] <- ifelse(!invalidLevel, ifelse(topPower, paste0("\\bftab ", round(z *100)) , round(z *100)), "")

}


first <- !duplicated(temp[[2]])
temp[[2]][!first] <- ""
temp[[2]][first] <- paste0("\\multirow{", 3, "}{*}{\\rotatebox[origin=c]{0}{\\bftab ", c("gamma", "laplace", "lognormal", "mixed", "uniform", "weibull"), "}}")
temp[[1]][-1] <- ""
temp[[1]][which(first)] <- "\\cline{2-15}"
temp[[1]][1] <- paste0("\\hline\\multirow{", 18, "}{*}{\\rotatebox[origin=c]{90}{$n \\approx p^{5/4}$}}")

print(xtable::xtable(temp, digits = 0),
      sanitize.text.function = force, include.rownames = F,
      add.to.row = list(pos = list(-1),
                        command = c("\\multicolumn{3}{|c|}{} & \\multicolumn{6}{c|}{Size} & \\multicolumn{6}{c|}{Power} \\\\\n")
      )
)

## n = p^2

temp <- perfTab[which(perfTab$pow == 2 & perfTab$p %in% c(10, 20, 45)), ]


for(j in 1:nrow(temp)){
  z <- as.numeric(unlist(temp[j, 4:9]))
  invalidLevel <- z > (.1 + abs(qnorm(.05/multTestAdj)) * sqrt(.1 * .9 /  500))
  temp[j,  4:9] <- ifelse(invalidLevel ,paste0("\\bftab ", round(z *100)) , round(z *100))
  
  z <- as.numeric(unlist(temp[j, 10:15]))
  mz <- max(z * (1-invalidLevel))
  topPower <- (mz - z * (1-invalidLevel)) <= (1.65 * sqrt(mz * (1-mz) / 500 + z * (1-z) / 500))
  temp[j, 10:15] <- ifelse(!invalidLevel, ifelse(topPower, paste0("\\bftab ", round(z *100)) , round(z *100)), "")
  
}


first <- !duplicated(temp[[2]])
temp[[2]][!first] <- ""
temp[[2]][first] <- paste0("\\multirow{", 3, "}{*}{\\rotatebox[origin=c]{0}{\\bftab ", c("gamma", "laplace", "lognormal", "mixed", "uniform", "weibull"), "}}")
temp[[1]][-1] <- ""
temp[[1]][which(first)] <- "\\cline{2-15}"
temp[[1]][1] <- paste0("\\hline\\multirow{", 18, "}{*}{\\rotatebox[origin=c]{90}{$n \\approx p^{5/4}$}}")

print(xtable::xtable(temp, digits = 0),
      sanitize.text.function = force, include.rownames = F,
      add.to.row = list(pos = list(-1),
                        command = c("\\multicolumn{3}{|c|}{} & \\multicolumn{6}{c|}{Size} & \\multicolumn{6}{c|}{Power} \\\\\n")
      )
)


### Table in Appendix with alternative test functions ###

testFuncAlt <- aggregate(alt_M ~ TestFunc + p + distro + pow.i, data = out_DKW, FUN = mean)
testFuncNull <- aggregate(null_M ~ TestFunc + p + distro + pow.i, data = out_DKW, FUN = mean)
z <- reshape(testFuncAlt, idvar = c("p", "distro", "pow.i"), timevar = "TestFunc", direction = "wide")
z1 <- reshape(testFuncNull, idvar = c("p", "distro", "pow.i"), timevar = "TestFunc", direction = "wide")

testFuncComp <- cbind(z1[, c(3,2,1, 4:ncol(z1))],
                     z[, -c(1:3)])
colnames(testFuncComp) <- c("pow", "dist", "p", rep(c("TM", "RM", "ZM", "MC", "MS", "T", "R", "Z"), times =2)) 

temp <- testFuncComp[which(testFuncComp$pow == 1.25 & altTestFunc$p %in% c(10, 20, 45)), ]

multTestAdj <- 18
for(j in 1:nrow(temp)){
  z <- as.numeric(unlist(temp[j, 4:11]))
  invalidLevel <- z > (.1 + abs(qnorm(.05/multTestAdj)) * sqrt(.1 * .9 /  500))
  temp[j,  4:11] <- ifelse(invalidLevel ,paste0("\\bftab ", round(z *100)) , round(z *100))
  
  z <- as.numeric(unlist(temp[j, 12:19]))
  mz <- max(z * (1-invalidLevel))
  topPower <- (mz - z * (1-invalidLevel)) <= (1.65 * sqrt(mz * (1-mz) / 500 + z * (1-z) / 500))
  temp[j, 12:19] <- ifelse(!invalidLevel, ifelse(topPower, paste0("\\bftab ", round(z *100)) , round(z *100)), "")
  
}


first <- !duplicated(temp[[2]])
temp[[2]][!first] <- ""
temp[[2]][first] <- paste0("\\multirow{", 3, "}{*}{\\rotatebox[origin=c]{0}{\\bftab ", c("gamma", "laplace", "lognormal", "mixed", "uniform", "weibull"), "}}")
temp[[1]][-1] <- ""
temp[[1]][which(first)] <- "\\cline{2-19}"
temp[[1]][1] <- paste0("\\hline\\multirow{", 18, "}{*}{\\rotatebox[origin=c]{90}{$n \\approx p^{5/4}$}}")

print(xtable::xtable(temp, digits = 0),
      sanitize.text.function = force, include.rownames = F,
      add.to.row = list(pos = list(-1),
                        command = c("\\multicolumn{3}{|c|}{} & \\multicolumn{8}{c|}{Size} & \\multicolumn{8}{c|}{Power} \\\\\n")
      )
)


## n = p^2
temp <- testFuncComp[which(testFuncComp$pow == 2 & altTestFunc$p %in% c(10, 20, 45)), ]

multTestAdj <- 18
for(j in 1:nrow(temp)){
  z <- as.numeric(unlist(temp[j, 4:11]))
  invalidLevel <- z > (.1 + abs(qnorm(.05/multTestAdj)) * sqrt(.1 * .9 /  500))
  temp[j,  4:11] <- ifelse(invalidLevel ,paste0("\\bftab ", round(z *100)) , round(z *100))
  
  z <- as.numeric(unlist(temp[j, 12:19]))
  mz <- max(z * (1-invalidLevel))
  topPower <- (mz - z * (1-invalidLevel)) <= (1.65 * sqrt(mz * (1-mz) / 500 + z * (1-z) / 500))
  temp[j, 12:19] <- ifelse(!invalidLevel, ifelse(topPower, paste0("\\bftab ", round(z *100)) , round(z *100)), "")
  
}


first <- !duplicated(temp[[2]])
temp[[2]][!first] <- ""
temp[[2]][first] <- paste0("\\multirow{", 3, "}{*}{\\rotatebox[origin=c]{0}{\\bftab ", c("gamma", "laplace", "lognormal", "mixed", "uniform", "weibull"), "}}")
temp[[1]][-1] <- ""
temp[[1]][which(first)] <- "\\cline{2-19}"
temp[[1]][1] <- paste0("\\hline\\multirow{", 18, "}{*}{\\rotatebox[origin=c]{90}{$n \\approx p^{5/4}$}}")

print(xtable::xtable(temp, digits = 0),
      sanitize.text.function = force, include.rownames = F,
      add.to.row = list(pos = list(-1),
                        command = c("\\multicolumn{3}{|c|}{} & \\multicolumn{8}{c|}{Size} & \\multicolumn{8}{c|}{Power} \\\\\n")
      )
)





#### Figure 2 of paper ####
# 
# Updated: 4/25/24 with max statistics
#

## Get Timing Results ##
sample.size <- 10
rep.runs <- 1
p.list <- seq(10, 20, by = 2)

d.list <- c("gamma")
param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list)
### Param grid Size 400

runInd <- 1
timeTab <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/timing/timingRes_max_",runInd, ".csv", sep = ""))[,-1]
missing <- c()
for(runInd in 2:nrow(param.grid)){
  if(file.exists(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/timing/timingRes_max_",runInd, ".csv", sep = ""))){
    temp <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/timing/timingRes_max_",runInd, ".csv", sep = ""))[,-1]
    timeTab <- rbind(timeTab, temp)
  } else {
    missing <- c(missing, runInd)
  }
}

if(length(missing) > 0){
  temp2 <- data.frame(param.grid[missing,1], rep(10000, 2), param.grid[missing,2],
                      temp[,-c(1:3)])
  temp2[,7] <- 1e8
  colnames(temp2) <- colnames(timeTab)
  
  timeP <- aggregate(time ~ p, FUN = median, data = rbind(timeTab, temp2))  
} else {
  timeP <- aggregate(time ~ p, FUN = median, data = timeTab)  
}

## Get results for everything else ##
sample.size <- 400
rep.runs <- 10
n.list <- c(500, 1000, 2500, 5000)
d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace", "mixed")
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), d.list)

runInd <- 1
outMaxTab <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/bnbMax10/bnbRes_Max_10_",runInd, ".csv", sep = ""))[,-1]
missing <- c()
for(runInd in 1:nrow(param.grid)){
  if(file.exists(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/bnbMax10/bnbRes_Max_10_",runInd, ".csv", sep = ""))){
    temp <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/bnbMax10/bnbRes_Max_10_",runInd, ".csv", sep = ""))[,-1]
    outMaxTab <- rbind(outMaxTab, temp)
  } else {
    missing <- c(missing, runInd)
  }
  
}

#

resTab <- aggregate(cbind(size, cover, ancest, pointEst, time) ~ n + distro + testFunc, dat = outMaxTab, FUN = mean)
# only use combined 
resTab <- resTab[which(resTab$testFunc == "comb"), -3]
resTab$size <- resTab$size / factorial(10)




### Code for Fig 2 ###
setEPS()
postscript("~/Dropbox/Apps/Overleaf/Confidence Sets for Causal Discovery/figures/confSets_max.eps", width = 9, height = 4)
par(mar = c(3, 3, 2, .5), mfrow = c(2,3), oma = c(0, 0, 0,7))

## Point Estimates
i <- 1
plot(resTab[which(resTab$distro == d.list[i]),3], ylim = c(0.05, .3), lty = 1,
     xaxt = "n", pch = "", ylab = "", xlab = "", yaxt = "n")
axis(side = 1, at = c(1:4), labels = c(500, 1000, 2500, 5000))
axis(side = 2, at = seq(0, .3, .05), labels = c("0", ".05", ".1", ".15", ".2", ".25", ".3"), las = 1, cex.axis = 1.2)
mtext("Point Estimate", side = 3, line = .2, pch = .9)
rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     col = "gray95") # Color
abline(h = seq(0, .6, by = .05), col = "white", lty = 1, lwd = 1.2)

for(i in 1:length(d.list)){
  lines(resTab[which(resTab$distro == d.list[i]), 6], type = "b", pch = i, col = i, lty = 1, lwd = 1.5)
}



## Size of Set

i <- 1
plot(resTab[which(resTab$distro == d.list[i]),3], ylim = c(min(resTab[,3])/2, max(resTab[,3])), col = i,
     xaxt = "n", pch = "", ylab = "", xlab = "", yaxt = "n", log = "y")
axis(side = 1, at = c(1:4), labels = c(500, 1000, 2500, 5000))
mtext("Size of Set", side = 3, line = .2, pch = .9)
axis(side = 2, at = c(.001, .005, .01, .02, .05, .1, .25, .5),
     labels = c(".001", ".005", ".01", ".02", ".05", ".1", ".25", ".5"),
     las = 1, cex.axis = 1.2)
# rect(par("usr")[1], par("usr")[3],
#      par("usr")[2], par("usr")[4],
#      col = "gray95") # Color
polygon(c(0,0,5,5), c(1e-5,5,5,1e-5), col = "gray95")

abline(h = c(.001, .005, .01, .02, .05, .1, .2, .5), col = "white", lty = 1, lwd = 1.2)

for(i in 1:length(d.list)){
  lines(resTab[which(resTab$distro == d.list[i]), 3], type = "b", pch = i, col = i, lty = 1, lwd = 1.5)
}


## Timing
i <- 1
plot(resTab[which(resTab$distro == d.list[i]),3], ylim = c(0, max(resTab$time)*1.1), col = i,
     xaxt = "n", pch = "", ylab = "", xlab = "", yaxt = "n")
axis(side = 1, at = c(1:4), labels = c(500, 1000, 2500, 5000))
axis(side = 2, at = seq(0, 2500, 500), labels = c("0", "5", "10", "15", "20", "25"), las = 1, cex.axis = 1.2)
mtext("Time (sec x 100)", side = 3, line = .2, pch = .9)

rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     col = "gray95") # Color
abline(h = seq(0, 2500, by = 500), col = "white", lty = 1, lwd = 1.2)


for(i in 1:length(d.list)){
  lines(resTab[which(resTab$distro == d.list[i]), 7], type = "b", pch = i, col = i, lty = 1, lwd = 1.5)
}



## Coverage
i <- 1
plot(resTab[which(resTab$distro == d.list[i]),3], ylim = c(.8, 1), col = i,
     xaxt = "n", pch = "", ylab = "", xlab = "", yaxt= "n")
axis(side = 1, at = c(1:4), labels = c(500, 1000, 2500, 5000))
axis(side = 2, at = seq(.8, 1, by = .05), labels = c(".8", ".85", ".9", ".95", "1"), las = 1, cex.axis = 1.2)
mtext("Coverage", side = 3, line = .2, pch = .9)

rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     col = "gray95") # Color
abline(h = seq(0, 1, by = .05), col = "white", lty = 1, lwd = 2)
abline(h = .9, col = "red", lty = 1, lwd = 2.5)

for(i in 1:length(d.list)){
  lines(resTab[which(resTab$distro == d.list[i]), 4], type = "b", pch = i, col = i, lty = 1, lwd = 1.5)
}


## Ancestral Relations
i <- 1
plot(resTab[which(resTab$distro == d.list[i]), 5], ylim = c(0, max(resTab[,5]) * 1.1), col = i,
     xaxt = "n", pch = "", ylab = "", xlab = "", yaxt= "n")
axis(side = 1, at = c(1:4), labels = c(500, 1000, 2500, 5000))
axis(side = 2, at = seq(0, .8, .1), labels = c("0", ".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8"), las = 1, cex.axis = 1.2)
mtext("Ancestral Relations", side = 3, line = .2, pch = .9)
rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     col = "gray95") # Color
abline(h = seq(0, 1, by = .1), col = "white", lty = 1, lwd = 1.2)

for(i in 1:length(d.list)){
  lines(resTab[which(resTab$distro == d.list[i]), 5], type = "b", pch = i, col = i, lty = 1, lwd = 1.5)
}






## Timing (p)

plot(timeP[,2], 
     xaxt = "n", pch = "", ylab = "", xlab = "", log = "y", yaxt = "n", ylim = c(1e3, 1.5e5))
axis(side = 1, at = c(1:6), labels = timeP[,1])
axis(side = 2, at = c(1000, 2000, 5000, 10000, 20000, 50000, 1.5e5),
     labels = c(1000, 2000, 5000, 10000, 20000, 50000, 1.5e5)/ 1000)
mtext("Time (sec x 1000)", side = 3, line = .2, pch = .9)
polygon(c(0,0,7,7), c(1e-5,1e6,1e6, 1e-5), col = "gray95")
abline(h = c(1000, 2000, 5000, 10000, 20000, 50000, 1.5e5), col = "white", lty = 1, lwd = 1.2)
lines(timeP[,2] , pch = 3, col = 3, type = "b")

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')

legend('topright',d.list,pch=1:6,bty='o', col = 1:6, inset = c(0, .07), cex = 1.4)

dev.off()

## Number of times conf set is empty (for reviewer comments)
zeroTab <- aggregate(size == 0 ~ n + distro, dat = outMaxTab, FUN = mean)
zeroWide <- reshape(zeroTab, idvar = "n", timevar = "distro", direction = "wide")
colnames(zeroWide) <- c("n", "Gamma", "Laplace", "Lognorm", "Mixed", "Unif", "Weibull")
print(xtable::xtable(zeroWide),include.rownames = F)





#### Figure 3 of paper ####
# confidence sets for non-linear SEMs

sample.size <- 400
rep.runs <- 2
n.list <- c(2500, 5000, 7500, 10000)
d.list <- c("gamma", "laplace")
func.list <- c("poly")
basis.list <- c("poly")
prob.list <- c(1/3)
p.list <- c(7)

param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), d.list, func.list, basis.list,
                          prob.list)


runInd <- 1
missing <- c()

for(runInd in 1:nrow(param.grid)){
  if(file.exists(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/nl_new/nl_exp_easy_max_",runInd, ".csv", sep = ""))){
  } else {
    missing <- c(runInd, missing)
  }
}
missing


runInd <- 1
outTabNL_ez <- data.frame(read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/nl_new/nl_exp_easy_max_", runInd, ".csv", sep = ""))[,-1])

for(runInd in 2:nrow(param.grid)){
    temp <- data.frame(read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/nl_new/nl_exp_easy_max_",runInd, ".csv", sep = ""))[,-1])
    outTabNL_ez <- rbind(outTabNL_ez, temp)

}

outTabNL_ez$cover[which(is.na(outTabNL_ez$cover))] <- F
polyRes <- aggregate(cbind(cover, size)~K + testFunc + distro + n, FUN = mean, data = outTabNL_ez)

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


runInd <- 1
outTabNL_hard <- data.frame(read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/nl_new/nl_exp_hard_", runInd, ".csv", sep = ""))[,-1])

for(runInd in 2:nrow(param.grid)){
  temp <- data.frame(read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/nl_new/nl_exp_hard_",runInd, ".csv", sep = ""))[,-1])
  outTabNL_hard <- rbind(outTabNL_hard, temp)
  
}


outTabNL_hard$cover[which(is.na(outTabNL_hard$cover))] <- F
sigRes <- aggregate(cbind(cover, size)~K + testFunc + distro + n, FUN = mean, data = outTabNL_hard)



setEPS()
postscript("~/Dropbox/Apps/Overleaf/Confidence Sets for Causal Discovery/figures/confSets_nl_max.eps", width = 9, height = 4)
par(mar = c(.6, 3.25, .5, .5), mfcol = c(2, 4), oma = c(4, 1, 3,0))
k <- 2
k.list <- c(2:5)

for(d in 1:length(d.list)){
  plot(polyRes[which(polyRes$testFunc =="comb" & polyRes$K == 2 & polyRes$distro == d.list[d]), "cover"], ylim = c(0, 1), lty = 1,
       xaxt = "n", pch = "", ylab = "", xlab = "", yaxt = "n")
  if(d == 1){
    
  mtext("Coverage", side = 2, line = 3, cex = 1)
  }
  axis(side = 2, at = seq(0, 1, .2), las = 1, cex.axis = 1.2)
  
  mtext(d.list[d], side = 3, line = .5, cex = 1)
  
  rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4],
       col = "gray95") # Color
  abline(h = seq(0, 1, by = .2), col = "white", lty = 1, lwd = 1)
  abline(h = .9, col = "red", lty = 1, lwd = 2)
  
  for(k in k.list){
    
    lines(polyRes[which(polyRes$testFunc =="comb" & polyRes$K == k & polyRes$distro == d.list[d]), "cover"],
          type = "b", pch = k, col = k, lty = 1, lwd = 1.5)
    
  }
  
  k <- 2
  plot(polyRes[which(polyRes$testFunc =="comb" & polyRes$K == 2 & polyRes$distro == d.list[d]), "size"] / factorial(7), ylim = c(1e-6, .25), lty = 1,
       xaxt = "n", pch = "", ylab = "", xlab = "", yaxt = "n", log = "y")
  axis(side = 1, at = c(1:4), labels = polyRes[which(polyRes$testFunc =="comb" & polyRes$K == 2 & polyRes$distro == d.list[d]), "n"], cex.axis = 1.2)
  
  if(d == 1){
  mtext("Size of Set", side = 2, line = 3, cex = 1)
  
 
  }
  axis(side = 2, at = c(1e-6, .0001, .005, .05, .25, .5),
       labels = c(expression(10^-6), expression(10^-4), ".005", ".05", ".25", ".5"),
       las = 1, cex.axis = 1.2)
  polygon(c(0,0,5,5), c(1e-10,5,5,1e-10), col = "gray95")
  
  abline(h = c(.000001, .0001, .005, .05, .25, .5), col = "white", lty = 1, lwd = 1)
  
  k.list <- c(2:5)
  for(k in k.list){
    
    lines(polyRes[which(polyRes$testFunc =="comb" & polyRes$K == k& polyRes$distro == d.list[d]), "size"] / factorial(7),
          type = "b", pch = k, col = k, lty = 1, lwd = 1.5)
    
  }
}

### Hard ###
k.list <- c(20, 40, 60)
for(d in 1:length(d.list)){
  plot(sigRes[which(sigRes$testFunc =="comb" & sigRes$K == 20 & sigRes$distro == d.list[d]), "cover"], ylim = c(0, 1), lty = 1,
       xaxt = "n", pch = "", ylab = "", xlab = "", yaxt = "n")
  if(d == 1){
    
    # mtext("Coverage", side = 2, line = 3.5, cex = 1)
  }
  
  axis(side = 2, at = seq(0, 1, .2), las = 1, cex.axis = 1.2)
  
  mtext(d.list[d], side = 3, line = .5, cex = 1)
  
  rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4],
       col = "gray95") # Color
  abline(h = seq(0, 1, by = .2), col = "white", lty = 1, lwd = 1)
  abline(h = .9, col = "red", lty = 1, lwd = 2)
  
  for(k in 1:length(k.list)){
    
    lines(sigRes[which(sigRes$testFunc =="comb" & sigRes$K == k.list[k] & sigRes$distro == d.list[d]), "cover"],
          type = "b", pch = k, col = k, lty = 1, lwd = 1.5)
    
  }
  
  plot(sigRes[which(sigRes$testFunc =="comb" & sigRes$K == 20 & sigRes$distro == d.list[d]), "size"] / factorial(7), ylim = c(.2, 1), lty = 1,
       xaxt = "n", pch = "", ylab = "", xlab = "", yaxt = "n")
  axis(side = 1, at = c(1:4), labels = sigRes[which(sigRes$testFunc =="comb" & sigRes$K == 20 & sigRes$distro == d.list[d]), "n"], cex.axis = 1.2)
  
  if(d == 1){
    # mtext("Size of Set", side = 2, line = 3.5, cex = 1)
    
    
  }
  
  axis(side = 2, at = seq(.2, 1, by = .2),
       labels = seq(.2, 1, by = .2),
       las = 1, cex.axis = 1.2)
  
  rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4],
       col = "gray95") # Color
  
  abline(h = seq(0, 1, by = .2), col = "white", lty = 1, lwd = 1)
  

  for(k in 1:length(k.list)){
    
    lines(sigRes[which(sigRes$testFunc =="comb" & sigRes$K == k.list[k] & sigRes$distro == d.list[d]), "size"] / factorial(7),
          type = "b", pch = k, col = k, lty = 1, lwd = 1.5)
    
  }
}



par(fig = c(0, .5, 0, 1), oma = c(0, 4, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
mtext("Polynomial", side = 3, line= -1.5)
legend('bottom', legend = c("K=", 2, 3, 4, 5), pch = c(NA,2:5), col = c(NA,2:5), xpd = TRUE, horiz = TRUE, cex = 1, seg.len=1)



par(fig = c(.5, 1, 0, 1), oma = c(0, 4, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
mtext("Sigmoid", side = 3, line= -1.5)
legend('bottom', legend = c("K=", 20, 40, 60), pch = c(NA,1:3), col = c(NA,1:3), xpd = TRUE, horiz = TRUE, cex = 1, seg.len=1)

dev.off()

#### Confidence Intervals which account for model selection ####
# Table 2

sample.size <- 400
rep.runs <- 5

n.list <- c(250, 500, 1000, 2000)
d.list <- c("laplace", "gamma")
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), d.list)

runInd <- 1
outTab <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/ciTest/ciTest_max_",runInd, ".csv", sep = ""))[,-1]
missing <- c()
for(runInd in 2:nrow(param.grid)){
  if(file.exists(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/ciTest/ciTest_max_",runInd, ".csv", sep = ""))){
    temp <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/ciTest/ciTest_max_",runInd, ".csv", sep = ""))[,-1]
    outTab <- rbind(outTab, temp)
  } else {
    missing <- c(missing, runInd)
  }

}


coverage <- aggregate(cbind(cov_ms, length_ms, cov_point, length_point, pointEst) ~ p + n + distro, 
                      data = outTab, FUN = mean)
lengthsMed <- aggregate(cbind(length_ms, length_point) ~ p + n + distro,
                        data = outTab, FUN = median)

covTab <- cbind(coverage[1:4, c(2,4,6, 5,7)], lengthsMed[1:4, 4:5], coverage[1:4, 8],
                coverage[5:8, c(4,6, 5,7)], lengthsMed[5:8, 4:5], coverage[5:8, 8])
numformat <- function(val) { sub("^(-?)0.", "\\1.", sprintf("%.2f", val)) }
covTab[, -c(1)] <- ifelse(c(unlist(covTab[, -c(1)])) < 1,
                            c(unlist(sapply(covTab[,-c(1)], numformat))), round(c(unlist(covTab[, -c(1)])), 1))

# How much larger is proposed vs naive
lengthsMed[1:4, 4]/ lengthsMed[1:4, 5]
lengthsMed[5:8, 4]/ lengthsMed[5:8, 5] 

print(xtable::xtable(covTab), include.rownames = F)


#### Plots for Data Analysis portion ####

orderingComp <- function(o1, o2){
  p <- length(o1)
  temp <- 0
  
  for(i in 1:(p-1)){
    d1 <- o1[(i+1):p]
    d2 <- o2[which(o2 == o1[i]):p]
    temp <- temp + sum(!(d1 %in% d2))
  }
  return(temp)
}


library(plot.matrix)
famaData <- read.csv("~/Dropbox/confSetGraphs/code/rPkg/data_analysis/data/fama12_23.csv")
# Only include data from 2019-2022
Y <- famaData[which(substr(famaData[,1], 1, 4) >= 2019), -1]
# centered (but not scaled) version which will be used to estimate causal effect later
Y.centered <- scale(Y, scale = F)
# centered and scaled which will be used for point estimate
Y <- scale(Y)
# The sum of squares for distance to calculate frechent mean 
sumSq <- as.matrix(read.csv("~/Dropbox/confSetGraphs/code/rPkg/data_analysis/results/sumSq_max_4_2019_23_95.csv"))
# The confidence set with alpha =.1
fama2019 <- as.matrix(read.csv("~/Dropbox/confSetGraphs/code/rPkg/data_analysis/results/fama_max_res_4_2019_23.csv"))




# Point estimate
pointEst <- causalXtreme::direct_lingam_search(Y)
colnames(Y)[pointEst]

# Frechet Mean
frechet_mean <- fama2019[which.min(sumSq), -1]
colnames(Y)[frechet_mean]

### Figure 3 ###

# Make Ancestral matrix
A <- cdcs::getAncest(fama2019)
colnames(A) <- rownames(A) <- colnames(Y)
ZZ <- matrix(0, 12, 12)
ZZ[lower.tri(ZZ, diag = F)] <- 1
symA <- A + t(ZZ - A)
mat <- data.frame(rep(colnames(Y), times = 12), rep(colnames(Y), each = 12), c(symA))
names(mat) <- c("Descendant", "Ancestor", "Proportion")

# Rearrange rows and columns to align with frechet mean
symA.rearrange <- symA[colnames(Y)[frechet_mean], colnames(Y)[frechet_mean]]
colnames(Y)[frechet_mean]

# Make names shorter
rownames(symA.rearrange) <- colnames(symA.rearrange) <- c("Utl", "Enrg", "Whl", "Drb", "Fin",
                                                          "Hlth","NoDur", "Tel", "Chem","BsEq", "Mfg","Other")


# Get Distances from Frechet Mean
distanceDistMean <- rep(0, nrow(fama2019))
distanceDistPoint <- rep(0, nrow(fama2019))
for(i in 1:nrow(fama2019)){
  distanceDistMean[i] <- orderingComp(fama2019[i, -1], fama2019[which.min(sumSq), -1])
  distanceDistPoint[i] <- orderingComp(fama2019[i, -1], pointEst)
}

p1 <- hist(distanceDistMean)  
p2 <- hist(distanceDistPoint) 

pdf("~/Dropbox/Apps/Overleaf/Confidence Sets for Causal Discovery/figures/frechetMean_max23.pdf",
           width = 8, height = 3)
par(oma = c(0, 0, 0, 2), mar = c(5, 4, 2, 2), mfrow = c(1,2))

plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,42), ylim = c(0, 2500), xlab = "", ylab = "", main = "")  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,42), add=T)  # second
mtext("Distance from Mean and Point Est", line = 2, side = 1, cex = .9)
# mtext("Freq", line = 2, side = 2, cex = .7)
# abline(v = orderingComp(fama2019[which.min(sumSq), -1], pointEst), col = "red")
legend("topleft", pch = 15, col = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)), legend = c("Mean", "Pt Est"), ncol = 2,
       cex = .7)

plot(symA.rearrange, col = c("white",RColorBrewer::brewer.pal(n = 9, "Blues")),
     cex.axis = .8, las = 2, xlab = "", ylab = "",
     main = "", key=list(cex.axis=.8, tick=FALSE, side = 4),
     spacing.key = c(1,.5, -.2))
mtext("Descendant", side = 2, line = 3, cex = .9)
mtext("Ancestor", side = 1, line = 3, cex = .9)


dev.off()

# ### Get info for text ###

# Set contains 1/45000 of possible 12! orderings
factorial(12) / nrow(fama2019[which(fama2019[,1] > .05), ])


# Creaste Ancest MAt
ancestMat <- cdcs::getAncest(fama2019[which(fama2019[,1] > .05), ])
colnames(ancestMat) <- rownames(ancestMat) <- names(Y)
ancestMat <- round(ancestMat, 2)
ancestMat[upper.tri(ancestMat, diag = T)] <- ""

xtable::xtable(ancestMat)





colnames(Y)
colnames(Y)[pointEst]
colnames(Y)[frechet_mean]

# Energy (4) onto Durables (2)
treatment <- 4
outcome <- 2
colnames(Y)[c(treatment, outcome)]
out <- cdcs::ci_modSelect(fama2019, treatment, outcome, effectType = "total", alpha = .05, Y.centered)
outReg <- lm(Y.centered[, outcome, drop = F] ~  Y.centered[, pointEst[1:which(pointEst == treatment) ]] -1)
confint(outReg, level = .9)
out
#Model selection
# [-0.134171423060265, 0.374767154913707]
# No Model selection
# -0.09639873  0.02100853



# Chems (5) onto manuf (3)
treatment <- 5
outcome <- 3
colnames(Y)[c(treatment, outcome)]
out <- cdcs::ci_modSelect(fama2019, treatment, outcome, effectType = "total", alpha = .05, Y.centered)
out

#[0, 0]
#[0.268213464470642, 0.413489013608875]
#[0.979916461717943, 1.09345786280352]


### Getting CI's per edge (reviewer comments) ###
lengthMat <- matrix(0, 12, 12)
adjMat <- matrix(0, 12, 12)
for(j in 1:12){
  for(k in 1:12){
    if(j != k){
      temp <- cdcs::ci_modSelect(fama2019, k, j, effectType = "total", alpha = .05, Y.centered)
      lengthMat[j, k] <- temp$le
      adjMat[j, k] <- temp$num
    }
  }
}


colnames(lengthMat) <- rownames(lengthMat) <-
  colnames(adjMat) <- rownames(adjMat) <- colnames(Y)


lengthMat.rearrange <- lengthMat[colnames(Y)[frechet_mean], colnames(Y)[frechet_mean]]
rownames(symA.rearrange) <- colnames(symA.rearrange) <- c("Utl", "Enrg", "Hlth", "NoDr", "Drb",
                                                          "Chem","Whl", "Fin", "Mfg","BusEq", "Tel","Oth")


adjMat.rearrange <- adjMat[colnames(Y)[frechet_mean], colnames(Y)[frechet_mean]]
rownames(symA.rearrange) <- colnames(symA.rearrange) <- c("Utl", "Enrg", "Hlth", "NoDr", "Drb",
                                                          "Chem","Whl", "Fin", "Mfg","BusEq", "Tel","Oth")



library(plot.matrix)
setEPS()
postscript("~/Dropbox/Apps/Overleaf/Confidence Sets for Causal Discovery/figures/ciHeat_max23.eps",
           width = 8, height = 3)
par(oma = c(0, 0, 0, 2), mar = c(5, 4, 2, 5), mfrow = c(1,2))
plot(lengthMat.rearrange, col = c("white",RColorBrewer::brewer.pal(n = 9, "Blues")),
     cex.axis = .8, las = 2, xlab = "", ylab = "",
     main = "", key=list(cex.axis=.8, tick=FALSE, side = 4),
     spacing.key = c(1,.5, -.2))
mtext("Descendant", side = 2, line = 3, cex = .9)
mtext("Ancestor", side = 1, line = 3, cex = .9)
mtext("Length of CI", side = 3, line = 0, cex = .9)

par(mar = c(5, 5, 2, 2))
plot(adjMat.rearrange, col = c("white",RColorBrewer::brewer.pal(n = 9, "Blues")),
     cex.axis = .8, las = 2, xlab = "", ylab = "",
     main = "", key=list(cex.axis=.8, tick=FALSE, side = 4),
     spacing.key = c(1,.5, -.2))
mtext("Descendant", side = 2, line = 3, cex = .9)
mtext("Ancestor", side = 1, line = 3, cex = .9)
mtext("Num of Adjustment Sets", side = 3, line = 0, cex = .9)
dev.off()





