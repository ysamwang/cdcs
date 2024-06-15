runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)


run.onceBnb <- function(p, n, distro, bs = 400, parent_prob = 1/3, verbose = F, cutoff = NULL){
  
  dat <- cdcs::rDAG(p, n, parent_prob = parent_prob, lowScale = .8,
                    highScale = 1, edgeVar = n^(-1/10),
                    dist = distro, uniqueTop = T)
  
    
    Y <- scale(dat$Y)
    
    
    outlingamDirect <- causalXtreme::direct_lingam_search(Y)

    k.list <- c(1:3)
    
    gFunc <- function(x, k){
      if(k == 1){
        return(scale(x^2))
      } else if (k ==2 ){
        return(scale(x^3))
      } else if (k ==3 ){
        return(scale(abs(x)^(2.5) * sign(x)))
      }
    }
      
    
    G <- array(0, dim = c(n , length(k.list), p))
    
    for(j in 1:p){
      for(k in 1:length(k.list)){
        G[ , k, j] <- gFunc(Y[, j], k.list[k])  
      }
    }
    
    time.rec <- system.time(out <- cdcs::brandAndBound(Y, G, bs = bs, withinAgg = 2, aggType = 2, alpha = .1,
                                         pValueAgg = "tippet", verbose = verbose))[3]
    
    rec <- c(sum(out$pValue > .1), all(out[1, -1] == 1:p),
                mean(cdcs::getAncest(out[which(out$pValue > .1),-1])[lower.tri(dat$B)] == 1 ), time.rec, all(outlingamDirect == 1:p))

    names(rec) <- c("size", "cover", "ancest", "time", "pointEst")
    
  return(rec)
}

##################
library(parallel)
library(cdcs)

sample.size <- 500
rep.runs <- 25
n.list <- c(500, 1000, 2500, 5000)
d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace", "mixed")
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), d.list)
### Param grid Size 400


p <- 10
n <- param.grid[runInd, 1]
distro <- param.grid[runInd, 2]

cl <- makeCluster(7)
clusterExport(cl, ls())
out <- t(parSapply(cl, 1:rep.runs, function(x){run.onceBnb(p, n, distro)}))


outTab <- data.frame(p, n, distro, out)


write.csv(outTab, paste("../results/bnb10/bnbRes10_",runInd, ".csv", sep = ""))

stopCluster(cl)

#### Analysis ####

###Timing ###
sample.size <- 10
rep.runs <- 1
p.list <- seq(10, 20, by = 2)

d.list <- c("gamma")
param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list)
### Param grid Size 400

runInd <- 1
outTab <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/timing/timingRes_new_large",runInd, ".csv", sep = ""))[,-1]
missing <- c()
for(runInd in 2:nrow(param.grid)){
  if(file.exists(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/timing/timingRes_new_large",runInd, ".csv", sep = ""))){
    if(runInd > 40){
      temp <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/timing/timingRes_new_large",runInd, ".csv", sep = ""))
    } else {
      temp <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/timing/timingRes_new_large",runInd, ".csv", sep = ""))[,-1]
    }

    outTab <- rbind(outTab, temp)
  } else {
    missing <- c(missing, runInd)
  }

}

timeP <- aggregate(time ~ p, FUN = median, data = outTab)
# aggregate(time ~ p, FUN = length, data = outTab)
# # 
# 
# 

sample.size <- 500
rep.runs <- 25
n.list <- c(500, 1000, 2500, 5000)
d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace", "mixed")
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), d.list)
### Param grid Size 400
# ### Param grid Size 400
p <- 10


runInd <- 1
outTab <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/bnb10/bnbRes10_",runInd, ".csv", sep = ""))[,-1]
missing <- c()
for(runInd in 2:nrow(param.grid)){
    if(file.exists(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/bnb10/bnbRes10_",runInd, ".csv", sep = ""))){
      temp <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/bnb10/bnbRes10_",runInd, ".csv", sep = ""))[,-1]
      outTab <- rbind(outTab, temp)
    } else {
      missing <- c(missing, runInd)
    }

}
#
resTab <- aggregate(cbind(size, cover, ancest, pointEst, time) ~ n + distro, dat = outTab, FUN = mean)
resTab$size <- resTab$size / factorial(10)

## Number of times conf set is empty
zeroTab <- aggregate(size == 0 ~ n + distro, dat = outTab, FUN = mean)
zeroWide <- reshape(zeroTab, idvar = "n", timevar = "distro", direction = "wide")
colnames(zeroWide) <- c("n", "Gamma", "Laplace", "Lognorm", "Mixed", "Unif", "Weibull")
xtable::xtable(zeroWide)
# 
# 
# 
# 
setEPS()
postscript("~/Dropbox/Apps/Overleaf/Confidence Sets for Causal Discovery/figures/confSets.eps", width = 9, height = 4)
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
     labels = c(".001", ".005", ".01", ".02", ".05", ".1", ".2", ".5"),
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
axis(side = 2, at = seq(0, .5, .1), labels = c("0", ".1", ".2", ".3", ".4", ".5"), las = 1, cex.axis = 1.2)
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
     xaxt = "n", pch = "", ylab = "", xlab = "", log = "y", yaxt = "n", ylim = c(6e2, 80000))
axis(side = 1, at = c(1:6), labels = timeP[,1])
axis(side = 2, at = c(1000, 2000, 5000, 10000, 20000, 50000, 80000),
     labels = c(1000, 2000, 5000, 10000, 20000, 50000, 80000)/ 1000)
mtext("Time (sec x 1000)", side = 3, line = .2, pch = .9)
polygon(c(0,0,7,7), c(1e-5,1e6,1e6, 1e-5), col = "gray95")
abline(h = c(1000, 2000, 5000, 10000, 20000, 50000, 80000), col = "white", lty = 1, lwd = 1.2)
lines(timeP[,2] , pch = 3, col = 3, type = "b")

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')

legend('topright',d.list,pch=1:6,bty='o', col = 1:6, inset = c(0, .07), cex = 1.4)
# 
# 
dev.off()
# 
# 
# tempTab <- cbind(nullTab2, altTab2)
# 
# for(j in 1:nrow(tempTab)){
#   z <- as.numeric(unlist(tempTab[j, 3:8]))
#   invalidLevel <- z > (.1 + 1.96 * sqrt(.1 * .9 /  500))
#   tempTab[j, 3:8] <- ifelse(invalidLevel ,paste0("\\bftab ", z *1000) , z *1000)
#   
#   z <- as.numeric(unlist(tempTab[j, 9:14])) 
#   mz <- max(z * (1-invalidLevel))
#   topPower <- (mz - z * (1-invalidLevel)) < (1.96 * sqrt(mz * (1-mz) / 500 + z * (1-z) / 500))
#   tempTab[j, 9:14] <- ifelse(topPower, paste0("\\bftab ", z *1000) , z *1000)
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
# 
# 
# tempTab <- resTab[, c(1:2,4)]
# 
# outTab <- tempTab[which(tempTab$distro == "gamma"), 3]
# for(i in c("laplace", "lognorm", "mixed", "unif", "weibull")){
#   outTab <- cbind(outTab, tempTab[which(tempTab$distro == i), 3])
# }
# colnames(outTab) <- c("gamma", "laplace", "lognorm", "mixed", "unif", "weibull")
# rownames(outTab) <- unique(tempTab$n)
# 
# xtable::xtable(outTab)
# 
# 
# tempTab <- resTab[, c(1:2,7)]
# 
# outTab <- tempTab[which(tempTab$distro == "gamma"), 3]
# for(i in c("laplace", "lognorm", "mixed", "unif", "weibull")){
#   outTab <- cbind(outTab, tempTab[which(tempTab$distro == i), 3])
# }
# colnames(outTab) <- c("gamma", "laplace", "lognorm", "mixed", "unif", "weibull")
# rownames(outTab) <- unique(tempTab$n)
# 
# print(xtable::xtable(outTab, digits = 0))

