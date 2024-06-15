runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)



run.onceTime <- function(p, n, distro, bs = 200, parent_prob = 1/3, verbose = F,
                        cutoff = NULL){
  
  dat <- cdcs::rDAG(p, n, parent_prob = parent_prob, lowScale = .8,
                    highScale = 1, edgeVar = n^(-1/10),
                    dist = distro, uniqueTop = T)
  
    
    Y <- scale(dat$Y)
    
    
    outlingamDirect <- causalXtreme::direct_lingam_search(Y)

    
    thresh <- function(x, cutoff){
      sign(x) * pmin(abs(x), cutoff)
    }
    
    fy <- function(x, method = 1){
      cutoff <- log(length(x)) * 3
      
      
      if( method == .5){
        return(thresh(sign(x) * abs(x)^(.5), cutoff))
        
      } else if( method == 1){
        return(thresh(abs(x), cutoff))
        
      } else if (method == 1.5){
        return(thresh(sign(x) * abs(x)^(1.5), cutoff))
        
      } else if (method == 2){
        return(thresh(x^2, cutoff))
        
      } else if (method == 2.5){
        return(thresh(sign(x) * x^2, cutoff))
        
      } else if (method == 3){
        return(thresh(sign(x) * abs(x)^(2.5), cutoff))
        
      } else if (method == 3.5){
        return(thresh(x^3, cutoff))
        
      }
      
    }
    
    
    ## Fill up the G array
    G <- array(0, dim = c(n, 7, p))
    for(pp in 1:p){
        G[, , pp] <- cbind(fy(Y[,pp], method = .5),
                           fy(Y[,pp], method = 1),
                           fy(Y[,pp], method = 1.5),
                           fy(Y[,pp], method = 2),
                           fy(Y[,pp], method = 2.5),
                           fy(Y[,pp], method = 3),
                           fy(Y[,pp], method = 3.5)) 
    }
    
    try.list <- list(c(4,6, 7))
    
    
    rec <- c()
    for(i in 1:length(try.list)){
      out_time <- system.time(
      out <- cdcs::brandAndBound(Y, G[, try.list[[i]], , drop = F], bs = bs, withinAgg = 2, aggType = 2, alpha = .1,
                                         pValueAgg = "tippet", verbose = verbose))
      temp <- c(sum(out$pValue > .1), all(out[1, -1] == 1:p),
                mean(cdcs::getAncest(out[which(out$pValue > .1),-1])[lower.tri(dat$B)] == 1 ), out_time[3])
      
      names(temp) <- paste(paste(try.list[[i]], collapse = ""), "_", c("size", "cover", "ancest", "time"), sep = "")
      rec <- c(rec, temp)
    }
    rec <- c(rec, all(outlingamDirect == 1:p))
    names(rec)[length(rec)] <- "PointEst"
  
  return(rec)
}

##################
library(cdcs)


sample.size <- 20
rep.runs <- 1

p.list <- c(10, 11, 12, 13)

d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace")
param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list)
### Param grid Size 400


 
n <- 5000
p <- param.grid[runInd, 1]
distro <- param.grid[runInd, 2]


iters <- 5
out <- t(sapply(1:iters, function(x){run.onceTime(p, n, distro)}))


outTab <- data.frame(rep(p, iters), rep(n, iters), rep(distro, iters), out)

colnames(outTab) <- c("p","n", "distro", c("size", "truth", "time", "ancest", "pointEst"))


write.csv(outTab, paste("../results/timing/timingResNew_",runInd, ".csv", sep = ""))


#### Analysis ####
# sample.size <- 250
# rep.runs <- 50
# 
# n.list <- c(1000, 2500, 5000)
# d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace")
# param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), d.list)
# ### Param grid Size 400
# 
# i <- 1
# dat1 <- data.frame(read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/bnbNew/bnb3/bnbRes23_agg23_", i, ".csv", sep = ""))[,-1])
# dat2 <- data.frame(read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/bnbNew/bnb3/bnbRes234_", i, ".csv", sep = ""))[,-1])
# dat3 <- data.frame(read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/bnbNew/bnb3/bnbRes2345_", i, ".csv", sep = ""))[,-1])
# 
# 
# for(i in 21:nrow(param.grid)){
#   dat1 <- rbind(dat1, data.frame(read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/bnbNew/bnb3/bnbRes23_agg23_", i, ".csv", sep = ""))[,-1]))
#   dat2 <- rbind(dat2, data.frame(read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/bnbNew/bnb3/bnbRes234_", i, ".csv", sep = ""))[,-1]))
#   
#   dat3 <- rbind(dat3,
#                 data.frame(read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/bnbNew/bnb3/bnbRes2345_", i, ".csv", sep = ""))[,-1]))
#   
# }
# 
# 
# 
# colnames(dat1) <- colnames(dat2) <- colnames(dat3) <- c("p", "n", "distro" ,paste(rep(c("twoFisher", "twoTippet", "threeFisher", "threeTippet"), each = 4),
#       rep(c("size", "truth", "time", "ancest"), times = 4), sep = "_"), "pointEst")
# 
# #
# res1 <- cbind(aggregate(twoFisher_size ~ p + n + distro, data = dat1, FUN = mean),
#              aggregate(twoTippet_size ~ p + n + distro, data = dat1, FUN = mean)[,4],
#              aggregate(threeFisher_size ~ p + n + distro, data = dat1, FUN = mean)[,4],
#              aggregate(threeTippet_size ~ p + n + distro, data = dat1, FUN = mean)[,4],
#              aggregate(twoFisher_truth ~ p + n + distro, data = dat1, FUN = mean)[,4],
#              aggregate(twoTippet_truth ~ p + n + distro, data = dat1, FUN = mean)[,4],
#              aggregate(threeFisher_truth ~ p + n + distro, data = dat1, FUN = mean)[,4],
#              aggregate(threeTippet_truth ~ p + n + distro, data = dat1, FUN = mean)[,4],
#              aggregate(twoFisher_ancest ~ p + n + distro, data = dat1, FUN = mean)[,4],
#              aggregate(twoTippet_ancest ~ p + n + distro, data = dat1, FUN = mean)[,4],
#              aggregate(threeFisher_ancest ~ p + n + distro, data = dat1, FUN = mean)[,4],
#              aggregate(threeTippet_ancest ~ p + n + distro, data = dat1, FUN = mean)[,4],
#              aggregate(pointEst ~ p + n + distro, data = dat1, FUN = mean)[,4])
# 
# 
# res2 <- cbind(aggregate(twoFisher_size ~ p + n + distro, data = dat2, FUN = mean),
#               aggregate(twoTippet_size ~ p + n + distro, data = dat2, FUN = mean)[,4],
#               aggregate(threeFisher_size ~ p + n + distro, data = dat2, FUN = mean)[,4],
#               aggregate(threeTippet_size ~ p + n + distro, data = dat2, FUN = mean)[,4],
#               aggregate(twoFisher_truth ~ p + n + distro, data = dat2, FUN = mean)[,4],
#               aggregate(twoTippet_truth ~ p + n + distro, data = dat2, FUN = mean)[,4],
#               aggregate(threeFisher_truth ~ p + n + distro, data = dat2, FUN = mean)[,4],
#               aggregate(threeTippet_truth ~ p + n + distro, data = dat2, FUN = mean)[,4],
#               aggregate(twoFisher_ancest ~ p + n + distro, data = dat2, FUN = mean)[,4],
#               aggregate(twoTippet_ancest ~ p + n + distro, data = dat2, FUN = mean)[,4],
#               aggregate(threeFisher_ancest ~ p + n + distro, data = dat2, FUN = mean)[,4],
#               aggregate(threeTippet_ancest ~ p + n + distro, data = dat2, FUN = mean)[,4],
#               aggregate(pointEst ~ p + n + distro, data = dat2, FUN = mean)[,4])
# 
# 
# res3 <- cbind(aggregate(twoFisher_size ~ p + n + distro, data = dat3, FUN = mean),
#               aggregate(twoTippet_size ~ p + n + distro, data = dat3, FUN = mean)[,4],
#               aggregate(threeFisher_size ~ p + n + distro, data = dat3, FUN = mean)[,4],
#               aggregate(threeTippet_size ~ p + n + distro, data = dat3, FUN = mean)[,4],
#               aggregate(twoFisher_truth ~ p + n + distro, data = dat3, FUN = mean)[,4],
#               aggregate(twoTippet_truth ~ p + n + distro, data = dat3, FUN = mean)[,4],
#               aggregate(threeFisher_truth ~ p + n + distro, data = dat3, FUN = mean)[,4],
#               aggregate(threeTippet_truth ~ p + n + distro, data = dat3, FUN = mean)[,4],
#               aggregate(twoFisher_ancest ~ p + n + distro, data = dat3, FUN = mean)[,4],
#               aggregate(twoTippet_ancest ~ p + n + distro, data = dat3, FUN = mean)[,4],
#               aggregate(threeFisher_ancest ~ p + n + distro, data = dat3, FUN = mean)[,4],
#               aggregate(threeTippet_ancest ~ p + n + distro, data = dat3, FUN = mean)[,4],
#               aggregate(pointEst ~ p + n + distro, data = dat3, FUN = mean)[,4])
# 
# #
# #
# colnames(res1) <- colnames(res2) <- colnames(res3) <- c("p", "n", "distro" ,paste(rep(c("twoFisher", "twoTippet", "threeFisher", "threeTippet"), times = 3),
#                                              rep(c("size", "truth", "ancest"), each = 4), sep = "_"), "pointEst")
# 
# 
# res1[, -c(1:3)] <- round(res1[, -c(1:3)], 2)
# res2[, -c(1:3)] <- round(res2[, -c(1:3)], 2)
# res3[, -c(1:3)] <- round(res3[, -c(1:3)], 2)
# 
# res1[, 4:7] / factorial(8)
# res2[, 4:7] / factorial(8)
# res3[, 4:7] / factorial(8)
# 
# res1


# # #### Analysis ####
sample.size <- 50
rep.runs <- 10

p.list <- c(10, 11, 12, 13)

#d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace")
d.list <- c("lognorm", "gamma", "weibull")
param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list)
# ### Param grid Size 400
#
i <- 1
dat1 <- data.frame(read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/timing/timingRes_", i, ".csv", sep = ""))[,-1])


for(i in 2:nrow(param.grid)){
  fn <- paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/timing/timingRes_", i, ".csv", sep = "")
  if(file.exists(fn)){
    dat1 <- rbind(dat1, data.frame(read.csv(fn)[,-1]))
  }
}




res3 <- cbind(aggregate(dat1$X47_size ~ p + n + distro, data = dat1, FUN = mean),
              aggregate(dat1$X47_cover ~ p + n + distro, data = dat1, FUN = mean)[,4],
              aggregate(dat1$X47_ancest ~ p + n + distro, data = dat1, FUN = mean)[,4],
              aggregate(dat1$X47_time ~ p + n + distro, data = dat1, FUN = mean)[,4],
              aggregate(dat1$PointEst ~ p + n + distro, data = dat1, FUN = mean)[,4])
              
names(res3) <- c("p", "n", "distro", "size", "cover", "ancest", "time", "pointEst")

cbind(res3[,(1:3)], round(res3[,4] / factorial(res3$p) * 100, 3), res3[, -c(1:4)])


d.list <- c("gamma", "lognorm", "weibull")
setEPS()
postscript("~/Dropbox/Apps/Overleaf/Confidence Sets for Causal Discovery/figures/confSets.eps", width = 6.5, height = 3)


j <- 1
distro <- d.list[j]
ind <- 7

plot(res3[which(res3[,3] == distro), ind], xlab = "p",
     ylim = c(200, 4000), type = "b", col =j, xaxt = "n",
     ylab = "Time (sec)", main = "Timing",
     lwd = 2, pch = 19, xlim = c(.5, 4.5))
axis(side = 1, at = 1:length(p.list), labels = p.list)


for(j in 2:length(d.list)){
  distro <- d.list[j]
  lines(res3[which(res3[,3] == distro), ind], col = j,
        type = "b", lwd = 2, pch = 19)
  lines(res3[which(res3[,3] == distro), 16], col = j, lty = 3)
  
}
legend("top", col =1:length(d.list), pch = 19, legend = d.list, ncol = 3)
dev.off()


png("~/Dropbox/presentations/2022/ocis/figures/timingRec.png", width = 600, height = 400)
par(mfrow = c(1,2))
j <- 1
distro <- d.list[j]
ind <- 4

plot(res3[which(res3[,3] == distro), ind] / factorial(res3[which(res3[,3] == distro), 1]), xlab = "p",
     ylim = c(0, .03), type = "b", col =j, xaxt = "n",
     ylab = "Prop of orderings in CS", main = "Size of Sets",
     lwd = 2, pch = 19, xlim = c(.5, 4.5))
axis(side = 1, at = 1:length(p.list), labels = p.list)


for(j in 2:length(d.list)){
  distro <- d.list[j]
  lines(res3[which(res3[,3] == distro), ind] / factorial(res3[which(res3[,3] == distro), 1]), col = j,
        type = "b", lwd = 2, pch = 19)
  lines(res3[which(res3[,3] == distro), 16], col = j, lty = 3)
  
}
legend("top", col =1:length(d.list), pch = 19, legend = d.list, ncol = 3, cex = .6)


j <- 1

ind <- 8
plot(res3[which(res3[,3] == distro), ind], xlab = "p",
     ylim = c(0, 1), type = "b", col =j, xaxt = "n",
     ylab = "Prop of Recovery", main = "Prop of Recovery",
     lwd = 2, pch = 19, xlim = c(.5, 4.5))
axis(side = 1, at = 1:length(p.list), labels = p.list)


for(j in 1:length(d.list)){
  distro <- d.list[j]
  lines(res3[which(res3[,3] == distro), ind], col = j,
        type = "b", lwd = 2, pch = 19)
  lines(res3[which(res3[,3] == distro), 16], col = j, lty = 3)
  
}
legend("top", col =1:length(d.list), pch = 19, legend = d.list, ncol = 3, cex = .6)
dev.off()

