library(kableExtra)
library(magrittr)




sample.size <- 1000
rep.runs <- 100
p.list <- c(10, 25, 50, 100, 150, 200)
d.list <- c("gauss", "unif", "lognorm", "gamma", "weibull", "laplace")
param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list)

err.ind <- c()
for(j in 1:dim(param.grid)[1]){
  mtry <- try(read.csv(paste("test/results/levelRes/levelRes_", j, ".csv", sep = "")), silent = T)
  if (class(mtry) == "try-error") {
    err.ind <- c(err.ind, j)
  }
}




j <- 1
dat <- read.csv(paste("test/results/levelRes/levelRes_", j, ".csv", sep = ""))
for(j in 2:dim(param.grid)[1]){
  dat <- rbind(dat, read.csv(paste("test/results/levelRes/levelRes_", j, ".csv", sep = "")))
}


out <- cbind(aggregate(totalpVal ~ p + distro, FUN = function(x){mean(x < .1)}, data = dat),
             aggregate(sumpVal ~ p + distro, FUN = function(x){mean(x < .1)}, data = dat)[ , 3],
             aggregate(minpVal ~ p + distro, FUN = function(x){mean(x < .1)}, data = dat)[ , 3])
             

out[, 4] <- out[, 3] + 2 * sqrt(out[, 3] * (1- out[, 3]) / 1000)
out[, 5] <- out[, 3] - 2 * sqrt(out[, 3] * (1- out[, 3]) / 1000)

names(dat)

hist()

distro <- "laplace"
p <- 150

laplace.list <- which(param.grid[, 1] == p & param.grid[, 2] == distro)
par(mfrow = c(2,5))
for(j in laplace.list){
   hist(read.csv(paste("test/results/levelRes/levelRes_", j, ".csv", sep = ""))[, 5], main = j)
}



j <- 1
dat <- read.csv(paste("test/results/levelRes/levelRes_", j, ".csv", sep = ""))
for(j in 2:dim(param.grid)[1]){
  dat <- rbind(dat, read.csv(paste("test/results/levelRes/levelRes_", j, ".csv", sep = "")))
}


outTotal <- aggregate(totalpVal ~ p + n + distro, FUN = function(x){mean(x < .1)}, data = dat)
outTotal[, 5] <- outTotal[, 4] + 1.96 * sqrt(outTotal[, 4] * (1- outTotal[, 4]) / 1000)
outTotal[, 6] <- outTotal[, 4] - 1.96 * sqrt(outTotal[, 4] * (1- outTotal[, 4]) / 1000)


outSum <- aggregate(sumpVal ~ p + n + distro, FUN = function(x){mean(x < .1)}, data = dat)
outSum[, 5] <- outSum[, 4] + 1.96 * sqrt(outSum[, 4] * (1 - outSum[, 4]) / 1000)
outSum[, 6] <- outSum[, 4] - 1.96 * sqrt(outSum[, 4] * (1- outSum[, 4]) / 1000)


outMin <- aggregate(minpVal ~ p + n + distro, FUN = function(x){mean(x < .1)}, data = dat)
outMin[, 5] <- outMin[, 4] + 1.96 * sqrt(outMin[, 4] * (1 - outMin[, 4]) / 1000)
outMin[, 6] <- outMin[, 4] - 1.96 * sqrt(outMin[, 4] * (1 - outMin[, 4]) / 1000)



#####################
sample.size <- 1000
rep.runs <- 100
p.list <- c(10, 25, 50, 100, 150, 200)
d.list <- c("gauss", "unif", "lognorm", "gamma", "weibull", "laplace")
param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list)

err.ind <- c()
for(j in 1:dim(param.grid)[1]){
  mtry <- try(read.csv(paste("test/results/levelRes/levelRes_", j, ".csv", sep = "")), silent = T)
  if (class(mtry) == "try-error") {
    err.ind <- c(err.ind, j)
  }
}




##################
sample.size <- 1000
rep.runs <- 50

n.list <- c(500, 1000, 1500, 2500)
d.list <- c("unif", "lognorm", "gamma", "weibull", "laplace")
param.grid <- expand.grid(rep(n.list, sample.size / rep.runs), d.list)
### Param grid Size 400

err.ind <- c()
for(j in 1:dim(param.grid)[1]){
  mtry <- try(read.csv(paste("test/results/bnb/bnbRes_", j, ".csv", sep = "")), silent = T)
  if (class(mtry) == "try-error") {
    err.ind <- c(err.ind, j)
  }
}

j <- 1
dat <- read.csv(paste("test/results/bnb/bnbRes_", j, ".csv", sep = ""))
for(j in 2:dim(param.grid)[1]){
  dat <- rbind(dat, read.csv(paste("test/results/bnb/bnbRes_", j, ".csv", sep = "")))
}


resSize <- cbind(aggregate(pointEst ~ n + distro, data = dat, FUN = mean),
                 aggregate(oneFisher3_size ~ n + distro, data = dat, FUN = median)[,3],
                 aggregate(oneFisher23_size ~ n + distro, data = dat, FUN = median)[,3],
             aggregate(oneTippet3_size ~ n + distro, data = dat, FUN = median)[,3],
             aggregate(oneTippet23_size ~ n + distro, data = dat, FUN = median)[,3],
             aggregate(twoFisher3_size ~ n + distro, data = dat, FUN = median)[,3],
             aggregate(twoFisher23_size ~ n + distro, data = dat, FUN = median)[,3],
             aggregate(twoTippet3_size ~ n + distro, data = dat, FUN = median)[,3],
             aggregate(twoTippet23_size ~ n + distro, data = dat, FUN = median)[,3])
             

names(resSize) <- c("n", "Distribution", "Corr. Est", "oneF3", "oneF23",
                    "oneT3", "oneT23",
                    "twoF3", "twoF23",
                    "twoT3", "twoT23")

resSize[, -c(1:3)] <- resSize[, -c(1:3)] /factorial(8)

kable(resSize, format = "latex", digits = 2, booktabs = T, align = "c") %>%
  kable_styling() %>%
  add_header_above(c(" " = 3, "Fisher" = 2, "Tippet" = 2, "Fisher" = 2, "Tippet" = 2)) %>%
  add_header_above(c(" " = 3, "T1" = 4, "T2" = 4)) %>%
  column_spec(1, bold=T) %>%
  collapse_rows(columns = 1:2, latex_hline = "major", valign = "middle")


resSizeThreeOnly <- resSize[, c(2,1,3, 4, 6, 8, 10)]
kable(resSizeThreeOnly, format = "latex", digits = 2, booktabs = T, align = "c") %>%
  kable_styling() %>%
  add_header_above(c(" " = 3, "T1" = 2, "T2" = 2)) %>%
  column_spec(1, bold=T) %>%
  collapse_rows(columns = 1:2, latex_hline = "major", valign = "middle")
  


resTruth <- cbind(aggregate(oneFisher3_truth ~ n + distro, data = dat, FUN = mean),
                 aggregate(oneTippet3_truth ~ n + distro, data = dat, FUN = mean)[,3],
                 aggregate(oneFisher23_truth ~ n + distro, data = dat, FUN = mean)[,3],
                 aggregate(oneTippet23_truth ~ n + distro, data = dat, FUN = mean)[,3],
                 aggregate(twoFisher3_truth ~ n + distro, data = dat, FUN = mean)[,3],
                 aggregate(twoTippet3_truth ~ n + distro, data = dat, FUN = mean)[,3],
                 aggregate(twoFisher23_truth ~ n + distro, data = dat, FUN = mean)[,3],
                 aggregate(twoTippet23_truth ~ n + distro, data = dat, FUN = mean)[,3])

names(resTruth) <- c("n", "distro", "oneF3", "oneT3",
                    "oneF23", "oneT23",
                    "twoF3", "twoT3",
                    "twoF23", "twoT23")

resTruth

kable(resTruth, format = "latex", digits = 2, booktabs = T, align = "c") %>%
  kable_styling() %>%
  add_header_above(c(" " = 2, "Fisher" = 2, "Tippet" = 2, "Fisher" = 2, "Tippet" = 2)) %>%
  add_header_above(c(" " = 2, "T1" = 4, "T2" = 4)) %>%
  column_spec(1, bold=T) %>%
  collapse_rows(columns = 1:2, latex_hline = "major", valign = "middle")




resTime <- cbind(aggregate(oneFisher3_time ~ n + distro, data = dat, FUN = median),
                  aggregate(oneTippet3_time ~ n + distro, data = dat, FUN = median)[,3],
                  aggregate(oneFisher23_time ~ n + distro, data = dat, FUN = median)[,3],
                  aggregate(oneTippet23_time ~ n + distro, data = dat, FUN = median)[,3],
                  aggregate(twoFisher3_time ~ n + distro, data = dat, FUN = median)[,3],
                  aggregate(twoTippet3_time ~ n + distro, data = dat, FUN = median)[,3],
                  aggregate(twoFisher23_time ~ n + distro, data = dat, FUN = median)[,3],
                  aggregate(twoTippet23_time ~ n + distro, data = dat, FUN = median)[,3])

names(resTime) <- c("n", "distro", "F3", "T3",
                     "F23", "T23",
                     "F3", "T3",
                     "F23", "T23")




kable(resTime, format = "latex", digits = 2, booktabs = T, align = "c") %>%
  kable_styling() %>%
  add_header_above(c(" " = 2, "Fisher" = 2, "Tippet" = 2, "Fisher" = 2, "Tippet" = 2)) %>%
  add_header_above(c(" " = 2, "T1" = 4, "T2" = 4)) %>%
  column_spec(1, bold=T) %>%
  collapse_rows(columns = 1:2, latex_hline = "major", valign = "middle")







###################
sample.size <- 1000
rep.runs <- 100
p.list <- c(10, 25, 50, 100, 150, 200)
d.list <- c("gamma", "gammaCor", "laplace", "laplaceCor", "unif", "unifCor")
param.grid <- expand.grid(rep(p.list, sample.size / rep.runs), d.list)

err.ind <- c()
for(j in 1:dim(param.grid)[1]){
  mtry <- try(read.csv(paste("test/results/lmgd/lmgdRes_", j, ".csv", sep = "")), silent = T)
  if (class(mtry) == "try-error") {
    err.ind <- c(err.ind, j)
  }
}


err.ind


j <- 1
dat <- t(read.csv(paste("test/results/lmgd/lmgdRes_", j, ".csv", sep = ""), stringsAsFactors = F))
tempNames <- dat[1, ]
dat <- dat[-1,]
dat <- data.frame(param.grid[j, ], data.matrix(dat), stringsAsFactors = F)
colnames(dat) <- c("p", "distro", tempNames)

for(j in 2:dim(param.grid)[1]){
  newDat <- data.frame(param.grid[j, ], data.matrix(t(read.csv(paste("test/results/lmgd/lmgdRes_", j, ".csv", sep = ""),
                                                        stringsAsFactors = F)[, -1])))
  colnames(newDat) <- c("p", "distro", tempNames)
  dat <- rbind(dat, newDat)
}


dat[, -c(1,2)] <- apply(dat[, -c(1,2)], MAR = 2, as.numeric)


colnames(dat)
res <- cbind(aggregate(dkwMoment1 < .1 ~ p + distro, data = dat, FUN = mean),
             aggregate(dkwMoment2 < .1 ~ p + distro, data = dat, FUN = mean)[, 3],
             aggregate(dkwHSIC < .1 ~ p + distro, data = dat, FUN = mean)[, 3],
             aggregate(ssMoment1 < .1 ~ p + distro, data = dat, FUN = mean)[, 3],
             aggregate(ssMoment2 < .1 ~ p + distro, data = dat, FUN = mean)[, 3],
             aggregate(ssHSIC < .1 ~ p + distro, data = dat, FUN = mean)[, 3])
colnames(res) <- c("p", "distro", tempNames)
res[which(res$p > 100), c(5,8)] <- NA

sqrt(.1 * .9 / 1000) * 2

d.list.title <- c("Gamma", "Corr. Gamma", "Laplace", "Corr. Laplace", "Uniform", "Corr. Uniform")




setEPS()
postscript("test/plots/lmgd.eps", width = 6.5, height = 5)
par(mfcol = c(2,3), mar = c(.5, 1, 1, 0), oma = c(4.5, 1.5, .3, .3))
for(d.ind in 1:length(d.list)){
  
  d <- d.list[d.ind]
  title <- d.list.title[d.ind]
  
  plot_dat <- res[which(res$distro == d), ]
  plot(plot_dat[, 3], ylim = c(0, .32), type = "p", pch = NA, ylab = "Empirical Coverage",
       xlab = "p", xaxt = "n", yaxt = "n")
  
  
  
  if(d.ind %% 2 == 0) {
    axis(1, at = 1:6, labels = p.list)
  }
  if(d.ind < 3){
    axis(2, at = seq(0, .4, by = .1), labels = seq(0, .4, by = .1))
  }
  mtext(title, cex = .8)
  
    
  


  abline(h = seq(0, .4, by = .1), lwd = 1, lty = 3, col = "gray")
  abline(h = .1, lwd = 2, col = "black")
  
  se <- sqrt(plot_dat[,  3] * (1- plot_dat[,  3]) / 1000)
  segments(1:6, plot_dat[, 3] + 1.96 *se, 1:6, plot_dat[,  3] - 1.96 * se, col = "#6699CC")
  lines(plot_dat[,  3], type = "b", col = "#6699CC", pch = 19)
  
  
  se <- sqrt(plot_dat[,  5] * (1- plot_dat[,  5]) / 1000)
  segments(1:6, plot_dat[, 5] + 1.96 *se, 1:6, plot_dat[,  5] - 1.96 * se, col = "#117733")
  lines(plot_dat[,  5], type = "b", col = "#117733", pch = 19)
  
  se <- sqrt(plot_dat[,  6] * (1- plot_dat[,  6]) / 1000)
  segments(1:6, plot_dat[, 6] + 1.96 *se, 1:6, plot_dat[,  6] - 1.96 * se, col = "#DDCC77")
  lines(plot_dat[,  6], type = "b", col = "#DDCC77", pch = 19)
  
  se <- sqrt(plot_dat[,  8] * (1- plot_dat[,  8]) / 1000)
  segments(1:6, plot_dat[, 8] + 1.96 *se, 1:6, plot_dat[,  8] - 1.96 * se, col = "#661100")
  lines(plot_dat[,  8], type = "b", col = "#661100", pch = 19)
  
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')

legend("bottom", ncol = 4, col = c("#6699CC", "#117733", "#DDCC77", "#661100"),
       legend = c("DKW-Moment", "DKW-HSIC", "SS-Moment", "SS-HSIC"),
       pch = 19, lty = 1, xpd =T, inset = c(0, .01))

dev.off()

