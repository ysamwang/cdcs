runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)

set.seed(100 + runInd)

sim.size <- 10

res100 <-  matrix(0, sim.size, ncol = 6)
res1000 <- matrix(0, sim.size, ncol = 6)

n <- 100
for(i in 1:sim.size){
  
  Y1_train <- rgamma(n/2, 1, 1) -1
  Y2_train <- .5 * Y1_train + rgamma(n/2, 1, 1) - 1  
  
  Y1_test <- rgamma(n/2, 1, 1) -1
  Y2_test <- .5 * Y1_test + rgamma(n/2, 1, 1) - 1 
  
  modCor <- RcppArmadillo::fastLm(X = cbind(rep(1, n/2), Y1_train), y = Y2_train)$coef
  modIncor <- RcppArmadillo::fastLm(X = cbind(rep(1, n/2), Y2_train), y = Y1_train)$coef
  
  resCor <- Y2_test - modCor[1] - modCor[2] * Y1_test
  resIncor <- Y1_test - modIncor[1] - modIncor[2] * Y2_test
  
  
  res100[i, 1] <- dHSIC::dhsic.test(resCor, Y1_test, method = "bootstrap")$p
  res100[i, 2] <- TauStar::tauStarTest(resCor, Y1_test)$p
  res100[i, 4] <- dHSIC::dhsic.test(resIncor, Y2_test, method = "bootstrap")$p
  res100[i, 5] <- TauStar::tauStarTest(resIncor, Y2_test)$p

  
  if(i %% 10 == 0){
    cat("i")
    cat(" ")
  }
  
}

n <- 1000
for(i in 1:sim.size){
  
  
  Y1_train <- rgamma(n/2, 1, 1) -1
  Y2_train <- .5 * Y1_train + rgamma(n/2, 1, 1) - 1  
  
  Y1_test <- rgamma(n/2, 1, 1) -1
  Y2_test <- .5 * Y1_test + rgamma(n/2, 1, 1) - 1 
  
  modCor <- RcppArmadillo::fastLm(X = cbind(rep(1, n/2), Y1_train), y = Y2_train)$coef
  modIncor <- RcppArmadillo::fastLm(X = cbind(rep(1, n/2), Y2_train), y = Y1_train)$coef
  
  resCor <- Y2_test - modCor[1] - modCor[2] * Y1_test
  resIncor <- Y1_test - modIncor[1] - modIncor[2] * Y2_test
  
  
  res1000[i, 1] <- dHSIC::dhsic.test(resCor, Y1_test, method = "bootstrap")$p
  res1000[i, 2] <- TauStar::tauStarTest(resCor, Y1_test)$p
  res1000[i, 4] <- dHSIC::dhsic.test(resIncor, Y2_test, method = "bootstrap")$p
  res1000[i, 5] <- TauStar::tauStarTest(resIncor, Y2_test)$p
  
  
  
  if(i %% 10 == 0){
    cat("i")
    cat(" ")
  }
  
}

colnames(res100) <- colnames(res1000) <- c("dhsicSize", "tauSize", "propSize",
                                           "dhsicPower", "tauPower", "propPower")

write.csv(res100, paste("~/confSetsNew/results/ex1/res100_split_bs_", runInd, ".csv", sep = ""), row.names = F)
write.csv(res1000, paste("~/confSetsNew/results/ex1/res1000_split_bs_", runInd, ".csv", sep = ""), row.names = F)


###

runInd <- 1
dat100 <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/ex1/res100_", runInd, ".csv", sep = ""))
dat1000 <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/ex1/res1000_", runInd, ".csv", sep = ""))
dat100ss <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/ex1/res100_split_", runInd, ".csv", sep = ""))
dat1000ss <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/ex1/res1000_split_", runInd, ".csv", sep = ""))

for(runInd in 2:100){
  dat100 <- rbind(dat100,
    read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/ex1/res100_", runInd, ".csv", sep = "")))

  dat1000 <- rbind(dat1000,
                   read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/ex1/res1000_", runInd, ".csv", sep = "")))
  dat100ss <- rbind(dat100ss,
                    read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/ex1/res100_split_", runInd, ".csv", sep = "")))
  dat1000ss <- rbind(dat1000ss,
                     read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/ex1/res1000_split_", runInd, ".csv", sep = "")))
}


a <- .1
tab <- cbind(rbind(c(colMeans(dat100 < a)[1:2], colMeans(dat100ss < a)[1:2], colMeans(dat100 < a)[3]),
            c(colMeans(dat1000 < a)[1:2], colMeans(dat1000ss < a)[1:2], colMeans(dat1000 < a)[3])),
      rbind(c(colMeans(dat100 < a)[4:5], colMeans(dat100ss < a)[4:5], colMeans(dat100 < a)[6]),
      c(colMeans(dat1000 < a)[4:5], colMeans(dat1000ss < a)[4:5], colMeans(dat1000 < a)[6])))
xtable::xtable(tab)




runInd <- 1
dat100 <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/ex1/res100_", runInd, ".csv", sep = ""))
dat1000 <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/ex1/res1000_", runInd, ".csv", sep = ""))
dat100ss <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/ex1/res100_split_bs_", runInd, ".csv", sep = ""))
dat1000ss <- read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/ex1/res1000_split_bs_", runInd, ".csv", sep = ""))

for(runInd in 2:100){
  dat100 <- rbind(dat100,
                  read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/ex1/res100_", runInd, ".csv", sep = "")))
  
  dat1000 <- rbind(dat1000,
                   read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/ex1/res1000_", runInd, ".csv", sep = "")))
  dat100ss <- rbind(dat100ss,
                    read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/ex1/res100_split_bs_", runInd, ".csv", sep = "")))
  dat1000ss <- rbind(dat1000ss,
                     read.csv(paste("~/Dropbox/confSetGraphs/code/rPkg/simResults/ex1/res1000_split_bs_", runInd, ".csv", sep = "")))
}


tab <- cbind(rbind(c(colMeans(dat100 < .1)[1:2], colMeans(dat100ss < .1)[1:2], colMeans(dat100 < .1)[3]),
                   c(colMeans(dat1000 < .1)[1:2], colMeans(dat1000ss < .1)[1:2], colMeans(dat1000 < .1)[3])),
             rbind(c(colMeans(dat100 < .1)[4:5], colMeans(dat100ss < .1)[4:5], colMeans(dat100 < .1)[6]),
                   c(colMeans(dat1000 < .1)[4:5], colMeans(dat1000ss < .1)[4:5], colMeans(dat1000 < .1)[6])))
xtable::xtable(tab)


tab <- cbind(rbind(c(colMeans(dat100 < .1)[1:2], colMeans(dat100ss < .1)[1:2], colMeans(dat100 < .1)[3]),
                   c(colMeans(dat1000 < .1)[1:2], colMeans(dat1000ss < .1)[1:2], colMeans(dat1000 < .1)[3])),
             rbind(c(colMeans(dat100 < .1)[4:5], colMeans(dat100ss < .1)[4:5], colMeans(dat100 < .1)[6]),
                   c(colMeans(dat1000 < .1)[4:5], colMeans(dat1000ss < .1)[4:5], colMeans(dat1000 < .1)[6])))
xtable::xtable(tab)