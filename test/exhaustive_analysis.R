files <- list.files("test/results/exhaustive/")
dat <- read.csv(paste("test/results/exhaustive/", files[1], sep = ""))
for(f in files[-c(1, length(files))]){
  dat <- rbind(dat, read.csv(paste("test/results/exhaustive/", f, sep = "")))
}

p <- 7
names(dat)
outSize <- cbind(aggregate(totalSetSize10_23 ~ n + p + distro, data = dat, FUN = function(x){mean(x)/factorial(p)}) ,
                 aggregate(totalSetSize10_2 ~ n + p + distro, data = dat, FUN = function(x){mean(x)/factorial(p)})[, 4],
                 aggregate(totalSetSize10_3 ~ n + p + distro, data = dat, FUN = function(x){mean(x)/factorial(p)})[, 4]) 
names(outSize) <-  c("n", "p", "distro", "twoThree", "two", "three")

aggregate(pointEst ~ n + p + distro, data = dat, FUN = mean)

outInc <- cbind(aggregate(totalpValTrue_3 ~ n + p + distro, data = dat, FUN = function(x){mean(x > .1)}) ,
                 aggregate(totalpValTrue_2 ~ n + p + distro, data = dat, FUN = function(x){mean(x > .1)})[, 4],
                 aggregate(totalpValTrue_23 ~ n + p + distro, data = dat, FUN = function(x){mean(x > .1)})[, 4]) 


names(outInc) <-  c("n", "p", "distro", "twoThree", "two", "three")

outEst <- aggregate(pointEst ~ n + p + distro, data = dat, FUN = function(x){mean(x)}) 
plot(outEst[, 4], outSize[, 4])
