## Fama and French 12 industry data
## Accessed Mar 2024 from Kenneth French's data library 
## https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/12_Industry_Portfolios_daily_CSV.zip
##

# load the data and subset to data from 2019-2023
famaData <- read.csv("data/fama12_23.csv")
Y <- famaData[which(substr(famaData[,1], 1, 4) >= 2019), -1]
# scale and center the data
Y <- scale(Y)


#### Compute the confidence set ####
## This takes roughly 
n <- nrow(Y)
p <- ncol(Y)

# Form the test functions
set.seed(100)
J <- 2
G1 <- array(0, dim = c(n, 4, p) )
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

G4 <- abind::abind(G1, G2, along = 2)

# Compute the 95% confidence set 
out <- cdcs::brandAndBound(Y, G4, bs = 800, aggType = 3, alpha = .05,
                            pValueAgg = "tippet")

# Size of 95% confidence set 
factorial(12) / nrow(out4[which(out[,1] >= .05), ])
table(out[which(out[,1] >= .05) ,2])

# Size of 90% confidence set 
factorial(12) / nrow(out[which(out[,1] >= .1), ])
table(out[which(out4[,1] >= .05), 2])


# write.csv(out4, "data/fama_max_res_4_2019_23.csv", row.names = F)



#### Calculating the Frechet Mean ####

# function to calculate the distance between orderings
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

## read in the previously saved output from 95% confidence set computed above 
fama2019 <- read.csv("data/fama_max_res_4_2019_23.csv")
fama2019 <- fama2019[which(fama2019[,1] >= .05), ]
fama2019 <- as.matrix(data.frame(fama2019[, -1]))

## calculate the sum of squared distances for each ordering to every other ordering
## in the confidence set
sumSq <- rep(0, nrow(fama2019))
for(i in 1:(nrow(fama2019))){
  temp <- 0
  for(j in setdiff(1:nrow(fama2019), i)){
    temp <- temp + orderingComp(fama2019[i,], fama2019[j, ])^2
  }

  sumSq[i] <- temp

  if(i %% 100 == 0){
    cat(i)
    cat("/")
    cat(nrow(fama2019))
    cat(" ")
    cat(i / nrow(fama2019))
    cat("\n")
  }
}

# write.csv(sumSq, "data/sumSq_max_4_2019_23_95.csv", row.names = F)


#### Make Plots ####

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

# pdf("~/Dropbox/Apps/Overleaf/Confidence Sets for Causal Discovery/figures/frechetMean_max23.pdf",
    # width = 8, height = 3)
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

# dev.off()


### Get info for text ###

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

#Model selection
#[0, 0]
#[0.268213464470642, 0.413489013608875]
#[0.979916461717943, 1.09345786280352]


### Appendix A: getting adj sets per edge ###
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

# Reorder the matrices according to the frechet mean ordering
lengthMat.rearrange <- lengthMat[colnames(Y)[frechet_mean], colnames(Y)[frechet_mean]]
adjMat.rearrange <- adjMat[colnames(Y)[frechet_mean], colnames(Y)[frechet_mean]]
rownames(lengthMat.rearrange) <- colnames(lengthMat.rearrange) <- c("Utl", "Enrg", "Whl", "Drb", "Fin",
                                                          "Hlth","NoDur", "Tel", "Chem","BsEq", "Mfg","Other")

rownames(adjMat.rearrange) <- colnames(adjMat.rearrange) <- c("Utl", "Enrg", "Whl", "Drb", "Fin",
                                                                    "Hlth","NoDur", "Tel", "Chem","BsEq", "Mfg","Other")


library(plot.matrix)
# setEPS()
# postscript("~/Dropbox/Apps/Overleaf/Confidence Sets for Causal Discovery/figures/ciHeat_max23.eps",
#            width = 8, height = 3)
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
# dev.off()


