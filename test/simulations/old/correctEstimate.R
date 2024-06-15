parent_prob <- .5
p <- 6
n <- 1000
distro <- "gamma"
lowE <- 0
highE <- n^(-1/8)

ss <- 100
rec <- matrix(0,nrow = ss, ncol = 3)
recOrd <- matrix(0,nrow = ss, ncol = 2)

# initDat <- cdcs::rDAG(p, n, parent_prob = parent_prob, lowScale = .5,
#                   highScale = 1, lowEdge = n^(-1/2), highEdge = n^(-1/4), dist = distro, uniqueTop = T)


for(i in 1:ss){
  dat <- cdcs::rDAG(p, n, parent_prob = parent_prob, lowScale = .8,
                    highScale = 1, edgeVar = highE, dist = distro, uniqueTop = T)#, BInput = initDat$B, scalesInput = initDat$scales)
  
  outlingamICA <- pcalg::lingam(dat$Y)
  outlingamDirectD <- cdcs::directLiNGAM(dat$Y, verbose = F, metric = "dhsic")
  outlingamDirectT <- cdcs::directLiNGAM(dat$Y, verbose = F, metric = "tStar")
  
  rec[i, 1] <- sum(abs(outlingamICA$Bpruned[upper.tri(outlingamICA$Bpruned)])) == 0
  rec[i, 2] <- all(outlingamDirectD == 1:p)
  rec[i, 3] <- all(outlingamDirectT == 1:p)
  
  cat(i)
  cat("\n")
  cat(colMeans(rec[1:i, , drop = F]))
  cat("\n")
}

colMeans(rec[1:i, ])

