### Plots for Figure 1 ###

set.seed(105)
n <- 5000
Y1 <- runif(n, -1, 1)
Y2 <- .5 * Y1 + runif(n, -1, 1) - 1

ind <- sample(n, size = 300, replace = F)

setEPS()
postscript("~/Dropbox/Apps/Overleaf/Confidence Sets for Causal Discovery/figures/causalDisc1.eps",
           width = 8, height = 2)
par(mfrow = c(1,3), mar = c(2, 4, 4, 0))
plot(Y1[ind], Y2[ind], main = "",
     pch = 19, cex = .4, xlab = "n", ylab = "",
     xaxt = "n", yaxt = "n")
mtext("Raw Data", line = .25)
mtext(expression(Y[1]), side = 1, line = 1)
mtext(expression(Y[2]), side = 2, line = .5)

plot(Y1[ind], lm(Y2 ~ Y1)$res[ind], main = "",
     pch = 19, cex = .4, xlab = "n", ylab = "",
     xaxt = "n", yaxt = "n")
mtext("Correct Model", line = .25)
mtext(expression(Y[1]), side = 1, line = 1)
mtext("Residuals", side = 2, line = .5)


plot(Y2[ind], lm(Y1 ~ Y2)$res[ind], main = "",
     pch = 19, cex = .4, xlab = "n", ylab = "",
     xaxt = "n", yaxt = "n")
mtext("Incorrect Model", line = .25)
mtext(expression(Y[2]), side = 1, line = 1)
mtext("Residuals", side = 2, line = .5)
dev.off()
