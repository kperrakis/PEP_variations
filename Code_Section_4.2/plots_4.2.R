resultsS1 <- dget('results_g_independent.txt')
resultsS2 <- dget('results_g_correlated.txt')

# -------------------------------------------------
# Figure 8 - line-plots of true probs
# -------------------------------------------------

method.names <- c('DR-PEP','CR-PEP','g-prior','hyper-g','hyper-g/n')
nmethods <- 5
nsize0 <- c(30, 50, 100, 500, 750)
nn <- length(nsize0)
color <- c(1:4, 6)
ppoints <- 1:5

filename <- 'ex_grow_p_true_prob.eps'
postscript(filename,  width = 12, height = 6,horizontal = F)
par(mfrow = c(1,2))

k <-
  1; plot(
    resultsS1$prob_true_model[[k]][1:nn], pch = ppoints[k], type = "b", cex.lab =
      0.9, cex.main = 0.95, main = 'Independent covariates',lwd = 2, ylim = c(0,1), xlab =
      'Sample size and number of covariates', ylab = 'Posterior probability of the true model', col =
      color[k], axes = F
  )
axis(2)
axis(
  1, at = 1:nn, labels = paste(
    paste('n=', nsize0, sep = ''), paste(' p=', resultsS1$p, sep = ''), sep =
      '\n'
  )[1:nn], lwd.ticks = 0, cex.axis = 0.8
)
for (k in 2:5) {
  lines(
    resultsS1$prob_true_model[[k]][1:nn], pch = ppoints[k], lwd = 2, col = color[k], lty =
      color[k], type = 'b'
  )
}
legend(
  'topleft', lty = color, pch = ppoints, col = color, legend = method.names[1:5],  bty =
    'n', lwd = 2
)

k <-
  1; plot(
    resultsS2$prob_true_model[[k]][1:nn], pch = ppoints[k], type = "b", cex.lab =
      0.9, cex.main = 0.95, main = 'Highly-correlated covariates', lwd = 2, ylim =
      c(0,1), xlab = 'Sample size and number of covariates', ylab = 'Posterior probability of the true model', col =
      color[k], axes = F
  )
axis(2)
axis(
  1, at = 1:nn, labels = paste(
    paste('n=', nsize0, sep = ''), paste(' p=', resultsS2$p, sep = ''), sep =
      '\n'
  )[1:nn], lwd.ticks = 0, cex.axis = 0.8
)
for (k in 2:5) {
  lines(
    resultsS2$prob_true_model[[k]][1:nn], pch = ppoints[k], lwd = 2, col = color[k], lty =
      color[k], type = 'b'
  )
}
legend(
  'topleft', lty = color, pch = ppoints, col = color, legend = method.names[1:5],  bty =
    'n', lwd = 2
)

graphics.off()

##########################
nsize0 <- c(30, 50, 100, 500, 750)
nn <- length(nsize0)n
names.n <-
  paste('(n=', nsize0, ", p=", resultsS1$p, ")", sep = '')[1:nn]
nm <- 5
##########################

# -------------------------------------------------
# Figure 9 - incl probs of X for each sample size
# -------------------------------------------------

filename <- paste('ex3_inc_probs.eps',sep = '')
postscript(filename,  horizontal = F, width = 12, height = 8)

# 4x2
pos <- matrix(
  c(1,1,2,3,4,4,5,5,6,6), nrow = 5, ncol = 2, byrow = TRUE
)
layout(pos,  heights = c(0.05, 0.23, 0.23, 0.23, 0.23))
par(lheight = .8)
par(mar = c(0,0,0,0))
plot.new()
color <- c(1:4, 6)
ppoints <- 14 + color
ppoints[nmethods] <- 25
legend(
  'center', pch = ppoints, col = color, legend = method.names,  bty =
    'n', ncol = 5,xpd = T, cex = 2.0, pt.bg = color[nmethods]
)

par(bty = 'l', mar = c(5, 4, 2, 2) + 0.1)
for (j in 1:nn) {
  pp <- as.vector(resultsS1$inc.probs[[1]][[j]])
  xindex <- 1:length(pp)
  plot(
    xindex - 0.4, pp, pch = ppoints[1], cex = 1.5, ylab = "", xlab = "", axes = F, main = names.n[j], ylim =
      c(0,1), xlim = c(0.4, length(pp) + 0.6), cex.main = 1.5
  )
  axis(2, pos = 0.4)
  Xindex <- parse(text = paste("X[", xindex, "]", sep = ""))
  axis(
    1, at = c(0.4,xindex, length(pp) + 0.6), labels = c("", Xindex, ""), pos =
      0, lwd.ticks = 0, cex.axis = 1.8
  )
  
  for (m in (2:nm)) {
    pp <- as.vector(resultsS1$inc.probs[[m]][[j]])
    points(
      xindex - 0.4 + 0.2 * (m - 1), pp, pch = ppoints[m], cex = 1.5, col = color[m], bg = color[nmethods]
    )
  }
  
  abline(v = xindex + 0.5, lty = 2)
}
mtext(
  'Posterior Inclusion Probability',side = 2,line = -1.5,outer = TRUE,cex =
    1.3
)
mtext(
  'Covariates \n (independent scenario)',side = 1,line = -1.1,outer = TRUE,cex =
    1.3
)

graphics.off()

# -------------------------------------------------
# Figure 10 - incl probs of X for each sample size
# -------------------------------------------------

filename <- paste('ex4_inc_probs.eps',sep = '')
postscript(filename,  horizontal = F, width = 12, height = 10)

# 4x2
pos <- matrix(
  c(1,1,2,3,4,4,5,5,6,6), nrow = 5, ncol = 2, byrow = TRUE
)
layout(pos,  heights = c(0.05, 0.23, 0.23, 0.23, 0.23))
par(lheight = .8)
par(mar = c(0,0,0,0))
plot.new()
legend(
  'center', pch = ppoints, col = color, legend = method.names[1:5],  bty =
    'n', ncol = 5,xpd = T, cex = 2.0, pt.bg = color[nmethods]
)

par(bty = 'l', mar = c(5, 4, 2, 2) + 0.1)
for (j in 1:nn) {
  pp <- as.vector(resultsS2$inc.probs[[1]][[j]])
  xindex <- 1:length(pp)
  plot(
    xindex - 0.4, pp, pch = ppoints[1], cex = 1.5, ylab = "", xlab = "", axes = F, main = names.n[j], ylim =
      c(0,1), xlim = c(0.4, length(pp) + 0.6), cex.main = 1.5
  )
  axis(2, pos = 0.4)
  Xindex <- parse(text = paste("X[", xindex, "]", sep = ""))
  axis(
    1, at = c(0.4,xindex, length(pp) + 0.6), labels = c("", Xindex, ""), pos =
      0, lwd.ticks = 0, cex.axis = 1.8
  )
  
  for (m in (2:nm)) {
    pp <- as.vector(resultsS2$inc.probs[[m]][[j]])
    points(
      xindex - 0.4 + 0.2 * (m - 1), pp, pch = ppoints[m], cex = 1.5, col = color[m], bg = color[nmethods]
    )
  }
  
  abline(v = xindex + 0.5, lty = 2)
}
mtext(
  'Posterior Inclusion Probability',side = 2,line = -1.5,outer = TRUE,cex =
    1.3
)
mtext(
  'Covariates \n (highly-correlated scenario)',side = 1,line = -1.1,outer =
    TRUE,cex = 1.3
)

graphics.off()