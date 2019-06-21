#---------------------
# PLOTS IN PAPER
#---------------------

results_normal <- dget('results_correlated_normal.txt')
results_laplace <- dget('results_correlated_laplace.txt')
results_a_normal <- dget('results_correlated_a-normal.txt')
results_h_normal <- dget('results_correlated_h-normal1.txt')

source('reset.R')
library(vcd)

# -------------------------------------------------
# Figure 6 - boxplots of true probs
# -------------------------------------------------

nsize0 <- c(30, 50 , 100, 500, 1000)
nn <- length(nsize0)
color <- gray((1:6 + 3) / 9)

method.names <- c('DR-PEP','CR-PEP','g-prior','hyper-g','hyper-g/n')

filename <- 'ex2_true_prob.eps'
postscript(filename,  width = 12, height = 8,horizontal = F)

nmethods <- 5

par(mfrow = c(2,2))

# normal errors

k <-
  1;boxplot(
    results_normal$prob_true_model[[k]][,2:(nn + 1)], at = seq(k,nmethods *
                                                                 nn,nmethods), xlim = c(0,nmethods * nn + 1), names = rep('',nn),
    axes = F, xlab = 'Sample Size', ylab =
      'Posterior probability of the true model', col = color[k],ylim = c(0,1),main =
      '(a) Normal errors'
  )
axis(2)
axis(
  1, at = seq(3,nmethods * nn,nmethods), labels = paste('n=', nsize0,sep =
                                                          ''),lwd.ticks = FALSE, cex.axis = 0.95
)
for (k in 2:nmethods) {
  boxplot(
    results_normal$prob_true_model[[k]][,2:(nn + 1)], at = seq(k,nmethods *
                                                                 nn,nmethods), xlim = c(0,nmethods * nn + 1), names = rep('',nn), col = color[k], add =
      T ,ylim = c(0,1),axes = FALSE
  )
}
axis(1, labels = FALSE, lwd.ticks = FALSE)
abline(v = 5.5, lty = 3)
abline(v = 10.5, lty = 3)
abline(v = 15.5, lty = 3)
abline(v = 20.5, lty = 3)

# laplace errors

k <-
  1;boxplot(
    results_laplace$prob_true_model[[k]][,2:(nn + 1)], at = seq(k,nmethods *
                                                                  nn,nmethods), xlim = c(0,nmethods * nn + 1), names = rep('',nn),
    axes = F, xlab = 'Sample Size', ylab =
      'Posterior probability of the true model', col = color[k],ylim = c(0,1),main =
      '(b) Laplace errors'
  )
axis(2)
axis(
  1, at = seq(3,nmethods * nn,nmethods), labels = paste('n=', nsize0,sep =
                                                          ''),lwd.ticks = FALSE, cex.axis = 0.95
)
for (k in 2:nmethods) {
  boxplot(
    results_laplace$prob_true_model[[k]][,2:(nn + 1)], at = seq(k,nmethods *
                                                                  nn,nmethods), xlim = c(0,nmethods * nn + 1), names = rep('',nn), col = color[k], add =
      T ,ylim = c(0,1),axes = FALSE
  )
}
axis(1, labels = FALSE, lwd.ticks = FALSE)
abline(v = 5.5, lty = 3)
abline(v = 10.5, lty = 3)
abline(v = 15.5, lty = 3)
abline(v = 20.5, lty = 3)

# asymmetric normal errors

k <-
  1;boxplot(
    results_a_normal$prob_true_model[[k]][,2:(nn + 1)], at = seq(k,nmethods *
                                                                   nn,nmethods), xlim = c(0,nmethods * nn + 1), names = rep('',nn),
    axes = F, xlab = 'Sample Size', ylab =
      'Posterior probability of the true model', col = color[k],ylim = c(0,1),main =
      '(c) Asymmetric normal errors'
  )
axis(2)
axis(
  1, at = seq(3,nmethods * nn,nmethods), labels = paste('n=', nsize0,sep =
                                                          ''),lwd.ticks = FALSE, cex.axis = 0.95
)
for (k in 2:nmethods) {
  boxplot(
    results_a_normal$prob_true_model[[k]][,2:(nn + 1)], at = seq(k,nmethods *
                                                                   nn,nmethods), xlim = c(0,nmethods * nn + 1), names = rep('',nn), col = color[k], add =
      T ,ylim = c(0,1),axes = FALSE
  )
}
axis(1, labels = FALSE, lwd.ticks = FALSE)
abline(v = 5.5, lty = 3)
abline(v = 10.5, lty = 3)
abline(v = 15.5, lty = 3)
abline(v = 20.5, lty = 3)

# heteroscedastic normal errors

k <-
  1;boxplot(
    results_h_normal$prob_true_model[[k]][,2:(nn + 1)], at = seq(k,nmethods *
                                                                   nn,nmethods), xlim = c(0,nmethods * nn + 1), names = rep('',nn),
    axes = F, xlab = 'Sample Size', ylab =
      'Posterior probability of the true model', col = color[k],ylim = c(0,1),main =
      '(d) Heteroscedastic normal errors'
  )
axis(2)
axis(
  1, at = seq(3,nmethods * nn,nmethods), labels = paste('n=', nsize0,sep =
                                                          ''),lwd.ticks = FALSE, cex.axis = 0.95
)
for (k in 2:nmethods) {
  boxplot(
    results_h_normal$prob_true_model[[k]][,2:(nn + 1)], at = seq(k,nmethods *
                                                                   nn,nmethods), xlim = c(0,nmethods * nn + 1), names = rep('',nn), col = color[k], add =
      T ,ylim = c(0,1),axes = FALSE
  )
}
axis(1, labels = FALSE, lwd.ticks = FALSE)

abline(v = 5.5, lty = 3)
abline(v = 10.5, lty = 3)
abline(v = 15.5, lty = 3)
abline(v = 20.5, lty = 3)

reset()

legend(
  x = 0.1,y = 0, fill = color, legend = method.names[1:5],  bty = 'n', ncol =
    5,xpd = T
)

graphics.off()

###################################################

ylabel <- character(10)
ylabel[1] <- expression(paste('Covariate ',X[1],sep = ''))
ylabel[2] <- expression(paste('Covariate ',X[2],sep = ''))
ylabel[3] <- expression(paste('Covariate ',X[3],sep = ''))
ylabel[4] <- expression(paste('Covariate ',X[4],sep = ''))
ylabel[5] <- expression(paste('Covariate ',X[5],sep = ''))
ylabel[6] <- expression(paste('Covariate ',X[6],sep = ''))
ylabel[7] <- expression(paste('Covariate ',X[7],sep = ''))
ylabel[8] <- expression(paste('Covariate ',X[8],sep = ''))
ylabel[9] <- expression(paste('Covariate ',X[9],sep = ''))
ylabel[10] <- expression(paste('Covariate ',X[10],sep = ''))

nsize0 <- c(30, 50 , 100, 500, 1000)
color <- gray((1:6 + 3) / 9)

###################################################


# -------------------------------------------------
# Figure 7 - boxplots of incl probs for each X
# -------------------------------------------------

filename <- paste('ex2_inc_probs_normal.eps',sep = '')
postscript(filename,  horizontal = FALSE)
nmethods <- 5

pos <- matrix(c(1,1,2:11), nrow = 6, ncol = 2, byrow = TRUE)
layout(pos, heights = c(0.05, rep(0.95 / 5, 5)))

par(mar = c(0,0,0,0))
plot.new()
legend(
  'center', fill = color, legend = method.names[1:5],  bty = 'n', ncol = 5,xpd =
    T
)

par(bty = 'l', mar = c(5, 4, 1, 2) + 0.1)

for (j in 1:10) {
  k <-
    1;boxplot(
      results_normal$inc.probs[[k]][j,,2:(nn + 1)], at = seq(k,nmethods * nn,nmethods), xlim =
        c(0,nmethods * nn + 1), names = rep('',nn),
      axes = F, xlab = '', ylab = '', col =
        color[k], ylim = c(0,1.0),inset = c(0,-0.2),
      main = ylabel[j]
    )
  axis(2)
  axis(
    1, at = seq(3,nmethods * nn,nmethods), labels = paste('n=', nsize0,sep =
                                                            '')
  )
  for (k in 2:5) {
    boxplot(
      results_normal$inc.probs[[k]][j,,2:(nn + 1)], at = seq(k,nmethods * nn,nmethods), xlim =
        c(0,nmethods * nn + 1), names = rep('',nn), col = color[k], add = T
    )
  }
  abline(v = 5.5, lty = 3)
  abline(v = 10.5, lty = 3)
  abline(v = 15.5, lty = 3)
  abline(v = 20.5, lty = 3)
}
mtext(
  'Posterior Inclusion Probability',side = 2,line = -1,outer = TRUE,cex =
    0.8
)
mtext(
  'Sample Size',side = 1,line = -1,outer = TRUE,cex = 0.8
)

graphics.off()

# -------------------------------------------------
# Figure D.1 - boxplots of incl probs for each X
# -------------------------------------------------

filename <- paste('ex2_inc_probs_laplace.eps',sep = '')
postscript(filename,  horizontal = FALSE)
nmethods <- 5

pos <- matrix(c(1,1,2:11), nrow = 6, ncol = 2, byrow = TRUE)
layout(pos, heights = c(0.05, rep(0.95 / 5, 5)))

par(mar = c(0,0,0,0))
plot.new()
legend(
  'center', fill = color, legend = method.names[1:5],  bty = 'n', ncol = 5,xpd =
    T
)

par(bty = 'l', mar = c(5, 4, 1, 2) + 0.1)

for (j in 1:10) {
  k <-
    1;boxplot(
      results_laplace$inc.probs[[k]][j,,2:(nn + 1)], at = seq(k,nmethods * nn,nmethods), xlim =
        c(0,nmethods * nn + 1), names = rep('',nn),
      axes = F, xlab = '', ylab = '', col =
        color[k], ylim = c(0,1.0),inset = c(0,-0.2),
      main = ylabel[j]
    )
  axis(2)
  axis(
    1, at = seq(3,nmethods * nn,nmethods), labels = paste('n=', nsize0,sep =
                                                            '')
  )
  for (k in 2:5) {
    boxplot(
      results_laplace$inc.probs[[k]][j,,2:(nn + 1)], at = seq(k,nmethods * nn,nmethods), xlim =
        c(0,nmethods * nn + 1), names = rep('',nn), col = color[k], add = T
    )
  }
  abline(v = 5.5, lty = 3)
  abline(v = 10.5, lty = 3)
  abline(v = 15.5, lty = 3)
  abline(v = 20.5, lty = 3)
}
mtext(
  'Posterior Inclusion Probability',side = 2,line = -1,outer = TRUE,cex =
    0.8
)
mtext(
  'Sample Size',side = 1,line = -1,outer = TRUE,cex = 0.8
)

graphics.off()

# -------------------------------------------------
# Figure D.2 - boxplots of incl probs for each X
# -------------------------------------------------

filename <- paste('ex2_inc_probs_a_normal.eps',sep = '')
postscript(filename,  horizontal = FALSE)
nmethods <- 5

pos <- matrix(c(1,1,2:11), nrow = 6, ncol = 2, byrow = TRUE)
layout(pos, heights = c(0.05, rep(0.95 / 5, 5)))

par(mar = c(0,0,0,0))
plot.new()
legend(
  'center', fill = color, legend = method.names[1:5],  bty = 'n', ncol = 5,xpd =
    T
)

par(bty = 'l', mar = c(5, 4, 1, 2) + 0.1)

for (j in 1:10) {
  k <-
    1;boxplot(
      results_a_normal$inc.probs[[k]][j,,2:(nn + 1)], at = seq(k,nmethods * nn,nmethods), xlim =
        c(0,nmethods * nn + 1), names = rep('',nn),
      axes = F, xlab = '', ylab = '', col =
        color[k], ylim = c(0,1.0),inset = c(0,-0.2),
      main = ylabel[j]
    )
  axis(2)
  axis(
    1, at = seq(3,nmethods * nn,nmethods), labels = paste('n=', nsize0,sep =
                                                            '')
  )
  for (k in 2:5) {
    boxplot(
      results_a_normal$inc.probs[[k]][j,,2:(nn + 1)], at = seq(k,nmethods * nn,nmethods), xlim =
        c(0,nmethods * nn + 1), names = rep('',nn), col = color[k], add = T
    )
  }
  abline(v = 5.5, lty = 3)
  abline(v = 10.5, lty = 3)
  abline(v = 15.5, lty = 3)
  abline(v = 20.5, lty = 3)
}
mtext(
  'Posterior Inclusion Probability',side = 2,line = -1,outer = TRUE,cex =
    0.8
)
mtext(
  'Sample Size',side = 1,line = -1,outer = TRUE,cex = 0.8
)

graphics.off()

# -------------------------------------------------
# Figure D.3 - boxplots of incl probs for each X
# -------------------------------------------------

filename <- paste('ex2_inc_probs_h_normal.eps',sep = '')
postscript(filename,  horizontal = FALSE)
nmethods <- 5

pos <- matrix(c(1,1,2:11), nrow = 6, ncol = 2, byrow = TRUE)
layout(pos, heights = c(0.05, rep(0.95 / 5, 5)))

par(mar = c(0,0,0,0))
plot.new()
legend(
  'center', fill = color, legend = method.names[1:5],  bty = 'n', ncol = 5,xpd =
    T
)

par(bty = 'l', mar = c(5, 4, 1, 2) + 0.1)

for (j in 1:10) {
  k <-
    1;boxplot(
      results_h_normal$inc.probs[[k]][j,,2:(nn + 1)], at = seq(k,nmethods * nn,nmethods), xlim =
        c(0,nmethods * nn + 1), names = rep('',nn),
      axes = F, xlab = '', ylab = '', col =
        color[k], ylim = c(0,1.0),inset = c(0,-0.2),
      main = ylabel[j]
    )
  axis(2)
  axis(
    1, at = seq(3,nmethods * nn,nmethods), labels = paste('n=', nsize0,sep =
                                                            '')
  )
  for (k in 2:5) {
    boxplot(
      results_h_normal$inc.probs[[k]][j,,2:(nn + 1)], at = seq(k,nmethods * nn,nmethods), xlim =
        c(0,nmethods * nn + 1), names = rep('',nn), col = color[k], add = T
    )
  }
  abline(v = 5.5, lty = 3)
  abline(v = 10.5, lty = 3)
  abline(v = 15.5, lty = 3)
  abline(v = 20.5, lty = 3)
}
mtext(
  'Posterior Inclusion Probability',side = 2,line = -1,outer = TRUE,cex =
    0.8
)
mtext(
  'Sample Size',side = 1,line = -1,outer = TRUE,cex = 0.8
)

graphics.off()