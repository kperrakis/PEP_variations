###### generate data

library(mvtnorm)
library(rmutil)
library(twopiece)

included.vars <- c(1, 3, 4)
b <- c(0.3, 0.5, 1.0)
s <- 1.0
s_laplace <- 1 / sqrt(2)
s_1 <- (1 - 0.5)
s_2 <- (1 + 0.5)

# choose error distribution available options are: 'normal', 'laplace', 'a-normal' (asymmetric normal) 
# and 'h-normal' (heteroscedastic normmal)

errors <- 'normal'

sim <- list()
N <- c(15, 30, 50, 100, 500, 1000)
p <- 10
R <- 100
S <- matrix(0.9, p, p)
diag(S) <- rep(1, p)

set.seed(1)

for (k in 1:length(N))
{
  sim[[k]] <- list()
  n <- N[k]
  x <- rmvnorm(n, mean = rep(0, p), S)
  y <- matrix(nrow = n, ncol = R)
  for (i in 1:R)
  {
    if (errors == 'normal') {
      y[, i] <-  x[, included.vars] %*% b  + rnorm(n, 0, s)
    }
    if (errors == 'laplace') {
      y[, i] <-  x[, included.vars] %*% b  + rlaplace(n, 0, s_laplace)
    }
    if (errors == 'a-normal') {
      y[, i] <-  x[, included.vars] %*% b  + rtp3(n, 0, s_1, s_2, rnorm)
    }
    if (errors == 'h-normal') {
      n.errors <- rnorm(n, 0, s)
      h.errors <- exp(x[, included.vars] %*% b) * n.errors
      constant <- sqrt(c(var(n.errors) / var(h.errors)))
      h.errors <- h.errors * constant
      y[, i] <-  x[, included.vars] %*% b  + h.errors
    }
  }
  sim[[k]]$x <- x
  sim[[k]]$y <- y
}

names(sim) <- paste("n=", N, sep = "")
sim$N <- N
sim$true.effects <- included.vars
sim$true.model.index <- sum(2 ^ (p - included.vars)) + 1
sim$true.beta <- b
sim$true.s <- s

###### run models

library(BAS)
source('fullenumeration.R')

true.effects <- sim$true.effects
true.model_index <- sim$true.model.index

p <- ncol(sim[[1]]$x)
nsamples <- ncol(sim[[1]]$y)
Nsizes <- length(sim$N)

nmodels <- 2 ^ p
true.model.matrix <- matrix(0,
                            nrow = nmodels,
                            ncol = p,
                            byrow = T)
true.model.matrix[, true.effects] <- 1
index <- 1:2 ^ p

gamma_m <- integer.base.b(0:(nmodels - 1))

results <- list()

model_prob.drpep <- array(NA, dim = c(2 ^ p, nsamples, Nsizes))
inc_prob.drpep <- array(NA, dim = c(p, nsamples, Nsizes))
model_prob.crpep <- array(NA, dim = c(2 ^ p, nsamples, Nsizes))
inc_prob.crpep <- array(NA, dim = c(p, nsamples, Nsizes))

incl.probs2 <- array(NA, dim = c(p, nsamples, Nsizes))
incl.probs3 <- array(NA, dim = c(p, nsamples, Nsizes))
incl.probs4 <- array(NA, dim = c(p, nsamples, Nsizes))
dimnames(incl.probs2)[[1]] <- paste('X', 1:p, sep = '')
dimnames(incl.probs3)[[1]] <- paste('X', 1:p, sep = '')
dimnames(incl.probs4)[[1]] <- paste('X', 1:p, sep = '')

prob_true_model0 <- matrix(NA, nrow = nsamples, ncol = Nsizes)
prob_true_model1 <- matrix(NA, nrow = nsamples, ncol = Nsizes)
prob_true_model2 <- matrix(NA, nrow = nsamples, ncol = Nsizes)
prob_true_model3 <- matrix(NA, nrow = nsamples, ncol = Nsizes)
prob_true_model4 <- matrix(NA, nrow = nsamples, ncol = Nsizes)

for (sample_sizes in 1:Nsizes)
{
  print(sample_sizes)
  x.original <- sim[[sample_sizes]]$x
  x <- scale(x.original, center = T, scale = F)
  y100 <- sim[[sample_sizes]]$y
  y100 <- scale(y100, center = T, scale = F)
  
  # NOTE that we center X and y for our method. For the rest of the methods working with y or y_center it is the same
  
  n <- nrow(y100)
  
  # results DRPEP
  
  for (i in 1:nsamples) {
    print(noquote(paste("DRPEP Sample", i)))
    res1 <- full_enumeration(y100[, i], x, g = nrow(x) ^ 2)
    res2 <- res1 - max(res1)
    prob <- exp(res2) / sum(exp(res2))
    model_prob.drpep[, i, sample_sizes] <- prob
    for (j in 1:p)
      inc_prob.drpep[j, i, sample_sizes] <-
      sum(prob[gamma_m[, j] == 1])
    prob_true_model0[i, sample_sizes] <-
      model_prob.drpep[true.model_index, i, sample_sizes]
  }
  
  # results CRPEP
  
  for (i in 1:nsamples) {
    print(noquote(paste("CRPEP Sample", i)))
    res1 <-
      full_enumeration(y100[, i], x, g = nrow(x) ^ 2, crpep = T)
    res2 <- res1 - max(res1)
    prob <- exp(res2) / sum(exp(res2))
    model_prob.crpep[, i, sample_sizes] <- prob
    for (j in 1:p)
      inc_prob.crpep[j, i, sample_sizes] <-
      sum(prob[gamma_m[, j] == 1])
    prob_true_model1[i, sample_sizes] <-
      model_prob.crpep[true.model_index, i, sample_sizes]
  }
  
  # results g-prior + hyper-g + hyper-g/n
  
  current.data <- data.frame(y100[, 1], x)
  names(current.data) <- c('y', paste('X', 1:p, sep = ''))
  
  ranking2 <- numeric(nsamples)
  ranking3 <- numeric(nsamples)
  ranking4 <- numeric(nsamples)
  
  for (k in 1:nsamples) {
    print(noquote(paste("Others Sample", k)))
    current.data$y <- y100[, k]
    res.gprior <-
      bas.lm(
        y ~ . ,
        data = current.data,
        prior = "g-prior",
        alpha = n,
        modelprior = beta.binomial()
      )
    res.hyperg <-
      bas.lm(
        y ~ . ,
        data = current.data,
        prior = "hyper-g",
        alpha = 3,
        modelprior = beta.binomial()
      )
    res.hypergn <-
      bas.lm(
        y ~ . ,
        data = current.data,
        prior = "hyper-g-n",
        alpha = 3,
        modelprior = beta.binomial()
      )
    
    incl.probs2[, k, sample_sizes] <- res.gprior[[1]][-1]
    incl.probs3[, k, sample_sizes] <- res.hyperg[[1]][-1]
    incl.probs4[, k, sample_sizes] <- res.hypergn[[1]][-1]
    
    best.gprior <- t(summary(res.gprior, n.models = nmodels)[,-1])
    best.hyperg <- t(summary(res.hyperg, n.models = nmodels)[,-1])
    best.hypergn <-
      t(summary(res.hypergn, n.models = nmodels)[,-1])
    
    ranking2[k] <-
      index[apply(true.model.matrix  == best.gprior[, 2:(p + 1)], 1, all)]
    ranking3[k] <-
      index[apply(true.model.matrix  == best.hyperg[, 2:(p + 1)], 1, all)]
    ranking4[k] <-
      index[apply(true.model.matrix  == best.hypergn[, 2:(p + 1)], 1, all)]
    
    prob_true_model2[k, sample_sizes] <-
      best.gprior[ranking2[k], p + 3]
    prob_true_model3[k, sample_sizes] <-
      best.hyperg[ranking3[k], p + 3]
    prob_true_model4[k, sample_sizes] <-
      best.hypergn[ranking4[k], p + 3]
  }
}

# end of for over the 5 methods for generating data

results$inc.probs <-
  list(inc_prob.drpep,
       inc_prob.crpep,
       incl.probs2,
       incl.probs3,
       incl.probs4)
results$prob_true_model <-
  list(
    prob_true_model0,
    prob_true_model1,
    prob_true_model2,
    prob_true_model3,
    prob_true_model4
  )

###### export results for plots

dput(results,
     paste(getwd(), '/results_correlated_', errors, '.txt', sep = ''))
