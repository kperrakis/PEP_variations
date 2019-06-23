library(BAS)
source('fullenumeration.R')

# generating data

included.vars <- c(1,3,4)
b <- c(0.3, 0.5, 1.0)
s <- 1.0

set.seed(1)
sim <- list()
N <- c(30,50,100,500,750)
P <- c(4,5,7,13,15)
R <- 100
k <- length(N)
sim[[k]] <- list()
n <- N[k]
p <- P[k]
x <- matrix(rnorm(n * p), ncol = p, nrow = n)
y <- matrix(nrow = n, ncol = R)
for (i in 1:R)
{
  y[,i] <-  x[,included.vars] %*% b  + rnorm(n, 0, s)
}
sim[[k]]$x <- x
sim[[k]]$y <- y

for (k in 1:(length(N) - 1)) {
  sim[[k]] <- list()
  n <- N[k]
  p <- P[k]
  sim[[k]]$x <- x[1:n, 1:p]
  sim[[k]]$y <- y[1:n,]
}

names(sim) <- paste("n=", N, ", p=", P, sep = "")
sim$N <- N
sim$P <- P
sim$true.effects <- included.vars
sim$true.model.index <- sum(2 ^ (included.vars - 1)) + 1
sim$true.beta <- b
sim$true.s <- s

# running the methods

true.effects <- sim$true.effects

nsamples <- ncol(sim[[1]]$y)
Nsizes <- length(sim$N)

results <- list()
model_prob.drpep <- list()
inc_prob.drpep  <- list()
model_prob.crpep <- list()
inc_prob.crpep  <- list()
incl.probs2 <- list()
incl.probs3 <- list()
incl.probs4 <- list()

nvec <- numeric()
pvec <- numeric()

prob_true_model0 <- matrix(NA, nrow = nsamples, ncol = Nsizes)
prob_true_model1 <- matrix(NA, nrow = nsamples, ncol = Nsizes)
prob_true_model2 <- matrix(NA, nrow = nsamples, ncol = Nsizes)
prob_true_model3 <- matrix(NA, nrow = nsamples, ncol = Nsizes)
prob_true_model4 <- matrix(NA, nrow = nsamples, ncol = Nsizes)

for (q in 1:Nsizes)
{
  print(q)
  x.original <- sim[[q]]$x
  x <- scale(x.original, center = T, scale = F)
  y_c <- sim[[q]]$y
  y_c <- scale(y_c, center = T, scale = F)
  
  # NOTE that we center X and y for our method. For the rest of the methods working with y or y_center it is the same
  
  n <- nrow(y_c)
  p <- ncol(sim[[q]]$x)
  
  true.model_index <- sum(2 ^ (p - true.effects)) + 1
  
  nvec <- c(nvec, n)
  pvec <- c(pvec, p)
  
  index <- 1:2 ^ p
  nmodels <- 2 ^ p
  gamma_m <- integer.base.b(0:(nmodels - 1))
  
  true.model.matrix <-
    matrix(0, nrow = nmodels, ncol = p, byrow = T)
  true.model.matrix[,true.effects] <- 1
  
  model_prob.drpep[[q]] <- array(NA, dim = c(2 ^ p, nsamples))
  inc_prob.drpep[[q]]  <- array(NA, dim = c(p, nsamples))
  model_prob.crpep[[q]] <- array(NA, dim = c(2 ^ p, nsamples))
  inc_prob.crpep[[q]]  <- array(NA, dim = c(p, nsamples))
  
  incl.probs2[[q]] <- array(NA, dim = c(p, nsamples))
  incl.probs3[[q]] <- array(NA, dim = c(p, nsamples))
  incl.probs4[[q]] <- array(NA, dim = c(p, nsamples))
  dimnames(incl.probs2[[q]])[[1]] <- paste('X', 1:p, sep = '')
  dimnames(incl.probs3[[q]])[[1]] <- paste('X', 1:p, sep = '')
  dimnames(incl.probs4[[q]])[[1]] <- paste('X', 1:p, sep = '')
  
  # results DRPEP
  
  for (i in 1:nsamples) {
    print(noquote(paste("DRPEP Sample",i)))
    res1 <- full_enumeration(y_c[,i], x, g = nrow(x) ^ 2)
    res2 <- res1 - max(res1)
    prob <- exp(res2) / sum(exp(res2))
    model_prob.drpep[[q]][,i] <- prob
    for (j in 1:p)
      inc_prob.drpep[[q]][j,i] <- sum(prob[gamma_m[,j] == 1])
    prob_true_model0[i,q] <-
      model_prob.drpep[[q]][true.model_index,i]
  }
  
  # results CRPEP
  
  for (i in 1:nsamples) {
    print(noquote(paste("CRPEP Sample",i)))
    res1 <- full_enumeration(y_c[,i], x, g = nrow(x) ^ 2, crpep = T)
    res2 <- res1 - max(res1)
    prob <- exp(res2) / sum(exp(res2))
    model_prob.crpep[[q]][,i] <- prob
    for (j in 1:p)
      inc_prob.crpep[[q]][j,i] <- sum(prob[gamma_m[,j] == 1])
    prob_true_model1[i,q] <-
      model_prob.crpep[[q]][true.model_index,i]
  }
  
  # results g-prior + hyper-g + hyper-g/n
  
  current.data <- data.frame(y_c[,1], x)
  names(current.data) <- c('y', paste('X', 1:p, sep = ''))
  
  ranking2 <- numeric(nsamples)
  ranking3 <- numeric(nsamples)
  ranking4 <- numeric(nsamples)
  
  for (i in 1:nsamples) {
    print(noquote(paste("Others Sample",i)))
    current.data$y <- y_c[,i]
    res.gprior <-
      bas.lm(
        y ~ . , data = current.data, prior = "g-prior", alpha = n, modelprior =
          beta.binomial()
      )
    res.hyperg <-
      bas.lm(
        y ~ . , data = current.data, prior = "hyper-g", alpha = 3, modelprior =
          beta.binomial()
      )
    res.hypergn <-
      bas.lm(
        y ~ . , data = current.data, prior = "hyper-g-n", alpha = 3, modelprior =
          beta.binomial()
      )
    
    incl.probs2[[q]][,i] <- res.gprior[[1]][-1]
    incl.probs3[[q]][,i] <- res.hyperg[[1]][-1]
    incl.probs4[[q]][,i] <- res.hypergn[[1]][-1]
    
    best.gprior <- t(summary(res.gprior, n.models = nmodels)[,-1])
    best.hyperg <- t(summary(res.hyperg, n.models = nmodels)[,-1])
    best.hypergn <- t(summary(res.hypergn, n.models = nmodels)[,-1])
    
    ranking2[i] <-
      index[apply(true.model.matrix  == best.gprior[,2:(p + 1)],1,all)]
    ranking3[i] <-
      index[apply(true.model.matrix  == best.hyperg[,2:(p + 1)],1,all)]
    ranking4[i] <-
      index[apply(true.model.matrix  == best.hypergn[,2:(p + 1)],1,all)]
    
    prob_true_model2[i,q] <- best.gprior[ranking2[i], p + 3]
    prob_true_model3[i,q] <- best.hyperg[ranking3[i], p + 3]
    prob_true_model4[i,q] <- best.hypergn[ranking4[i], p + 3]
  }
  
  results$inc.probs <-
    list(inc_prob.drpep,inc_prob.crpep, incl.probs2,incl.probs3,incl.probs4)
  results$prob_true_model <-
    list(
      prob_true_model0,prob_true_model1,prob_true_model2,prob_true_model3,prob_true_model4
    )
  results$p <- pvec
  results$n <- nvec
}

# end of for-loop for generating the results

# summarizing and naming the results

methods <- c('DR-PEP', 'CR-PEP', 'g-prior', 'hyper-g', 'hyper-g/n')
nmethods <- length(methods)

for (k in 1:nmethods) {
  results$inc.probs[[k]] <-
    lapply(results$inc.probs[[k]], function (x)
      as.matrix(apply(x,1,mean)))
  results$prob_true_model[[k]] <-
    apply(results$prob_true_model[[k]],2,mean)
}

names(results$inc.probs) <- methods
names(results$prob_true_model) <- methods
sample.sizes <- paste('n=', sim$N[1:Nsizes], sep = '')
for (k in 1:nmethods) {
  names(results$inc.probs[[k]]) <- sample.sizes[1:Nsizes]
  names(results$prob_true_model[[k]])  <-
    c(sample.sizes[1:Nsizes], rep("", 5 - Nsizes))
}

# exporting results

dput(results, paste(getwd(), '/results_g_independent.txt', sep = ''))