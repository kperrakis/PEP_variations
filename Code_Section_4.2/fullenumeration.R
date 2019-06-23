library(MASS)
library(mvtnorm)

integer.base.b <- function(x, b = 2) {
  xi <- as.integer(x)
  if (any(is.na(xi) | ((x - xi) != 0)))
    print(list(ERROR = "x not integer", x = x))
  N <- length(x)
  xMax <- max(x)
  ndigits <- (floor(logb(xMax, base = 2)) + 1)
  Base.b <- array(NA, dim = c(N, ndigits))
  for (i in 1:ndigits) {
    #i <- 1
    Base.b[, ndigits - i + 1] <- (x %% b)
    x <- (x %/% b)
  }
  if (N == 1)
    Base.b[1,]
  else
    Base.b
}

prior_model_probabilities2 <- function(k_j, k)
{
  m <- k / 2
  res <-
    gamma(1 + k_j) * gamma((k - m) / m + k - k_j) / gamma(1 + (k - m) / m +
                                                            k) * gamma(1 + (k - m) / m) / gamma((k - m) / m)
  return(res)
}


# full enumeration PCEP = DR-PEP & CRPEP (default is DRPEP)
full_enumeration <-
  function(y,
           x,
           nstar = NULL,
           delta = NULL,
           g = NULL,
           refmodel = "NULL",
           crpep = F)
  {
    n <- nrow(x)
    p <- ncol(x)
    nmodels <- 2 ^ p
    result <- rep(0, nmodels)
    gamma <- integer.base.b(0:(nmodels - 1))
    gamma <- cbind(1, gamma)
    X <- cbind(1, x)
    
    if (is.null(nstar))
      nstar <- n
    
    if (is.null(delta))
      delta <- nstar
    
    if (is.null(g))
      g <- nstar
    
    w <- g / (g + delta)
    
    if (crpep)
      D <- 1
    else
      D <- delta
    
    if (nstar == n)
    {
      Xstar = X
    } else {
      index <- sample(1:n, nstar)
      Xstar <- X[index,]
    }
    
    if (refmodel == "NULL")
    {
      X0star <- matrix(1, nstar, 1)
    } else {
      X0star <- Xstar
    }
    P0star <- X0star %*% solve(t(X0star) %*% X0star) %*% t(X0star)
    L0star <- 1 / D * (diag(nstar) - (g / (g + D)) * P0star)
    
    for (i in 1:nmodels)
    {
      #if(floor(i/100)==i/100) print(i)
      
      Xgstar <-
        as.matrix(Xstar[, gamma[i,] == 1], nrow = nstar)
      Xg <- as.matrix(X[, gamma[i,] == 1], nrow = n)
      Pstar <- Xgstar %*% solve(t(Xgstar) %*% Xgstar) %*% t(Xgstar)
      A <- delta * L0star + w * Pstar
      V_inv <-
        1 / delta * (t(Xgstar) %*% Xgstar * 1 / w - t(Xgstar) %*% solve(A) %*% Xgstar)
      V <- solve(V_inv)
      TAU <- diag(n) + Xg %*% V %*% t(Xg)
      Sigma_tilde <- solve(V_inv + t(Xg) %*% Xg)
      SS <- sum(y ^ 2) - t(y) %*% Xg %*% Sigma_tilde %*% t(Xg) %*% y
      result[i] <-
        lgamma(n / 2) - n / 2 * log(pi) - 1 / 2 * determinant(TAU, logarithm = TRUE)$modulus -
        n / 2 * log(SS)
      
      #return marginal likelihood times prior probability (beta binomial), i.e. log posterior
      
      result[i] <-
        result[i] + log(prior_model_probabilities2((ncol(Xg) - 1), p))
    }
    
    return(result)
    
  }
