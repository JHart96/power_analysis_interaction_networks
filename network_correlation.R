# Simulate random graph data collection. Gamma is parameterised by a, b.
simulate_data <- function (n, p, D, a=10, b=1) {
  # Simulate true interaction rates.
  n <- 20 # Individuals in the network.
  p <- 0.3 # Network density.
  D <- 20 # Mean sampling time per dyad.
  A <- matrix(rbinom(n^2, 1, p) * runif(n^2), n, n)
  A <- A * upper.tri(A)
  A <- A + t(A)
  
  # Simulate observed interaction counts and sampling effort.
  X <- rpois(n^2, as.vector(A) * D)
  X <- matrix(X, n, n)
  X <- X * upper.tri(X)
  X <- X + t(X)
  D <- rgamma(n^2, a, b)
  D <- matrix(D, n, n)
  D <- D * upper.tri(D)
  D <- D + t(D)
  
  list(X=X, D=D, A=A)
}

# Functions for estimating the correlation between an estimated network and its unobserved, true network.
hmean <- function(x) {
  length(x)/sum(1/x)
}

# Uses Equation 3 to estimate the correlation between true and estimated network given key properties of data.
compute_rho <- function(S, I=NULL, mu=NULL, d=NULL) {
  if (is.null(I)) {
    I <- mu * hmean(d)
  }
  (S * sqrt(I))/(sqrt(1 + I * S^2))
}

# Log-likelihood function of negative binomial with gamma parameterisation.
lk_nbinom <- function(par, x, d) {
  if (min(par) < 0) {
    return(-Inf)
  }
  r <- par[1]
  p <- par[2]/(par[2] + d)
  sum(dnbinom(x, r, p, log=TRUE))
}

# Perform numerical maximum likelihood estimation for point estimates of parameters.
numerical_mle <- function(x, d) {
  target <- function(par) lk_nbinom(par, x, d)
  parameters <- optim(c(1, 1), function(par) -target(par))$par
  parameters
}

# Estimates the correlation between true and estimated network given count data x and sampling data d.
estimate_correlation <- function(x, d) {
  parameters <- numerical_mle(x, d)
  
  a_ <- parameters[1]
  b_ <- parameters[2]
  
  S_ <- 1/sqrt(a_)
  mu_ <- a_/b_
  
  I_ <- mu_ * hmean(d)
  
  rho_ <- compute_rho(S_, I_)
  
  list(cor=rho_, S=S_, mu=mu_, I=I_, a=a_, b=b_)
}

# Gibbs sampler. target is function accepting parameter vector. initial is starting vector, to be accepted by target.
gibbs <- function(target, initial, iterations=10000, warmup=2000, thin=100) {
  k <- length(initial)
  chain <- matrix(0, iterations + warmup, k)
  
  chain[1, ] <- initial
  
  for (i in 2:(iterations + warmup)) {
    current <- chain[i - 1, ]
    candidate <- chain[i - 1, ]
    for (j in 1:k) {
      candidate[j] <- rnorm(1, mean=current[j], sd=0.5)
      A <- exp(target(candidate) - target(current))
      if (runif(1) < A) {
        current <- candidate
      } else {
        candidate <- current
      }
    }
    chain[i, ] <- current
  }
  return(chain[seq(warmup, iterations + warmup, thin), ])
}

# Log-likelihood function of negative binomial with gamma parameterisation designed for use with Gibbs sampler.
lk_nbinom_mcmc <- function(par, x, d, priors) {
  if (min(par) < 0) {
    return(-Inf)
  }
  
  r <- par[1]
  p <- par[2]/(par[2] + d)
  
  # Priors
  a1 <- priors$a1
  a2 <- priors$a2
  b1 <- priors$b1
  b2 <- priors$b2
  
  sum(dnbinom(x, r, p, log=TRUE)) + dgamma(par[1], a1, b1, log=TRUE) + dgamma(par[2], a2, b2, log=TRUE)
}

# Estimates the correlation between true and estimated network given count data x and sampling data d.
estimate_correlation_mcmc <- function(x, d) {
  parameters <- numerical_mle(x, d)
  
  priors <- list()
  
  priors$b1 <- 1/10 # Coefficient of dispersion = 10.
  priors$b2 <- 1/10
  priors$a1 <- parameters[1] * priors$b1
  priors$a2 <- parameters[2] * priors$b2
  
  target <- function(par) lk_nbinom_mcmc(par, x, d, priors)
  chain <- gibbs(target, c(parameters[1], parameters[2]), iterations=10000, warmup=1000, thin=1)
  
  a_ <- chain[, 1]
  b_ <- chain[, 2]
  
  S_ <- 1/sqrt(a_)
  mu_ <- a_/b_
  
  I_ <- mu_ * hmean(d)
  
  rho_ <- compute_rho(S_, I_)
    # correlation=quantile(rho_, probs=c(0.025, 0.5, 0.975))
  
  S_se <- sd(S_)
  rho_se <- sd(rho_)
  
  S_q <- quantile(S_, probs=c(0.025, 0.5, 0.975))
  rho_q <- quantile(rho_, probs=c(0.025, 0.5, 0.975))
  
  summary <- matrix(nrow=5, ncol=4)
  rownames(summary) <- c("Observed CV", "Interaction Rate", "Sampling Effort", "Social Differentiation", "Correlation")
  colnames(summary) <- c("Estimate", "SE", "Lower CI", "Upper CI")
  summary[1, 1] <- sd(x)/mean(x)
  summary[2, 1] <- mean(x/d)
  summary[3, 1] <- mean(x/d) * hmean(d)
  summary[4, ] <- c(S_q[2], S_se, S_q[1], S_q[3])
  summary[5, ] <- c(rho_q[2], rho_se, rho_q[1], rho_q[3])
  
  # list(param_summary=param_summary, cor=rho_, S=S_, mu=mu_, I=I_, a=a_, b=b_)
  summary <- signif(summary, 3)
  summary
}

# Estimates the correlation between true and estimated network given count data x and sampling data d.
estimate_correlation_interaction <- function(X, D, directed=FALSE) {
  if (directed) {
    x <- X[!diag(dim(X)[1])]
    d <- D[!diag(dim(D)[1])]
  } else {
    x <- X[upper.tri(X)]
    d <- D[upper.tri(D)]
  }
  return(estimate_correlation_mcmc(x, d))
}

# Uses the diminishing returns/elbow estimator to calculate optimal sampling effort.
estimate_dr <- function(S, rho_max=0.99, I_max=1000) {
  I <- seq(0, I_max, 0.01) # Choosing a maximum value here is difficult. It is the maximum mean number of interactions seen per dyad.
  
  rho_ <- compute_rho(S, I)
  
  rho_argmax <- which.min(abs(rho_ - rho_max))
  max_rho_ <- rho_[rho_argmax]
  max_I <- I[rho_argmax]

  theta <- atan2(max_rho_, max_I)
  co = cos(theta)
  si = sin(theta)
  rotation_matrix = matrix(c(co, -si, si, co), nrow=2, ncol=2)
  
  rotated_data <- cbind(rho_, I) %*% rotation_matrix
  
  rho_prime <- rho_[which.max(rotated_data[, 1])] # Optimal rho, by diminishing returns principle.
  
  I_prime <- I[which.max(rotated_data[, 1])] # Corresponding optimal sampling effort.
  
  list(I=I_prime, rho=rho_prime)
}