power_analysis_nodal <- function(n, r, social_differentiation, interaction_rate, sampling_times, metric_fn=NULL, num_iters=200) {
  
  if (is.null(metric_fn)) {
    metric_fn <- igraph::strength
  }
  
  if (length(sampling_times) == 1) {
    sampling_times <- rpois(n, sampling_times)
    sampling_times[sampling_times == 0] = 1
  }
  
  a <- 1/social_differentiation^2
  b <- 1/(social_differentiation^2 * interaction_rate)
  
  pvals <- rep(0, 1)
  
  for (i in 1:num_iters) {
    alpha <- rgamma(n^2, a, b) # Actual rate.
    net <- igraph::graph_from_adjacency_matrix(matrix(alpha, n, n), weighted=TRUE)
    
    metric <- metric_fn(net)
  
    effect_size <- r * sqrt(1/(1 - r))/(sd(metric) * sqrt(r + 1))
    
    trait <- 1 + effect_size * metric + rnorm(n, sd=1)
    
    alpha_hat <- rpois(n^2, alpha*sampling_times)/sampling_times # Sampled rate.
    
    # Build network
    net_est <- igraph::graph_from_adjacency_matrix(matrix(alpha_hat, n, n), weighted=TRUE)
    metric_est <- metric_fn(net_est)
    
    fit_summary <- node_regression(metric_est, trait)
    pvals[i] <- fit_summary$pval
  }
  power <- mean(pvals < 0.05)
  list(n=n, r=r, social_differentiation=social_differentiation, interaction_rate=interaction_rate, power=power)
}

power_analysis_elbow <- function(S, rho_max, I_max=1000) {
  I <- seq(0, I_max, 0.01) # Maximum sampling effort.
  
  rho_ <- (S * sqrt(I))/(sqrt(1 + I * S^2))
  
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
  
  list(sampling_effort=I_prime, correlation=rho_prime)
}
