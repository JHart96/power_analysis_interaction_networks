source("network_correlation.R")
source("power_analysis.R")

# Simulate data with 20 nodes, 30% network density, and 25 mean units of observation time per dyad.
sim_data <- simulate_data(20, 0.3, 20)

X <- sim_data$X # 20 x 20 symmetric matrix of integer observation counts.
D <- sim_data$D # 20 x 20 matrix of positive real-valued sampling times.
A <- sim_data$A # 20 x 20 matrix of true interaction rates (usually unknown in empirical studies).

# Calculate network correlation.
summary_obj <- estimate_correlation_interaction(X, D)
summary_obj

social_differentiation <- summary_obj[4, 1]
interaction_rate <- summary_obj[2, 1]
sampling_times <- D # Matrix of sampling times OR Single value of mean sampling times

# Calculate power of nodal regression for effect size r = 0.5
power_analysis_nodal(20, 0.5, social_differentiation, interaction_rate, sampling_times)

# Check our sampling is sufficient using the diminishing returns/elbow method:
power_analysis_elbow(social_differentiation, 0.99)
