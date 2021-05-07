# Accuracy and Power Analysis of Interaction Networks

This repository contains code to compute the correlation between an estimated social network and the true, unobserved social network for networks constructed from interaction rates. Simulations are also included to verify the method. This repository also contains code to compute power analyses on interaction networks.

## Example

Load in R files using `source`:

```
source("network_correlation.R")
source("power_analysis.R")
```

Simulate data with 20 nodes, 50% network density, and 20 mean units of observation time per dyad. Of course this won't be necessary if you're using your own data.
```{r}
sim_data <- simulate_data(20, 0.5, 20)

X <- sim_data$X # 20 x 20 symmetric matrix of integer observation counts.
D <- sim_data$D # 20 x 20 matrix of positive real-valued sampling times.
```

Use the two matrices to estimate the correlation between the sampled network and the true, underlying network. This may take a few seconds as it uses an MCMC algorithm to estimate the confidence interval of the correlation.
```{r}
summary_obj <- estimate_correlation_interaction(X, D)
summary_obj
```
```
                       Estimate     SE Lower CI Upper CI
Observed CV               1.900     NA       NA       NA
Interaction Rate          0.159     NA       NA       NA
Sampling Effort           1.500     NA       NA       NA
Social Differentiation    2.430 0.2420    2.030     2.98
Correlation               0.948 0.0131    0.918     0.97
```

If you want to conduct power analysis for nodal regression, extract social differentiation, interaction rate, and sampling times from the summary object and the data matrices. To conduct the power analysis, provide a `r` value. This is the correlation coefficient and reflects the true relationship between the response and predictors in the regression (the effect size we'd see with perfect sampling). 
```{r}
social_differentiation <- summary_obj[4, 1]
interaction_rate <- summary_obj[2, 1]
sampling_times <- D # Matrix of sampling times OR Single value of mean sampling times

# Calculate power of nodal regression for effect size r = 0.5
power_analysis_nodal(20, 0.5, social_differentiation, interaction_rate, sampling_times)
```
```
$n
[1] 20

$r
[1] 0.5

$social_differentiation
[1] 2.43

$interaction_rate
[1] 0.159

$power
[1] 0.62
```

This shows us that we could expect a power of 62% given the properties of the data and the true effect size r = 0.5.

If a different type of analysis is being conducted, the diminishing returns/elbow estimator could be used to determine if sufficient data are available:
```{r}
power_analysis_elbow(social_differentiation, rho_max=0.99) # Use rho_max=0.99 as in the paper.
```
```
$sampling_effort
[1] 0.72

$correlation
[1] 0.8997662
```

The elbow method says we need a sampling effort of 0.72 to reach the optimal level of correlation, which is roughly 90%. From our run of `estimate_correlation_interaction` we know that sampling effort is 1.5, more than enough to get the optimal level of correlation. This indicates that the network could be subset to create two networks if required.
