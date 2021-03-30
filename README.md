# Accuracy and Power Analysis of Interaction Networks

This repository contains code to compute the correlation between an estimated social network and the true, unobserved social network for networks constructed from interaction rates. Simulations are also included to verify the method. This repository also contains code to compute power analyses on interaction networks.

## Example

A working example similar to this one can be found in `srkw_data.Rmd`. The network correlation method has been implemented in the aninet package (https://github.com/MNWeiss/aninet/). In this example we use southern resident killer whale contact data from the aninet package.

First, load in aninet and load the square interaction matrix `X`, where `X[i, j]` is the number of times an interaction was seen between node `i` and `j`, and the square sampling matrix `D`, where `D[i, j]` is the amount of time spent sampling the pair `i` and `j` (units don't matter, as long as it's consistent).

```{r}
library(aninet)
X <- srkw_contact
D <- srkw_sampling
```

To compute the social differentiation and accuracy of this collected network, use the function `social_differentiation`. It is important to use the argument `method="Negative-binomial"` for interaction rate data. If you are using association data, where the number of associations are known per dyad, and the number of sampling periods per dyad is known, you will need to use `method="Whitehead"` or `method="Beta-binomial"`. The `directed` argument will only be used for interaction rate data, because interactions can be directed, unlike associations. It's important to set this correctly, otherwise it will result in incorrect estimates.

```{r}
summary <- social_differentiation(X, D, method="Negative-binomial", directed=FALSE)
summary
```

```
                       Estimate      SE Lower CI Upper CI
Observed CV               2.510      NA       NA       NA
Mean interactions         3.600      NA       NA       NA
Sampling Effort           3.300      NA       NA       NA
Social Differentiation    2.530 0.13600    2.280    2.830
Correlation               0.977 0.00235    0.972    0.982
```

This shows us that the observed coefficient of variation/social differentiation was 2.51, and the estimated true coefficient of variation/social differentiation was 2.53 with a 95% CI of [2.28, 2.83]. This network has a correlation of 97.7% with the true network, with a 95% CI of [97.2, 98.2].

If we want to conduct a post hoc power analysis on this network. We extract the number of individuals `n`, the social differentiation `s`, and the mean number of interactions per dyad `mu`. If the effect size is `r` between a covariate and network metric for the true network, the `power_nodal_interaction` will give us the power of the analysis given our amount of sampling. Note that this will currently only work for interaction rate networks.

```{r}
n <- dim(X)[1] # Number of individuals
s <- summary[4, 1] # Social differentiation
mu <- summary[3, 1] # Mean number of interactions

r <- 0.5 # Effect size in the true network.

power_nodal_interaction(n, r, s, mu, D, num_iters = 200)
```

```
$n
[1] 22

$r
[1] 0.5

$social_differentiation
[1] 2.54

$interaction_rate
[1] 3.597403

$power
[1] 0.71
```

If we aren't running nodal regression analysis, the power of the analysis will need to be calculated in a different way. As an alternative to this, we might want to use the elbow method to work out the optimal trade-off between sampling effort and correlation. We can use the `correlation_elbow_interaction` to do this. It takes the argument `rho_max`, which describes the maximum effective correlation we're biologically interested in. We choose `rho_max = 0.99` here, but it may be desirable to change this depending on the context.

```{r}
correlation_elbow_interaction(s, 0.99)
```

```
$sampling_effort
[1] 0.66

$correlation
[1] 0.8998976
```

This shows us that the ideal trade-off for `rho_max = 0.99` with this level of social differentiation is a sampling effort of 0.66 to achieve a correlation of just below 90%. This shows that the estimated network is likely a good representation of the true network.