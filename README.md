# Estimating Accuracy of Sampling with Interaction Rate Data - Derivation

If we assume the population association rates $\lambda_i$ are gamma-distributed, we have:
$$
\lambda_i \sim \text{Gamma}(\alpha, \beta).
$$
where $\alpha, \beta$ are the population-level parameters.

The observed number of interactions $X_i$ for the $i$-th dyad is likely to be Poisson-distributed, and will depend on the number of times the dyad were observed $d_i$. This gives us:
$$
X_i \sim \text{Poisson}(d_i \lambda_i) .
$$
Since $d_i$ is a constant, the distribution of $d_i \lambda_i$ is given by
$$
d_i \lambda_i \sim \text{Gamma}\left(\alpha, \frac{\beta}{d_i}\right).
$$
This kind of model is known as a Poisson-Gamma mixture, and can be described by the negative binomial in terms of the gamma parameters $\alpha, \beta$, and the sampling effort $d_i$:

$$
X_i \sim \text{NegBinomial}\left(\alpha, \frac{\beta}{\beta + d_i}\right).
$$

Typically when we estimate the association rate, we compute 
$$
\hat{\lambda_i} = \frac{X_i}{d_i}.
$$

Our aim is to compute the total variance in the estimated association rates $\hat{\lambda}_i$ explained by the total variance in the true association rates $\lambda_i$, that is:
$$
\text{EV} = \frac{\sum_{i=1}^n \text{Var}(\hat{\lambda}_i)}{\sum_{i=1}^n \text{Var}(\lambda_i)} = \frac{\hat{T}}{T}
$$
First, the variance of the true association rate is simply given by the variance of the gamma distribution:
$$
\text{Var}(\lambda_i) = \frac{\alpha}{\beta^2}.
$$

Summing over $n$ dyads we get
$$
T = \sum_{i=1}^n \text{Var}(\lambda_i) = \frac{n \alpha}{\beta^2}.
$$

The variance of the estimated association rates is slightly trickier, we use the variance of the negative binomial to get:
$$
\begin{aligned}
\text{Var}(\hat{\lambda_i})  &= \text{Var} \left( \frac{X_i}{d_i}\right) \\ &= \frac{1}{d_i^2} \text{Var}(X_i) \\ &= \frac{\alpha d_i (\beta + d_i)}{d_i^2\beta^2}.
\end{aligned}
$$

Again, summing over the $n$ dyads we get
$$
\begin{aligned}
\hat{T} &= \sum_{i=1}^n \text{Var}(\hat{\lambda}_i) \\ &=\sum_{i=1}^n \frac{\alpha d_i (\beta + d_i)}{d_i^2 \beta^2} \\ &= \frac{\alpha}{\beta^2} \sum_{i=1}^n \left(\frac{d_i^2 + d_i \beta}{d_i^2} \right) \\ &= \frac{\alpha}{\beta^2} \sum_{i=1}^n \left(1 + \frac{\beta}{d_i}\right) \\ &= \frac{\alpha}{\beta^2} \left(n + \beta \sum_{i=1}^n d_i^{-1}\right).
\end{aligned}
$$

Bringing these together, we get:
$$
\begin{aligned}
\text{EV} & = T\hat{T}^{-1} \\&= \frac{n\alpha}{\beta^2} \frac{\beta^2}{\alpha} \left(n + \beta \sum_{i=1}^n d_i^{-1}\right)^{-1} \\&= \frac{n}{n + \beta \sum_{i=1}^n d_i^{-1}}.
\end{aligned}
$$

The coefficient of variation/social differentiation of the interaction rate under this method is:
$$
S_\lambda = \text{CV}(\lambda) = \frac{\sigma(\lambda)}{\mu(\lambda)} = \frac{\alpha}{\beta^2} \frac{\beta}{\alpha} = \frac{1}{\beta}.
$$

So as $S_\lambda \rightarrow \infty$ we have $\text{EV} \rightarrow 1$, and as $S_\lambda \rightarrow 0$ we have $\text{EV} \rightarrow 0$. Therefore increases in social differentiation increase the explained variance. 