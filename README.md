# modelselection13

R package for multi-path AIC stepwise selection, resampling-based stability, and plausible model sets.

## Installation

``` r
# install.packages("remotes") 
remotes::install_github("R-4-Data-Science/Final_Project_13") 
library(modelselection13)
```

## Linear Regression

``` r
set.seed(1)
n <- 80
p <- 6
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
y <- 2*X[,1] - 1.5*X[,3] + rnorm(n)

# Build multi-path forest
forest <- build_paths(
  x = X, y = y,
  family = "gaussian",
  K = min(p, 10),
  eps = 1e-6,
  delta = 1,
  L = 50
)

# Compute stability
stab <- stability(
  x = X, y = y,
  B = 30,
  resample = "bootstrap",
  family = "gaussian",
  K = min(p, 10),
  eps = 1e-6,
  delta = 1,
  L = 50
)

# Plausible models
plaus <- plausible_models(
  forest = forest,
  pi = stab$pi,
  Delta = 2,
  tau = 0.6
)

plaus
```

## Logistic Regression

``` r
set.seed(2)
n <- 120
p <- 5
Xb <- matrix(rnorm(n * p), n, p)
colnames(Xb) <- paste0("x", 1:p)

eta <- 1.2*Xb[,1] - 1.0*Xb[,2]
prob <- 1/(1 + exp(-eta))
ybin <- rbinom(n, 1, prob)

# Build multi-path forest
forest_b <- build_paths(
  x = Xb, y = ybin,
  family = "binomial",
  K = min(p, 10),
  eps = 1e-6,
  delta = 1,
  L = 50
)

# Stability
stab_b <- stability(
  x = Xb, y = ybin,
  B = 30,
  resample = "bootstrap",
  family = "binomial",
  K = min(p, 10),
  eps = 1e-6,
  delta = 1,
  L = 50
)

# Plausible models
plaus_b <- plausible_models(
  forest = forest_b,
  pi = stab_b$pi,
  Delta = 2,
  tau = 0.6
)

plaus_b
```

## Optional extras: overlap heatmap and branching plot

### Overlap heatmap of plausible models

``` r
# Using 'plaus' from the linear example above
if (nrow(plaus) >= 2) {
  jac_mat <- plot_heatmap(
    plaus,
    reorder = TRUE    
  )
  jac_mat             
}
```

### Branching structure of the multi-path search

``` r
# Using 'forest' from the linear example above
plot_branching(
  forest,
  vertex_size = 18,
  label_cex   = 0.7
)
```
