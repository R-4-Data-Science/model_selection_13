## ----setup, message = FALSE---------------------------------------------------
library(modelselection13)
library(care)

set.seed(28)

## ----data-prep----------------------------------------------------------------
# Load the diabetes progression data from the 'care' package
data("efron2004", package = "care")

# Convert to a plain numeric matrix (remove AsIs)
X_base_mat <- as.matrix(unclass(efron2004$x))
dim(X_base_mat) 

# Convert to data.frame
X_base <- as.data.frame(X_base_mat)
dim(X_base) 

y <- as.numeric(efron2004$y[, 1])

# ---- Build second-order terms: quadratics + pairwise interactions ----

# Quadratic terms: x_j^2 for each baseline predictor
X_quad <- X_base^2
colnames(X_quad) <- paste0(colnames(X_base), "_sq")

# Pairwise interactions: x_j * x_k for j < k
p <- ncol(X_base)
interaction_list <- vector("list", length = p * (p - 1) / 2)
interaction_names <- character(length(interaction_list))

idx <- 1L
for (j in seq_len(p - 1L)) {
  for (k in (j + 1L):p) {
    interaction_list[[idx]] <- X_base[[j]] * X_base[[k]]
    interaction_names[idx] <- paste0(colnames(X_base)[j], ":", colnames(X_base)[k])
    idx <- idx + 1L
  }
}
X_int <- as.data.frame(interaction_list)
colnames(X_int) <- interaction_names

# Combine all predictors: main effects + quadratics + interactions
X_full <- cbind(X_base, X_quad, X_int)

dim(X_full)

# ---- Train/test split (e.g., 70/30) ----

n <- nrow(X_full)
train_prop <- 0.7
n_train <- floor(train_prop * n)

set.seed(2025)  # reproducible split
train_idx <- sample(seq_len(n), size = n_train, replace = FALSE)

X_train <- X_full[train_idx, , drop = FALSE]
y_train <- y[train_idx]

X_test  <- X_full[-train_idx, , drop = FALSE]
y_test  <- y[-train_idx]

# Quick sanity checks
dim(X_train)
dim(X_test)
length(y_train)
length(y_test)

## ----build-paths, message=FALSE-----------------------------------------------
# Multi-path AIC forward selection on the training set
set.seed(28)

K_val <- 10     # max number of variables to include
eps_val <- 1e-6 # minimum AIC improvement
delta_val <- 2  # allow small AIC near-ties
L_val <- 30     # limit models per step

paths <- build_paths(
  x = X_train,
  y = y_train,
  family = "gaussian",
  K = K_val,
  eps = eps_val,
  delta = delta_val,
  L = L_val
)

paths$meta
length(paths$frontiers)          # number of steps+1
paths$frontiers[[1]]             # step 0: empty model
paths$frontiers[[2]][1:5, ]      # first few models at step 1

## ----stability, message=FALSE-------------------------------------------------
set.seed(28)

stab <- stability(
  x = X_train,
  y = y_train,
  B = 50,
  resample = "bootstrap",
  family = "gaussian",
  K = K_val,     
  eps = eps_val,
  delta = delta_val,
  L = L_val
)

# Show first few stability scores
head(stab$pi)
summary(stab$pi)

## ----plausible-models---------------------------------------------------------
Delta_val <- 2        # AIC window
tau_val <- 0.4       # stability threshold

plaus <- plausible_models(
  forest = paths,
  pi = stab$pi,
  Delta = Delta_val,
  tau = tau_val,
  jaccard_threshold = 0.8
)

plaus

