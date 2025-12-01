#' @title Multi-path forward selection using AIC
#'
#' @description Builds a small "forest" of alternative forward-selection paths by
#' allowing near-ties in AIC at each step.
#'
#' @param x Matrix or data.frame of predictors (n x p).
#' @param y Response vector of length n.
#' @param family Model family: "gaussian" (lm) or "binomial" (glm).
#' @param K Maximum number of forward steps (max model size).
#' @param eps Minimum AIC improvement required to expand from a parent.
#' @param delta AIC tolerance for keeping near-tie child models (per parent).
#' @param L Maximum number of models to keep per level (after deduplication).
#'
#' @return A list of class "path_forest" with components:
#'   - frontiers: list of data.frames, one per step, describing models at that step.
#'   - aic_by_model: data.frame summarizing all models across all steps.
#'   - meta: list with input parameters and variable names.
#'
#' @export
build_paths <- function(x, y,
                        family = c("gaussian", "binomial"),
                        K = min(ncol(x), 10),
                        eps = 1e-6,
                        delta = 1,
                        L = 50) {

  ## ---- Input checks and preparation ----

  family <- match.arg(family)

  # Coerce x to data.frame and ensure column names
  if (is.matrix(x)) {
    x <- as.data.frame(x)
  } else {
    x <- as.data.frame(x)
  }

  if (is.null(colnames(x))) {
    colnames(x) <- paste0("x", seq_len(ncol(x)))
  }

  n <- NROW(x)
  p <- NCOL(x)

  if (length(y) != n) {
    stop("Length of y must match number of rows in x.")
  }

  # Build working data.frame with y + predictors
  df <- data.frame(y = y, x, check.names = FALSE)
  var_names <- colnames(x)

  # Helper: fit model and compute AIC
  fit_model_aic <- function(vars) {
    # vars: character vector of predictor names
    if (length(vars) == 0L) {
      # Intercept-only model
      form <- y ~ 1
    } else {
      form <- as.formula(
        paste("y ~", paste(vars, collapse = " + "))
      )
    }

    if (family == "gaussian") {
      fit <- stats::lm(form, data = df)
    } else {
      fit <- stats::glm(form, family = stats::binomial(), data = df)
    }
    aic_val <- stats::AIC(fit)

    list(fit = fit, aic = aic_val, formula = form)
  }

  # Helper: build a unique key for a model (sorted variable names)
  model_key <- function(vars) {
    # Always sort vars for a canonical representation
    if (length(vars) == 0L) {
      return("<empty>")
    }
    paste(sort(vars), collapse = " + ")
  }

  ## ---- Initialization: level 0 (empty model) ----

  # Fit the intercept-only model once
  empty_fit <- fit_model_aic(character(0))

  # We store current_models as a list of model descriptors
  # Each element is a list with:
  #   vars   : character vector of included predictors
  #   aic    : numeric AIC value
  #   key    : unique model key
  #   parent : key of parent model (NA for step 0)
  #   step   : step index (0 for empty)
  current_models <- list(
    list(
      vars   = character(0),
      aic    = empty_fit$aic,
      key    = model_key(character(0)),
      parent = NA_character_,
      step   = 0L,
      formula = empty_fit$formula
    )
  )

  # Frontiers: store level 0 explicitly if desired
  frontiers <- list()
  frontiers[[1]] <- data.frame(
    step   = 0L,
    model_id = seq_along(current_models),
    key    = vapply(current_models, function(m) m$key, character(1)),
    aic    = vapply(current_models, function(m) m$aic, numeric(1)),
    size   = vapply(current_models, function(m) length(m$vars), integer(1)),
    parent = vapply(current_models, function(m) m$parent, character(1)),
    vars   = I(lapply(current_models, function(m) m$vars)),
    stringsAsFactors = FALSE
  )

  # We also maintain a global table of all models encountered
  all_models <- frontiers[[1]]

  ## ---- Forward steps ----

  # Outer loop over up to K steps
  for (k in seq_len(K)) {

    # List to store all candidate children from all parents at this step
    candidate_children <- list()
    child_idx <- 1L

    # For each parent model at current frontier
    for (p_idx in seq_along(current_models)) {
      parent <- current_models[[p_idx]]
      parent_vars <- parent$vars
      parent_aic  <- parent$aic
      parent_key_ <- parent$key

      # Determine remaining variables that are not yet in parent
      remaining <- setdiff(var_names, parent_vars)
      if (length(remaining) == 0L) {
        # Cannot expand this parent further
        next
      }

      # For each remaining variable, build a candidate child model
      for (v in remaining) {
        child_vars <- sort(c(parent_vars, v))
        key <- model_key(child_vars)

        # Fit and compute AIC
        res <- fit_model_aic(child_vars)

        candidate_children[[child_idx]] <- list(
          step      = k,
          parent    = parent_key_,
          parent_aic = parent_aic,
          vars      = child_vars,
          key       = key,
          added     = v,
          aic       = res$aic,
          formula   = res$formula
        )
        child_idx <- child_idx + 1L
      }
    }

    # If no candidates at all, we stop
    if (length(candidate_children) == 0L) {
      break
    }

    # Turn candidate list into a data.frame for easier handling
    cand_df <- data.frame(
      step      = vapply(candidate_children, function(z) z$step, integer(1)),
      parent    = vapply(candidate_children, function(z) z$parent, character(1)),
      parent_aic= vapply(candidate_children, function(z) z$parent_aic, numeric(1)),
      key       = vapply(candidate_children, function(z) z$key, character(1)),
      added     = vapply(candidate_children, function(z) z$added, character(1)),
      aic       = vapply(candidate_children, function(z) z$aic,  numeric(1)),
      stringsAsFactors = FALSE
    )
    # Keep the vars/formula in parallel list (list columns are fine, but
    # here we keep them in a separate list to simplify handling).
    cand_vars    <- lapply(candidate_children, function(z) z$vars)
    cand_formula <- lapply(candidate_children, function(z) z$formula)

    ## ---- Apply delta/eps rules per parent ----

    # For each parent, we:
    #  1) find best child AIC,
    #  2) check improvement vs parent_aic,
    #  3) keep children within delta of that best AIC if improvement >= eps.
    keep_idx <- logical(nrow(cand_df))
    parent_keys_unique <- unique(cand_df$parent)

    for (pk in parent_keys_unique) {
      idx <- which(cand_df$parent == pk)
      parent_aic_val <- cand_df$parent_aic[idx[1]]  # all same

      # Best AIC among children of this parent
      best_child_aic <- min(cand_df$aic[idx])

      # Check if best child improves AIC by at least eps
      if ((parent_aic_val - best_child_aic) >= eps) {
        # Keep all near-ties within delta
        near_idx <- idx[ cand_df$aic[idx] <= (best_child_aic + delta) ]
        keep_idx[near_idx] <- TRUE
      } else {
        # No children kept for this parent
        next
      }
    }

    # Filter candidates
    if (!any(keep_idx)) {
      # No expanding parent had sufficient improvement -> stop algorithm
      break
    }

    cand_df_kept   <- cand_df[keep_idx, , drop = FALSE]
    cand_vars_kept <- cand_vars[keep_idx]
    cand_form_kept <- cand_formula[keep_idx]

    ## ---- Deduplicate across parents: same model key ----

    # Multiple parents may produce the same child model (same var set).
    # We deduplicate by key, keeping the one with smallest AIC.
    # We'll build an index of best per key.
    unique_keys <- unique(cand_df_kept$key)
    dedup_idx   <- integer(length(unique_keys))

    for (i in seq_along(unique_keys)) {
      k_ <- unique_keys[i]
      idx <- which(cand_df_kept$key == k_)
      # Among these, pick the smallest AIC
      best_idx <- idx[ which.min(cand_df_kept$aic[idx]) ]
      dedup_idx[i] <- best_idx
    }

    dedup_idx <- sort(dedup_idx)
    dedup_df   <- cand_df_kept[dedup_idx, , drop = FALSE]
    dedup_vars <- cand_vars_kept[dedup_idx]
    dedup_form <- cand_form_kept[dedup_idx]

    ## ---- Limit to at most L models by best AIC ----

    if (nrow(dedup_df) > L) {
      order_idx <- order(dedup_df$aic)
      keep_L    <- order_idx[seq_len(L)]
      dedup_df   <- dedup_df[keep_L, , drop = FALSE]
      dedup_vars <- dedup_vars[keep_L]
      dedup_form <- dedup_form[keep_L]
    }

    # Build current_models for next iteration
    current_models <- lapply(seq_len(nrow(dedup_df)), function(i) {
      list(
        vars    = dedup_vars[[i]],
        aic     = dedup_df$aic[i],
        key     = dedup_df$key[i],
        parent  = dedup_df$parent[i],
        step    = k,
        formula = dedup_form[[i]]
      )
    })

    # Build frontier data.frame for this step
    frontier_k <- data.frame(
      step     = k,
      model_id = seq_along(current_models),
      key      = vapply(current_models, function(m) m$key, character(1)),
      aic      = vapply(current_models, function(m) m$aic, numeric(1)),
      size     = vapply(current_models, function(m) length(m$vars), integer(1)),
      parent   = vapply(current_models, function(m) m$parent, character(1)),
      vars     = I(lapply(current_models, function(m) m$vars)),
      stringsAsFactors = FALSE
    )

    frontiers[[k + 1L]] <- frontier_k
    all_models <- rbind(all_models, frontier_k)
  }

  ## ---- Assemble return object ----

  meta <- list(
    family   = family,
    K        = K,
    eps      = eps,
    delta    = delta,
    L        = L,
    n        = n,
    p        = p,
    var_names = var_names
  )

  out <- list(
    frontiers    = frontiers,
    aic_by_model = all_models,
    meta         = meta
  )
  class(out) <- "path_forest"

  out
}
