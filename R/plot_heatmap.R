#' @title Overlap heatmap of plausible models
#'
#' @description Plot a heatmap of pairwise overlaps between plausible models,
#' typically obtained from \code{plausible_models()}. The overlap is measure
#' using the Jaccard similarity between sets of predictors:
#' \deqn{J(S_i, S_j) = |S_i \cap S_j| / |S_i \cup S_j|.}
#'
#' This visualization helps to see which plausible models are almost identical
#' (high overlap) and which provide genuinely different combinations of
#' predictors.
#'
#' @param plaus A data.frame returned by \code{plausible_models()}, containing
#'   at least a list-column \code{vars} with the predictors in each model.
#' @param reorder Logical; if \code{TRUE} (default), rows/columns are
#'   hierarchically clustered and reordered by similarity. If \code{FALSE},
#'   models are kept in their current order.
#' @param main Character string: title for the heatmap.
#' @param ... Additional arguments passed to \code{stats::heatmap()}.
#'
#' @return Invisibly returns the Jaccard similarity matrix used for plotting.
#'
#' @export
plot_heatmap <- function(plaus,
                         reorder = TRUE,
                         main = "Overlap Heatmap of Plausible Models",
                         ...) {

  ## ---- Basic checks ----

  if (is.null(plaus) || nrow(plaus) == 0L) {
    stop("plaus must be a non-empty data.frame of plausible models.")
  }

  if (!"vars" %in% names(plaus)) {
    stop("plaus must contain a list-column named 'vars' with predictor sets.")
  }

  models_vars <- plaus$vars
  M <- length(models_vars)

  if (M < 2L) {
    stop("Need at least two models to compute an overlap heatmap.")
  }

  # Use model_id or key as labels if available
  if ("model_id" %in% names(plaus)) {
    model_labels <- paste0("M", plaus$model_id)
  } else if ("key" %in% names(plaus)) {
    model_labels <- plaus$key
  } else {
    model_labels <- paste0("M", seq_len(M))
  }

  ## ---- Jaccard similarity matrix ----

  # Helper: Jaccard similarity between two character vectors
  jaccard <- function(a, b) {
    a <- unique(a); b <- unique(b)
    inter <- length(intersect(a, b))
    union <- length(union(a, b))
    if (union == 0L) return(0)
    inter / union
  }

  # Build M x M similarity matrix
  jac <- matrix(0, nrow = M, ncol = M)
  for (i in seq_len(M)) {
    for (j in seq_len(M)) {
      jac[i, j] <- jaccard(models_vars[[i]], models_vars[[j]])
    }
  }

  colnames(jac) <- model_labels
  rownames(jac) <- model_labels

  ## ---- Plot heatmap ----

  # If reorder = FALSE, suppress dendrograms and clustering
  if (!reorder) {
    stats::heatmap(jac,
                   Rowv = NA, Colv = NA,
                   scale = "none",
                   col = heat.colors(20),
                   main = main,
                   xlab = "Model", ylab = "Model",
                   ...)
  } else {
    # Default hierarchical clustering on (1 - similarity) as a distance
    dist_mat <- as.dist(1 - jac)
    stats::heatmap(jac,
                   Rowv = stats::hclust(dist_mat),
                   Colv = stats::hclust(dist_mat),
                   scale = "none",
                   col = heat.colors(20),
                   main = main,
                   xlab = "Model", ylab = "Model",
                   ...)
  }

  invisible(jac)
}
