#' @title Overlap heatmap of plausible models
#'
#' @description Plot a heatmap of pairwise overlaps between plausible models, measured by
#' Jaccard similarity between sets of predictors.
#'
#' @param plaus A data.frame returned by \code{plausible_models()}, containing
#'   at least a list-column \code{vars} with the predictors in each model.
#' @param reorder Logical; if \code{TRUE} (default), rows/cols are reordered
#'   by hierarchical clustering on (1 - similarity). If \code{FALSE}, models
#'   are kept in their current order.
#' @param main Title for the plot.
#' @param ... Further arguments passed to \code{image()}.
#'
#' @return Invisibly returns the Jaccard similarity matrix.
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

  ## ---- Compute Jaccard similarity matrix ----
  jaccard <- function(a, b) {
    a <- unique(a); b <- unique(b)
    inter <- length(intersect(a, b))
    union <- length(union(a, b))
    if (union == 0L) return(0)
    inter / union
  }

  jac <- matrix(0, nrow = M, ncol = M)
  for (i in seq_len(M)) {
    for (j in seq_len(M)) {
      jac[i, j] <- jaccard(models_vars[[i]], models_vars[[j]])
    }
  }

  ## ---- Optional reordering by clustering ----
  ord <- seq_len(M)
  if (reorder) {
    dist_mat <- stats::dist(1 - jac)
    hc <- stats::hclust(dist_mat)
    ord <- hc$order
  }

  jac_ord <- jac[ord, ord]
  labels_ord <- model_labels[ord]

  ## ---- Draw heatmap with image() ----
  # Create a simple color palette from light to dark
  pal <- grDevices::colorRampPalette(c("white", "orange", "red"))(50)

  # image() draws matrix with rows bottom->top, so we flip it vertically
  image(
    x = seq_len(M),
    y = seq_len(M),
    z = t(jac_ord[nrow(jac_ord):1, ]),
    col = pal,
    axes = FALSE,
    xlab = "Model",
    ylab = "Model",
    main = main,
    ...
  )

  # Axes: x-axis left-to-right, y-axis top-to-bottom matching matrix
  axis(1, at = seq_len(M), labels = labels_ord, las = 2, cex.axis = 0.7)
  axis(2, at = seq_len(M), labels = rev(labels_ord), las = 2, cex.axis = 0.7)

  box()

  invisible(jac)
}
