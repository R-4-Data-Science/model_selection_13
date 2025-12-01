#' @title Confusion matrix and classification metrics for logistic models
#'
#' @description Computes a confusion matrix and standard classification metrics at a given
#' probability cutoff (default 0.5) for binary outcomes. This is mainly intended
#' for evaluating logistic regression models in the vignette.
#'
#' The function expects observed binary outcomes \code{y} and predicted
#' probabilities \code{prob} for the positive class.
#'
#' @param y Observed binary response. Can be:
#'   \itemize{
#'     \item numeric (0/1),
#'     \item logical (\code{FALSE}/\code{TRUE}),
#'     \item factor with 2 levels (the second level is treated as "positive").
#'   }
#' @param prob Numeric vector of predicted probabilities for the positive class,
#'   of the same length as \code{y}.
#' @param cutoff Numeric cutoff in (0, 1) used to classify predictions as
#'   positive; default is 0.5.
#'
#' @return A list with components:
#'   \item{confusion}{2x2 matrix with rows = observed (0/1), columns = predicted (0/1),
#'         in the order: 0 = negative, 1 = positive.}
#'   \item{metrics}{Named numeric vector with:
#'         \code{prevalence}, \code{accuracy}, \code{sensitivity},
#'         \code{specificity}, \code{FDR}, \code{DOR}.}
#'
#' @details
#' The metrics are defined as:
#' \itemize{
#'   \item \code{prevalence} = (TP + FN) / N
#'   \item \code{accuracy}   = (TP + TN) / N
#'   \item \code{sensitivity} (recall, TPR) = TP / (TP + FN)
#'   \item \code{specificity} (TNR) = TN / (TN + FP)
#'   \item \code{FDR} (false discovery rate) = FP / (TP + FP)
#'   \item \code{DOR} (diagnostic odds ratio) =
#'         (TP / FN) / (FP / TN) = (TP * TN) / (FP * FN)
#' }
#' When a denominator is zero (e.g. no positives, no predicted positives, or a
#' cell count is zero), the corresponding metric is returned as \code{NA}.
#'
#' @export
confusion_metrics <- function(y, prob, cutoff = 0.5) {

  ## ---- Input checks ----

  if (length(y) != length(prob)) {
    stop("y and prob must have the same length.")
  }

  if (!is.numeric(prob)) {
    stop("prob must be a numeric vector of predicted probabilities.")
  }

  if (any(is.na(y)) || any(is.na(prob))) {
    stop("Missing values in y or prob are not allowed.")
  }

  if (cutoff <= 0 || cutoff >= 1) {
    stop("cutoff must be strictly between 0 and 1.")
  }

  # Convert y to numeric 0/1
  y_bin <- NULL

  if (is.factor(y)) {
    if (nlevels(y) != 2L) {
      stop("If y is a factor, it must have exactly 2 levels.")
    }
    # Treat the second level as positive (1)
    y_bin <- as.integer(y == levels(y)[2L])
  } else if (is.logical(y)) {
    y_bin <- as.integer(y)
  } else if (is.numeric(y)) {
    # We tolerate numeric 0/1; otherwise stop
    uniq <- sort(unique(y))
    if (!all(uniq %in% c(0, 1))) {
      stop("Numeric y must only contain 0 and 1.")
    }
    y_bin <- as.integer(y)
  } else {
    stop("y must be numeric (0/1), logical, or a 2-level factor.")
  }

  ## ---- Construct predicted classes ----

  # Predicted positive if prob >= cutoff
  y_hat <- as.integer(prob >= cutoff)

  ## ---- Compute confusion matrix components ----

  # Confusion counts:
  # TN = true negatives (y = 0, y_hat = 0)
  # FP = false positives (y = 0, y_hat = 1)
  # FN = false negatives (y = 1, y_hat = 0)
  # TP = true positives (y = 1, y_hat = 1)
  TN <- sum(y_bin == 0 & y_hat == 0)
  FP <- sum(y_bin == 0 & y_hat == 1)
  FN <- sum(y_bin == 1 & y_hat == 0)
  TP <- sum(y_bin == 1 & y_hat == 1)

  N <- TN + FP + FN + TP

  # Build a 2x2 confusion matrix with rows = observed, cols = predicted
  confusion <- matrix(
    c(TN, FP, FN, TP),
    nrow = 2, ncol = 2, byrow = TRUE,
    dimnames = list(
      observed  = c("0", "1"),
      predicted = c("0", "1")
    )
  )

  ## ---- Compute metrics ----

  # Prevalence: proportion of positives in the sample
  prevalence <- if (N > 0) (TP + FN) / N else NA_real_

  # Accuracy: proportion of correct predictions
  accuracy <- if (N > 0) (TP + TN) / N else NA_real_

  # Sensitivity (recall, TPR): TP / (TP + FN)
  denom_sens <- TP + FN
  sensitivity <- if (denom_sens > 0) TP / denom_sens else NA_real_

  # Specificity (TNR): TN / (TN + FP)
  denom_spec <- TN + FP
  specificity <- if (denom_spec > 0) TN / denom_spec else NA_real_

  # False discovery rate (FDR): FP / (TP + FP)
  denom_fdr <- TP + FP
  FDR <- if (denom_fdr > 0) FP / denom_fdr else NA_real_

  # Diagnostic odds ratio (DOR): (TP/FN) / (FP/TN) = (TP * TN) / (FP * FN)
  # Only defined when FP > 0 and FN > 0 and TN > 0 and TP > 0.
  if (TN > 0 && TP > 0 && FP > 0 && FN > 0) {
    DOR <- (TP * TN) / (FP * FN)
  } else {
    DOR <- NA_real_
  }

  metrics <- c(
    prevalence  = prevalence,
    accuracy    = accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    FDR         = FDR,
    DOR         = DOR
  )

  ## ---- Return ----

  list(
    confusion = confusion,
    metrics   = metrics
  )
}
