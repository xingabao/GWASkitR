#' Complete GWAS Statistics for Beta, OR, and Z
#'
#' This function completes missing statistical measures (Beta, Odds Ratio (OR), and Z-score)
#' in GWAS summary statistics. It requires a data frame with at least the standard error (SE)
#' and one of the statistics (Beta, OR, or Z) to compute the remaining two statistics.
#'
#' @param sumstats A data frame containing GWAS summary statistics. Must include a column for
#'   standard error (SE) and at least one of Beta, OR, or Z-score.
#' @param col.se Character string. The name of the column containing standard error (SE). Default is "se".
#' @param col.beta Character string or NULL. The name of the column containing Beta values. Default is NULL.
#' @param col.z Character string or NULL. The name of the column containing Z-scores. Default is NULL.
#' @param col.or Character string or NULL. The name of the column containing Odds Ratios (OR). Default is NULL.
#' @param col.p Character string or NULL. The name of the column to store computed p-values (if applicable). Default is NULL.
#' @param col.lp Logical. If TRUE and p-values are computed, apply log10 transformation to p-values. Default is FALSE.
#'
#' @return A data frame with completed statistics for Beta, OR, Z-score, and optionally p-values.
#'
#' @examples
#' \dontrun{
#' library(GWASkitR)
#'
#' # Example usage:
#' impute_beta_or_z(
#'   sumstats = GWASkitR::sumstats,
#'   col.se = 'se',
#'   col.beta = 'beta',
#'   col.p = 'pval',
#'   col.lp = TRUE
#' )
#' }
#'
#' @export
impute_beta_or_z <- function(
    sumstats,
    col.se,
    col.beta = NULL,
    col.z = NULL,
    col.or = NULL,
    col.p = NULL,
    col.lp = FALSE,
    verbose = TRUE
) {

  # Load the logger package
  if (!requireNamespace("logger", quietly = TRUE)) {
    stop("Package 'logger' is required but not installed.")
  }

  # Initialize logger
  if (verbose) logger::log_threshold(logger::INFO)

  # Validate input data frame
  if (!is.data.frame(sumstats)) {
    if (verbose) logger::log_error("Input 'sumstats' must be a data frame.")
    stop("Input 'sumstats' must be a data frame.")
  }

  if (verbose) logger::log_info("Starting imputation of GWAS statistics.")

  # Check if SE column exists in the data frame
  if (!is.valid(col.se) || !col.se %in% names(sumstats)) {
    if (verbose) logger::log_error("Standard error column '{col.se}' not found in data frame.")
    stop("Standard error column not found in data frame.")
  }

  # Validate that exactly one of Beta, Z, or OR is provided
  valid.cols <- sum(sapply(list(col.beta, col.z, col.or), function(x) !is.null(x) && nchar(as.character(x)) > 0))
  if (valid.cols != 1) {
    if (verbose) logger::log_error("Exactly one of 'col.beta', 'col.z', or 'col.or' must be provided.")
    stop("Exactly one of 'col.beta', 'col.z', or 'col.or' must be provided.")
  }

  # Log which statistic is being used as the input
  if (is.valid(col.beta)) {
    if (verbose) logger::log_info("Using Beta values from column '{col.beta}' for imputation.")
    if (!col.beta %in% names(sumstats)) {
      if (verbose) logger::log_error("Beta column '{col.beta}' not found in data frame.")
      stop("Beta column not found in data frame.")
    }
    # Compute Z-score and OR from Beta
    sumstats$Z <- sumstats[[col.beta]] / sumstats[[col.se]]
    sumstats$OR <- exp(sumstats[[col.beta]])
  } else if (is.valid(col.z)) {
    if (verbose) logger::log_info("Using Z-scores from column '{col.z}' for imputation.")
    if (!col.z %in% names(sumstats)) {
      if (verbose) logger::log_error("Z-score column '{col.z}' not found in data frame.")
      stop("Z-score column not found in data frame.")
    }
    # Compute Beta and OR from Z-score
    sumstats$beta <- sumstats[[col.z]] * sumstats[[col.se]]
    sumstats$OR <- exp(sumstats$beta)
  } else if (is.valid(col.or)) {
    if (verbose) logger::log_info("Using Odds Ratios from column '{col.or}' for imputation.")
    if (!col.or %in% names(sumstats)) {
      if (verbose) logger::log_error("OR column '{col.or}' not found in data frame.")
      stop("OR column not found in data frame.")
    }
    # Compute Beta and Z-score from OR
    sumstats$beta <- log(sumstats[[col.or]])
    sumstats$Z <- sumstats$beta / sumstats[[col.se]]
  }

  # Optionally compute p-values if col.p is provided
  # if (is.valid(col.p)) {
  # if (verbose) logger::log_info("Computing p-values and storing in column '{col.p}'.")
  # sumstats[[col.p]] <- 2 * (1 - pnorm(abs(sumstats$Z)))
  # Apply log10 transformation if col.lp is TRUE
  if (col.lp) {
    if (verbose) logger::log_info("Applying log10 transformation to p-values.")
    sumstats$LP <- -log10(sumstats[[col.p]])
  }
  # }

  # Warn about potential NA values in computed columns
  if (any(is.na(sumstats$Z)) || any(is.na(sumstats$OR)) || any(is.na(sumstats$beta))) {
    if (verbose) logger::log_warn("Some computed values (Z, OR, or Beta) contain NA due to missing or invalid input data.")
  }

  if (verbose) logger::log_info("Completed imputation of GWAS statistics.")
  return(sumstats)
}
