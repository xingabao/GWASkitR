#' Impute Standard Errors for GWAS Summary Statistics
#'
#' This function imputes standard errors (SE) for GWAS summary statistics using
#' either beta coefficients or odds ratios (OR) along with p-values. It calculates
#' the SE based on the provided statistics and stores the result in a specified column.
#'
#' @param sumstats A data frame containing GWAS summary statistics.
#' @param col.p Character. The name of the column in \code{sumstats} containing p-values.
#' @param col.beta Character, optional. The name of the column in \code{sumstats} containing beta coefficients.
#'   Provide this if SE should be imputed from beta values.
#' @param col.or Character, optional. The name of the column in \code{sumstats} containing odds ratios.
#'   Provide this if SE should be imputed from odds ratios.
#' @param col.se.set Character. The name of the column in \code{sumstats} where the imputed standard errors
#'   will be stored. Defaults to 'se'.
#' @param verbose Logical. If TRUE, log messages will be printed to track the progress and
#'   potential issues during imputation. Defaults to TRUE.
#'
#' @return A data frame with the original summary statistics and an additional column
#'   for the imputed standard errors.
#'
#' @details
#' Exactly one of \code{col.beta} or \code{col.or} must be provided to impute standard errors.
#' The function computes Z-scores from p-values and uses them to derive SE based on either
#' beta coefficients or odds ratios. If p-values or other input data are invalid or missing,
#' NA values may appear in the output.
#'
#' @examples
#' # Example with beta values
#' sumstats <- data.frame(pval = c(1e-5, 0.01), beta = c(0.5, -0.3))
#' result <- impute_se(sumstats, col.p = "pval", col.beta = "beta")
#'
#' # Example with odds ratios
#' sumstats_or <- data.frame(pval = c(1e-5, 0.01), or = c(1.2, 0.8))
#' result_or <- impute_se(sumstats_or, col.p = "pval", col.or = "or")
#'
#' # Example with built-in dataset
#' impute_se(sumstats = GWASkitR::sumstats, col.p = 'pval', col.beta = 'beta', col.se.set = 'se2')
#'
#' @importFrom logger log_threshold log_info log_warn log_error
#' @importFrom glue glue
#' @export
impute_se <- function(
    sumstats,
    col.p,
    col.beta = NULL,
    col.or = NULL,
    col.se.set = 'se',
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

  # Validate that exactly one of Beta or OR is provided
  valid.cols <- sum(sapply(list(col.beta, col.or), function(x) !is.null(x) && nchar(as.character(x)) > 0))
  if (valid.cols != 1) {
    if (verbose) logger::log_error("Exactly one of 'col.beta' or 'col.or' must be provided.")
    stop("Exactly one of 'col.beta' or 'col.or' must be provided.")
  }

  # Check if SE column exists in the data frame
  if (!is.valid(col.p) || !col.p %in% names(sumstats)) {
    if (verbose) logger::log_error("Standard error column '{col.p}' not found in data frame.")
    stop("Standard error column not found in data frame.")
  }

  # Check if SE column exists in the data frame
  if (!is.valid(col.se.set)) {
    if (verbose) logger::log_error("The name '{col.se.set}' is invalid. Please rename it.")
    stop(glue::glue("The name '{col.se.set}' is invalid. Please rename it."))
  }

  # Log input data dimensions for debugging
  if (verbose) logger::log_info(glue::glue("Input data frame has {nrow(sumstats)} rows and {ncol(sumstats)} columns."))

  # Log which statistic is being used as the input
  if (is.valid(col.beta)) {
    if (verbose) logger::log_info("Using Beta values from column '{col.beta}' for imputation.")
    if (!col.beta %in% names(sumstats)) {
      if (verbose) logger::log_error("Beta column '{col.beta}' not found in data frame.")
      stop("Beta column not found in data frame.")
    }
    # Compute Z-score and OR from Beta
    sumstats$Z <- qnorm(sumstats[[col.p]]/2, lower.tail = FALSE) * sign(sumstats[[col.beta]])
    sumstats[[col.se.set]] <- abs(sumstats[[col.beta]] / sumstats$Z)
  } else if (is.valid(col.or)) {
    if (verbose) logger::log_info("Using Odds Ratios from column '{col.or}' for imputation.")
    if (!col.or %in% names(sumstats)) {
      if (verbose) logger::log_error("OR column '{col.or}' not found in data frame.")
      stop("OR column not found in data frame.")
    }
    # Compute Beta and Z-score from OR
    sumstats$beta <- log(sumstats[[col.or]])
    sumstats$Z <- qnorm(sumstats[[col.p]]/2, lower.tail = FALSE) * sign(sumstats[[col.beta]])
    sumstats[[col.se.set]] <- abs(sumstats[[col.beta]] / sumstats$Z)
  }

  # Warn about potential NA values in computed columns
  if (any(is.na(sumstats[[col.se.set]]))) {
    if (verbose) logger::log_warn("Some computed values ({col.se.set}) contain NA due to missing or invalid input data.")
  }

  if (verbose) logger::log_info("Completed imputation of GWAS statistics.")
  return(sumstats)
}


