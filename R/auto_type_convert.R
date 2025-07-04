#' Automatically Convert Data Frame Column Types
#'
#' This function inspects the first \code{n_check} rows of each column in a data frame
#' to infer the most appropriate data type (e.g., numeric, integer, logical) and converts
#' the columns accordingly. It is particularly useful for datasets where columns are
#' initially read as character but should be of a different type.
#'
#' @param df A data frame whose columns require type conversion.
#' @param exclusion Character vector. Names of columns to exclude from type conversion.
#'   Default is an empty vector \code{c()}, meaning no columns are excluded.
#' @param n_check Integer. The number of initial rows to use for type inference.
#'   Default is 1000. If the data frame has fewer rows, all rows are used.
#' @param verbose Logical. If \code{TRUE}, prints progress and debugging messages
#'   using the \code{logger} package. Default is \code{FALSE}.
#'
#' @return A data frame with columns converted to their inferred data types.
#'
#' @details
#' The function uses the \code{type.convert} utility from base R to infer column types
#' based on the content of the first \code{n_check} rows. It then applies the inferred
#' types to the entire column, handling potential warnings or errors during conversion.
#' Supported types include numeric, integer, and logical. Columns already matching the
#' inferred type are left unchanged.
#'
#' @examples
#' \dontrun{
#' library(GWASkitR)
#'
#' # Example using a sample dataset
#'
#' # Check the structure of sumstats
#' sumstats <- GWASkitR::sumstats
#' sumstats$se <- as.character(sumstats$se)
#' sumstats$eaf <- as.character(sumstats$eaf)
#' str(sumstats)
#'
#' df <- GWASkitR:::auto_type_convert(df = sumstats, exclusion = c("se"), verbose = TRUE)
#'
#' # Check the structure of the converted data frame
#' str(df)
#' }
#'
#' @importFrom logger log_info log_debug log_warn log_error
auto_type_convert <- function(df, exclusion = c(), n_check = 1000, verbose = TRUE) {

  if (verbose) logger::log_info("Starting automatic type conversion for data frame with {ncol(df)} columns and {nrow(df)} rows.")

  # Select the first n_check rows (or all rows if fewer) for type inference
  n_rows_to_check <- min(nrow(df), n_check)
  df_head <- df[seq_len(n_rows_to_check), , drop = FALSE]
  if (verbose) logger::log_debug("Using the first {n_rows_to_check} rows for type inference.")

  # Use type.convert to infer the data type for each column
  # type.convert tries to convert character vectors to logical, integer, numeric, etc.
  # as.is = TRUE prevents automatic conversion of character vectors to factors
  col_types <- sapply(df_head, function(x) class(type.convert(x, as.is = TRUE)))

  if (verbose) logger::log_debug("Inferred column types: {paste(names(col_types), col_types, sep=':', collapse=', ')}")

  # Define a helper function to perform conversion to the target type
  auto_convert <- function(x, target_type, colname) {
    if (target_type == "numeric") {
      if (verbose) logger::log_debug("Converting column '{colname}' to numeric.")
      suppressWarnings(as.numeric(x))
    } else if (target_type == "integer") {
      logger::log_debug("Converting column '{colname}' to integer.")
      if (verbose) suppressWarnings(as.integer(x))
    } else if (target_type == "logical") {
      if (verbose) logger::log_debug("Converting column '{colname}' to logical.")
      suppressWarnings(as.logical(x))
    } else {
      if (verbose) logger::log_debug("No conversion needed for column '{colname}', keeping as {target_type}.")
      x # Return original vector if no conversion is needed
    }
  }

  # Loop through each column and apply the inferred conversion
  for (col in names(df)) {
    if (col %in% exclusion) {
      if (verbose) logger::log_info("Column '{col}' is excluded from conversion: skipping due to inclusion in [{paste(exclusion, collapse = ', ')}].")
      next
    }
    old_class <- class(df[[col]])
    target_type <- col_types[col]
    # Only convert if the class is different from the inferred type
    if (!(target_type %in% old_class)) {
      if (verbose) logger::log_info("Column '{col}': converting from {old_class} to {target_type}.")
      tryCatch({
        df[[col]] <- auto_convert(df[[col]], target_type, col)
      }, warning = function(w) {
        if (verbose) logger::log_warn("Column '{col}': conversion produced warnings: {conditionMessage(w)}")
        df[[col]]
      }, error = function(e) {
        if (verbose) logger::log_error("Column '{col}': conversion failed: {conditionMessage(e)}")
        df[[col]]
      })
    } else {
      if (verbose) logger::log_debug("Column '{col}': already of type {target_type}, no conversion applied.")
    }
  }

  if (verbose) logger::log_info("Automatic type conversion complete.")
  return(df)
}
