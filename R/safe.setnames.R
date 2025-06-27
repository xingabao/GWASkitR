#' Safely Rename Columns in a Data Table
#'
#' This function safely renames a column in a data table using \code{data.table::setnames},
#' but only if the old and new column names are different. It provides optional verbose
#' logging to track renaming operations.
#'
#' @param dt Data table, the input data table whose column needs to be renamed.
#' @param old Character, the current name of the column to be renamed.
#' @param new Character, the new name for the column.
#' @param verbose Logical, whether to log information about the renaming operation (default: FALSE).
#'
#' @return Invisible, the modified data table with the renamed column (if applicable).
#'
#' @examples
#' \dontrun{
#'   library(data.table)
#'   dt <- data.table(a = 1:3)
#'   safe.setnames(dt, old = "a", new = "b", verbose = TRUE)
#'   print(dt)
#' }
#'
#' @importFrom logger log_info
#' @importFrom data.table setnames
#'
safe.setnames <- function(dt, old, new, verbose = FALSE) {
  # Check if old and new names are different before renaming
  if (!identical(old, new)) {
    # Log the renaming operation if verbose is TRUE
    if (verbose) {
      logger::log_info("Renaming column '{old}' to '{new}'.")
    }
    # Perform the renaming using data.table::setnames
    data.table::setnames(dt, old = old, new = new)
  } else {
    # Log a message if no renaming is needed (only if verbose is TRUE)
    if (verbose) {
      logger::log_info("No renaming needed; old and new column names are identical: '{old}'.")
    }
  }
  # Return the data table invisibly
  invisible(dt)
}
