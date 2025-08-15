#' Check if a list of column names is valid
#'
#' This function checks whether the input column names (passed as separate arguments) are non-empty,
#' do not contain NA values, and do not contain empty strings.
#'
#' @param ... One or more column names as character values.
#'
#' @return Logical. Returns `TRUE` if all provided column names are not NA or empty strings;
#' returns `FALSE` otherwise, or if no arguments are supplied.
#'
#' @examples
#' is.valid("a", "b")           # TRUE
#' is.valid("a", "")            # FALSE
#' is.valid()                   # FALSE
#' is.valid(NA, "b")            # FALSE
#' is.valid("a", NULL)          # FALSE
#'
is.valid <- function(...) {
  # Collect all arguments into a list to preserve NULL values
  args <- list(...)

  # Return FALSE if no arguments are supplied
  if (length(args) == 0) return(FALSE)

  # Check if any argument is NULL
  if (any(sapply(args, is.null))) return(FALSE)

  # Convert arguments to character vector for further checks
  cols <- unlist(args)

  # Check each element: not NA and not an empty string
  all(!is.na(cols) & cols != "")
}

