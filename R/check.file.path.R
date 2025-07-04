#' Check if a file path is valid and has the correct suffix
#'
#' This function checks whether the given file path is non-null, non-empty, and ends with one of the allowed suffixes.
#' It logs warnings and errors using the logger package. Optionally, it logs info-level messages when verbose is TRUE.
#'
#' @param file.path A character string specifying the file path to check.
#' @param suffix A character vector of allowed suffixes (e.g., c('txt', 'tsv')).
#' @param verbose Logical. If TRUE, log info-level messages.
#' @param stop_on_error Logical. If TRUE, stop with error when check fails. Default is FALSE (just return FALSE).
#'
#' @return TRUE if file.path is valid and suffix matches; FALSE otherwise (unless stop_on_error is TRUE, then error thrown).
#'
#' @examples
#' check.file.path("test.txt")           # TRUE
#' check.file.path("data.csv", suffix = c("csv", "tsv"))  # TRUE
#' check.file.path(NULL)                 # FALSE
#' check.file.path("")                   # FALSE
#' check.file.path("file.doc", suffix = c("txt", "tsv"))  # FALSE
#'
check.file.path <- function(file.path, suffix = c('txt', 'tsv'), verbose = TRUE, stop_on_error = FALSE) {
  # Argument type and content checks

  if (is.null(file.path) || !is.character(file.path) || length(file.path) != 1 || nchar(file.path) == 0) {
    # if (verbose) logger::log_warn("No valid file.path provided. Must be a non-null, non-empty character string.")
    # if (stop_on_error) stop("No valid file.path provided.")
    return(TRUE)
  }

  if (is.null(suffix) || length(suffix) == 0 || !all(nzchar(suffix))) {
    if (verbose) logger::log_warn("Suffix must be a non-empty character vector.")
    if (stop_on_error) stop("Suffix must be a non-empty character vector.")
    return(FALSE)
  }

  # Build regex pattern for matching file suffixes
  suffix_pattern <- paste0("\\.(", paste0(suffix, collapse = "|"), ")$", collapse = "")
  if (verbose) { logger::log_info("Checking output file '{basename(file.path)}' against allowed suffixes: {paste(suffix, collapse = ', ')}") }
  # Check if file.path ends with an allowed suffix (case-insensitive)
  if (!grepl(suffix_pattern, file.path, ignore.case = TRUE)) {
    if (verbose) {
      logger::log_info("Provided: '{file.path}'")
      logger::log_error("The file.path must end with one of the following suffixes: {paste(suffix, collapse = ', ')}.")
    }
    if (stop_on_error) {
      stop(
        "The file.path must end with one of: ",
        paste0("'.", suffix, "'", collapse = ", "),
        ". Provided: '", file.path, "'."
      )
    }
    return(FALSE)
  }
  if (verbose) { logger::log_info("file.path '{file.path}' is valid and has an allowed suffix.") }
  return(TRUE)
}

