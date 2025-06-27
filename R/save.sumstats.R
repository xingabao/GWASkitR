#' Save summary statistics data to a text or TSV file
#'
#' This function saves a data frame or data table (such as summary statistics) to a `.txt` or `.tsv` file.
#' It checks the file extension, creates the target directory if necessary, writes the table in tab-delimited format,
#' and logs each important step using the `logger` package.
#'
#' @param sumstats.dat A data frame or data table containing summary statistics to be saved.
#' @param save.path File path where the summary statistics will be saved. Must end with `.txt` or `.tsv`.
#'
#' @return No return value. The function is called for its side effect of writing to disk and logging.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(a = 1:3, b = c("x", "y", "z"))
#' save.sumstats(df, "output/sumstats.txt")
#' }
#'
#' @importFrom data.table fwrite
#' @importFrom logger log_info log_warn log_error
#'
save.sumstats <- function(sumstats.dat, save.path) {
  # Check if the input data is provided
  if (is.null(sumstats.dat)) {
    logger::log_warn("Input data (sumstats.dat) is NULL. Nothing will be saved.")
    return(invisible(NULL))
  }

  # Check if save.path is provided
  if (is.null(save.path)) {
    logger::log_warn("No save.path provided. Cannot save summary statistics.")
    return(invisible(NULL))
  }

  # Check if save.path ends with .txt or .tsv (case-insensitive)
  if (!grepl("\\.(txt|tsv)$", save.path, ignore.case = TRUE)) {
    logger::log_error("The save.path must end with '.txt' or '.tsv'. Provided: '{save.path}'")
    stop("Error: The save.path must end with '.txt' or '.tsv'.")
  }

  # Extract directory from the file path
  save.dir <- dirname(save.path)

  # Create the directory if it doesn't already exist (including parent directories)
  if (!dir.exists(paths = save.dir)) {
    logger::log_info("Directory '{save.dir}' does not exist. Creating now...")
    dir.create(path = save.dir, recursive = TRUE)
    logger::log_info("Directory '{save.dir}' created successfully.")
  }

  # Try to write the data as a tab-delimited file
  tryCatch({
    logger::log_info("Saving formatted summary statistics to '{save.path}' ...")
    data.table::fwrite(
      as.data.frame(sumstats.dat),
      file = save.path,
      sep = "\t",
      row.names = FALSE
    )
    logger::log_info("Formatted summary statistics successfully saved to '{save.path}'.")
  }, error = function(e) {
    logger::log_error("Failed to save formatted summary statistics to '{save.path}': {e$message}")
    stop("Failed to save summary statistics: ", e$message)
  })
}
