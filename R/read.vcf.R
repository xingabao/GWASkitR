#' Read and Parse a VCF File for GWAS Summary Data
#'
#' This function reads a Variant Call Format (VCF) file containing GWAS summary
#' statistics, extracts relevant data, processes FORMAT fields, and optionally
#' saves the parsed data for each sample. Special handling is provided for
#' log-transformed p-values (LP, interpreted as -log10 p-value).
#'
#' @param vcfFile Character string. The path to the input VCF file. The file must
#'   exist and be readable.
#' @param multiple Logical. If \code{TRUE}, parses data for all samples in the VCF
#'   file. If \code{FALSE}, only the first sample is processed. Default is \code{FALSE}.
#' @param remove_qual_filter_info Logical. If \code{TRUE}, removes the \code{QUAL},
#'   \code{FILTER}, and \code{INFO} columns from the output. Default is \code{TRUE}.
#' @param save.path Character string or \code{NULL}. Directory path to save parsed
#'   data for each sample as a tab-separated text file. If \code{NULL}, no files are
#'   saved. Default is \code{NULL}.
#' @param verbose Logical. If \code{TRUE}, prints progress and debugging messages
#'   using the \code{logger} package. Default is \code{FALSE}.
#'
#' @return If \code{multiple = FALSE}, returns a \code{data.frame} containing parsed
#'   data for the first sample. If \code{multiple = TRUE}, returns a named list of
#'   \code{data.frame} objects, with each element corresponding to a sample in the
#'   VCF file.
#'
#' @details
#' This function processes VCF files to extract GWAS summary statistics. It handles
#' the \code{FORMAT} field to split subfields into separate columns and converts
#' log-transformed p-values (LP) to standard p-values if detected. The function
#' supports both single-sample and multi-sample VCF files, with optional output
#' saving for each sample's data.
#'
#' @examples
#' \dontrun{
#' library(GWASkitR)
#'
#' # Path to example VCF file
#' vcf.file <- system.file("extdata", "ieu-a-2.vcf.gz", package = "GWASkitR")
#'
#' # Read and parse data for a single sample
#' sumstats.dt <- read.vcf(vcfFile = vcf.file, verbose = FALSE)
#'
#' # Display the first few rows
#' head(sumstats.dt)
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select all_of
#' @importFrom data.table fread setnames fwrite
#' @importFrom logger log_info log_debug log_warn log_error
#' @importFrom glue glue
#'
#' @export
read.vcf <- function(
    vcfFile,
    multiple = FALSE,
    remove_qual_filter_info = TRUE,
    save.path = NULL,
    verbose = FALSE
) {

  # Info log: Starting file loading
  if (verbose) logger::log_info("Starting to load GWAS summary data from VCF file:")
  if (verbose) logger::log_info("{vcfFile}")

  # Check file existence
  if (!file.exists(vcfFile)) {
    if (verbose) logger::log_error("File does not exist: {vcfFile}")
    stop("VCF file not found: ", vcfFile)
  }

  # Try reading VCF lines, catch errors
  lines <- tryCatch({
    data.table::fread(vcfFile, fill = TRUE, header = FALSE)$V1
  }, error = function(e) {
    if (verbose) logger::log_error("Failed to read VCF file: {vcfFile}")
    if (verbose) logger::log_error("Error message: {e$message}")
    stop(e)
  })

  # Debug log: VCF file line count
  if (verbose) logger::log_debug("Number of lines read from VCF: {length(lines)}")

  # Initialize vectors to store FORMAT lines, data lines, and header line
  format.lines <- c()
  data.lines <- c()
  header.line <- NULL

  # Iterate over all lines to separate format lines, header, and data lines
  for (l in lines) {
    if (startsWith(l, "##FORMAT")) {
      format.lines <- c(format.lines, l)
    } else if (startsWith(l, "#CHR")) {
      header.line <- l
    } else if (!startsWith(l, "##")) {
      data.lines <- c(data.lines, l)
    }
  }

  # Warn if header line is missing
  if (is.null(header.line)) {
    if (verbose) logger::log_error("Header line starting with '#CHR' not found in VCF file.")
    stop("Invalid VCF file: missing header line.")
  }

  # Debug log: header and FORMAT lines
  if (verbose) logger::log_debug("Header line: {header.line}")
  if (verbose) logger::log_debug("Number of FORMAT lines: {length(format.lines)}")

  # Process column names from the header line
  header.cols <- unlist(strsplit(header.line, "\t"))
  header.cols[1] <- "CHR"  # Ensure the first column is named "CHR"

  # Check column count
  if (length(header.cols) < 8) {
    if (verbose) logger::log_warn("Number of columns in header is suspiciously low ({length(header.cols)}).")
  }

  # Read data lines as a data.frame
  DF <- tryCatch({
    data.table::fread(text = paste(data.lines, collapse = "\n"), sep = "\t", header = FALSE)
  }, error = function(e) {
    if (verbose) logger::log_error("Failed to read VCF data lines.")
    if (verbose) logger::log_error("Error message: {e$message}")
    stop(e)
  })

  # Set column names for the data
  data.table::setnames(DF, header.cols)
  if (verbose) logger::log_debug("Column names set: {paste(header.cols, collapse=', ')}")

  # Optionally remove QUAL, FILTER, INFO columns
  if (remove_qual_filter_info) {
    cols_to_remove <- c("QUAL", "FILTER", "INFO")
    missing_cols <- setdiff(cols_to_remove, names(DF))
    if (length(missing_cols) > 0) {
      if (verbose) logger::log_warn("Some columns to remove not present in data: {paste(missing_cols, collapse=', ')}")
    }
    DF <- DF %>% dplyr::select(-any_of(cols_to_remove))
    if (verbose) logger::log_info("Removed columns: {paste(intersect(cols_to_remove, names(DF)), collapse=', ')}")
  }

  if (verbose) logger::log_info("Loaded '{basename(vcfFile)}' successfully.")

  # Determine sample columns and non-sample columns
  header.cols <- colnames(DF)
  if (!"FORMAT" %in% header.cols) {
    if (verbose) logger::log_error("FORMAT column not found in VCF data.")
    stop("FORMAT column is missing in the VCF data.")
  }
  sample.cols <- header.cols[(which(header.cols == 'FORMAT') + 1):length(header.cols)]
  albert.cols <- setdiff(header.cols, sample.cols)
  if (!multiple) { sample.cols = sample.cols[1] }

  if (verbose) logger::log_debug("Sample columns: {paste(sample.cols, collapse=', ')}")
  if (verbose) logger::log_debug("Non-sample columns: {paste(albert.cols, collapse=', ')}")

  # Check FORMAT definitions for log10 p-value field
  log10 <- FALSE
  for (ft in format.lines) {
    if (grepl("ID=LP", ft) && grepl("-log10", ft)) {
      log10 <- TRUE
      if (verbose) logger::log_info("Detected FORMAT field 'LP' as -log10(p-value). Will convert to normal p-value.")
      break
    }
  }

  data = list()
  DF.this <- data.frame()
  if ("FORMAT" %in% names(DF)) {
    # Parse the FORMAT field to get subfield names
    FTs <- unlist(strsplit(DF$FORMAT[1], ":"))
    if (length(FTs) == 0) {
      if (verbose) logger::log_warn("FORMAT field is empty or malformed.")
    } else {
      for (sample.col in sample.cols) {
        if (verbose) logger::log_info("Parsing sample column: {sample.col}")
        DF.this <- DF %>% dplyr::select(all_of(albert.cols))
        SEs <- lapply(DF[[sample.col]], function(x) unlist(strsplit(x, ":")))
        # Assign each subfield as a new column
        for (i in seq_along(FTs)) {
          DF.this[[FTs[i]]] <- sapply(SEs, function(x) if (length(x) >= i) x[i] else NA)
        }
        # Convert LP (-log10 p-value) to PVAL if necessary
        if (log10) {
          if (verbose) logger::log_info("Converting -log10 p-value (LP) to normal p-value for sample: {sample.col}")
          DF.this$PVAL <- 10 ^ (-1 * as.numeric(DF.this$LP))
        } else if ("LP" %in% names(DF.this)) {
          if (verbose) logger::log_info("Detected LP field, using it as normal p-value for sample: {sample.col}")
          DF.this$PVAL <- as.numeric(DF.this$LP)
        }
        # Remove FORMAT and sample columns from output
        if ("FORMAT" %in% names(DF.this)) {
          DF.this$FORMAT <- NULL
        }
        # Auto-convert column types using custom function
        DF.this <- tryCatch({
          auto_type_convert(df = DF.this, n_check = 1000, verbose = verbose)
        }, error = function(e) {
          if (verbose) logger::log_warn("auto_type_convert failed for sample {sample.col}: {e$message}")
          DF.this
        })
        data[[sample.col]] = DF.this
        if (verbose) logger::log_debug("Finished parsing sample: {sample.col}")
      }
    }
  }

  # Save each sample's data to file if save.path is specified
  if (length(data) > 0 && !is.null(save.path)) {
    if (!dir.exists(paths = save.path)) {
      dir.create(path = save.path, recursive = TRUE)
      if (verbose) logger::log_info("Created output directory: {save.path}")
    }
    for (dataname in names(data)) {
      out.file <- file.path(save.path, glue::glue('{dataname}.txt'))
      data.table::fwrite(data[[dataname]], file = out.file, sep = "\t", row.names = FALSE)
      if (verbose) logger::log_info("Saved parsed table for sample '{dataname}' to: {out.file}")
    }
  }

  # Return a data.frame for a single sample, otherwise a list of data.frames
  if (length(data) == 1) {
    if (verbose) logger::log_info("Returning parsed data for single sample.")
    return(data[[1]])
  }
  if (verbose) logger::log_info("Returning parsed data for {length(data)} samples as a list.")
  return(data)
}
