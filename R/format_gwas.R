#' Format GWAS Summary Statistics
#'
#' A generic function to read, validate, select, and rename columns in GWAS summary
#' statistics files. It supports both VCF and tabular formats, ensuring required
#' columns exist, applying user-specified renaming, optionally saving the processed
#' table, and logging each step using the \code{logger} package.
#'
#' @param sumstats An object to be formatted. For the default method, this can be a character
#'   string (path to a summary statistics file in VCF or tabular format) or a
#'   \code{data.frame}/\code{data.table} object containing the data.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return A \code{tibble} object containing the formatted summary statistics with
#'   selected and renamed columns.
#'
#' @details
#' This is a generic function for formatting GWAS summary statistics. The default
#' method processes data from VCF or tabular files, or directly from a data frame.
#' It validates the presence of specified input columns, selects only the required
#' columns, and renames them to user-defined output names. VCF files are read using
#' an internal \code{read.vcf} function, while tabular files are read using
#' \code{data.table::fread}. The function logs all steps and can save the formatted
#' data to disk if a \code{save.path} is provided.
#'
#' @examples
#' \dontrun{
#' library(GWASkitR)
#'
#' # Example 1: Format a tabular file
#' sumstats.file <- system.file("extdata", "GWAS-test.txt", package = "GWASkitR")
#' sumstats.dt <- format_gwas(
#'   sumstats = sumstats.file,
#'   col.chr = "chr",
#'   col.pos = "pos",
#'   col.snp = "SNP",
#'   col.ref = "other_allele",
#'   col.alt = "effect_allele",
#'   col.metrics = "beta",
#'   col.eaf = "eaf",
#'   col.se = "se",
#'   col.pval = "pval",
#'   col.nsample = "samplesize",
#'   col.chr.set = "CHR",
#'   col.pos.set = "POS",
#'   col.snp.set = "SNP",
#'   col.ref.set = "other_allele",
#'   col.alt.set = "effect_allele",
#'   col.metrics.set = "beta",
#'   col.eaf.set = "eaf",
#'   col.se.set = "se",
#'   col.pval.set = "pval",
#'   col.nsample.set = "nsamples",
#'   verbose = TRUE
#' )
#' head(sumstats.dt)
#'
#' # Example 2: Format a data frame
#' sumstats.dt <- format_gwas(
#'   sumstats = GWASkitR::sumstats,
#'   col.chr = "CHR",
#'   col.pos = "POS",
#'   col.snp = "SNP",
#'   col.ref = "other_allele",
#'   col.alt = "effect_allele",
#'   col.metrics = "beta",
#'   col.eaf = "eaf",
#'   col.se = "se",
#'   col.pval = "pval",
#'   col.chr.set = "CHR",
#'   col.pos.set = "POS",
#'   col.snp.set = "SNP",
#'   col.ref.set = "other_allele",
#'   col.alt.set = "effect_allele",
#'   col.metrics.set = "beta",
#'   col.eaf.set = "eaf",
#'   col.se.set = "se",
#'   col.pval.set = "pval",
#'   verbose = TRUE
#' )
#' head(sumstats.dt)
#'
#' # Example 3: Format a VCF file
#' vcf.file <- system.file("extdata", "ieu-a-2.vcf.gz", package = "GWASkitR")
#' sumstats.dt <- format_gwas(
#'   sumstats = vcf.file,
#'   col.chr = "CHR",
#'   col.pos = "POS",
#'   col.snp = "ID",
#'   col.ref = "REF",
#'   col.alt = "ALT",
#'   col.metrics = "ES",
#'   col.eaf = "AF",
#'   col.se = "SE",
#'   col.pval = "PVAL",
#'   col.nsample = "SS",
#'   col.chr.set = "CHR",
#'   col.pos.set = "POS",
#'   col.snp.set = "SNP",
#'   col.ref.set = "other_allele",
#'   col.alt.set = "effect_allele",
#'   col.metrics.set = "beta",
#'   col.eaf.set = "eaf",
#'   col.se.set = "se",
#'   col.pval.set = "pval",
#'   col.nsample.set = "nsample",
#'   verbose = TRUE
#' )
#' head(sumstats.dt)
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom data.table fread setnames
#' @importFrom dplyr select all_of as_tibble
#' @importFrom logger log_info log_warn log_error
#' @export
format_gwas <- function(sumstats, ...) {
  UseMethod("format_gwas")
}


#' @rdname format_gwas
#' @param col.chr Character string or \code{NULL}. Name of the chromosome column
#'   in the input data.
#' @param col.pos Character string or \code{NULL}. Name of the position column in
#'   the input data.
#' @param col.snp Character string or \code{NULL}. Name of the SNP identifier column
#'   in the input data.
#' @param col.ref Character string or \code{NULL}. Name of the reference allele
#'   column in the input data.
#' @param col.alt Character string or \code{NULL}. Name of the alternative allele
#'   column in the input data.
#' @param col.metrics Character string or \code{NULL}. Name of the effect size column
#'   in the input data (e.g., beta, OR, ES, Z).
#' @param col.eaf Character string or \code{NULL}. Name of the effect allele frequency
#'   column in the input data.
#' @param col.se Character string or \code{NULL}. Name of the standard error column
#'   in the input data.
#' @param col.pval Character string or \code{NULL}. Name of the p-value column in
#'   the input data.
#' @param col.ncase Character string or \code{NULL}. Name of the case sample size
#'   column in the input data.
#' @param col.nsample Character string or \code{NULL}. Name of the total sample size
#'   column in the input data.
#' @param col.chr.set Character string. Desired output column name for chromosome.
#'   Default is \code{"CHR"}.
#' @param col.pos.set Character string. Desired output column name for position.
#'   Default is \code{"POS"}.
#' @param col.snp.set Character string. Desired output column name for SNP.
#'   Default is \code{"SNP"}.
#' @param col.ref.set Character string. Desired output column name for reference allele.
#'   Default is \code{"other_allele"}.
#' @param col.alt.set Character string. Desired output column name for alternative allele.
#'   Default is \code{"effect_allele"}.
#' @param col.metrics.set Character string. Desired output column name for effect size.
#'   Default is the value of \code{col.metrics}.
#' @param col.eaf.set Character string. Desired output column name for effect allele frequency.
#'   Default is \code{"eaf"}.
#' @param col.se.set Character string. Desired output column name for standard error.
#'   Default is \code{"se"}.
#' @param col.pval.set Character string. Desired output column name for p-value.
#'   Default is \code{"pval"}.
#' @param col.ncase.set Character string. Desired output column name for case sample size.
#'   Default is \code{"ncase"}.
#' @param col.nsample.set Character string. Desired output column name for total sample size.
#'   Default is \code{"nsample"}.
#' @param save.path Character string or \code{NULL}. File path to save the formatted
#'   summary statistics as a tab-separated file. If \code{NULL}, no file is saved.
#'   Default is \code{NULL}.
#' @param verbose Logical. If \code{TRUE}, prints progress and debugging messages
#'   using the \code{logger} package. Default is \code{TRUE}.
#'
#' @export
format_gwas.default <- function(
    sumstats,
    col.chr = NULL,
    col.pos = NULL,
    col.snp = NULL,
    col.ref = NULL,
    col.alt = NULL,
    col.metrics = NULL,
    col.eaf = NULL,
    col.se = NULL,
    col.pval = NULL,
    col.ncase = NULL,
    col.nsample = NULL,
    col.chr.set = 'CHR',
    col.pos.set = 'POS',
    col.snp.set = 'SNP',
    col.ref.set = 'other_allele',
    col.alt.set = 'effect_allele',
    col.metrics.set = col.metrics,
    col.eaf.set = 'eaf',
    col.se.set = 'se',
    col.pval.set = 'pval',
    col.ncase.set = 'ncase',
    col.nsample.set = 'nsample',
    save.path = NULL,
    verbose = TRUE
) {

  sumstats.o <- sumstats
  check.file.path(file.path = save.path, suffix = c('tsv', 'txt'), verbose = verbose, stop_on_error = TRUE)

  # Start of the function
  if (verbose) start.time <- Sys.time()
  if (verbose) logger::log_info("Starting GWAS summary statistics formatting.")

  if (TRUE) {
    # If input is a file path, read the file
    if (is.character(sumstats) && file.exists(sumstats)) {
      if (any(endsWith(sumstats, suffix = c('vcf', 'vcf.gz')))) {
        sumstats <- read.vcf(vcfFile = sumstats, multiple = FALSE, remove_qual_filter_info = TRUE, save.path = NULL, verbose = FALSE, nrows = 1000)
      } else {
        sumstats <- as.data.frame(data.table::fread(sumstats, showProgress = FALSE, nrows = 1000))
      }
    } else if (is.character(sumstats)) {
      if (verbose) logger::log_error("The file '{sumstats}' does not exist.")
      stop("Input summary statistics file does not exist: ", sumstats)
    }
  }

  # Collect the required columns
  .cols <- c(col.chr, col.pos, col.snp, col.ref, col.alt, col.metrics, col.eaf, col.se, col.pval, col.ncase, col.nsample)

  # Remove NULLs from .cols to avoid checking for missing NULL columns
  .cols <- .cols[!is.null(.cols)]

  # Validate that required columns exist in the data
  missing.cols <- setdiff(.cols, names(sumstats))
  if (length(missing.cols) > 0) {
    if (verbose) logger::log_warn("These are the column names in the input data: {paste(names(sumstats), collapse = ', ')}")
    if (verbose) logger::log_info("Preview of the first few rows of input data:")
    print(utils::head(sumstats))
    if (verbose) logger::log_error("Missing required columns: {paste(missing.cols, collapse = ', ')}")
    stop("Missing columns: ", paste(missing.cols, collapse = ', '))
  } else if (length(.cols) == 0) {
    if (verbose) logger::log_warn("No columns exist in the data frame, please check.")
    if (verbose) logger::log_warn("These are the column names in the input data: {paste(names(sumstats), collapse = ', ')}")
    if (verbose) logger::log_info("Preview of the first few rows of input data:")
    print(utils::head(sumstats))
  }

  # If input is a file path, read the file
  sumstats <- sumstats.o
  if (is.character(sumstats) && file.exists(sumstats)) {
    if (any(endsWith(sumstats, suffix = c('vcf', 'vcf.gz')))) {
      if (verbose) logger::log_info("Detected VCF file: '{sumstats}'. Reading as VCF.")
      sumstats <- read.vcf(vcfFile = sumstats, multiple = FALSE, remove_qual_filter_info = TRUE, save.path = NULL, verbose = verbose)
      if (verbose) logger::log_info("VCF file loaded successfully.")
    } else {
      if (verbose) logger::log_info("Detected local file: '{sumstats}'. Reading as table.")
      sumstats <- as.data.frame(data.table::fread(sumstats, showProgress = verbose))
      if (verbose) logger::log_info("Local file loaded successfully.")
    }
  } else if (is.character(sumstats)) {
    if (verbose) logger::log_error("The file '{sumstats}' does not exist.")
    stop("Input summary statistics file does not exist: ", sumstats)
  } else {
    if (verbose) logger::log_info("Input provided as data.frame or data.table. Proceeding.")
  }

  if (verbose) logger::log_info("Selecting required columns: {paste(.cols, collapse=', ')}")
  sumstats.dt <- sumstats %>% dplyr::select(dplyr::all_of(.cols))

  # Rename columns if necessary, logging each operation
  if (is.valid(col.chr, col.chr.set)) safe.setnames(sumstats.dt, col.chr, col.chr.set, verbose = verbose)
  if (is.valid(col.pos, col.pos.set)) safe.setnames(sumstats.dt, col.pos, col.pos.set, verbose = verbose)
  if (is.valid(col.snp, col.snp.set)) safe.setnames(sumstats.dt, col.snp, col.snp.set, verbose = verbose)
  if (is.valid(col.ref, col.ref.set)) safe.setnames(sumstats.dt, col.ref, col.ref.set, verbose = verbose)
  if (is.valid(col.alt, col.alt.set)) safe.setnames(sumstats.dt, col.alt, col.alt.set, verbose = verbose)
  if (is.valid(col.metrics, col.metrics.set)) safe.setnames(sumstats.dt, col.metrics, col.metrics.set, verbose = verbose)
  if (is.valid(col.eaf, col.eaf.set)) safe.setnames(sumstats.dt, col.eaf, col.eaf.set, verbose = verbose)
  if (is.valid(col.se, col.se.set)) safe.setnames(sumstats.dt, col.se, col.se.set, verbose = verbose)
  if (is.valid(col.pval, col.pval.set)) safe.setnames(sumstats.dt, col.pval, col.pval.set, verbose = verbose)
  if (is.valid(col.ncase, col.ncase.set)) safe.setnames(sumstats.dt, col.ncase, col.ncase.set, verbose = verbose)
  if (is.valid(col.nsample, col.nsample.set)) safe.setnames(sumstats.dt, col.nsample, col.nsample.set, verbose = verbose)

  if (verbose) logger::log_info("Column selection and renaming complete.")

  # Optionally, save the formatted data to file
  if (!is.null(save.path)) {
    save.sumstats(sumstats.dt, save.path)
  } else {
    if (verbose) logger::log_info("No save.path provided. Data will not be written to disk.")
  }

  # Log completion and execution time
  if (verbose) execution.time <- Sys.time() - start.time
  if (verbose) logger::log_info("GWAS summary statistics formatting completed in {round(execution.time, 2)} minutes")

  gc()

  return(sumstats.dt %>% dplyr::as_tibble())
}


#' Format GWAS Summary Statistics for MiXeR Analysis
#'
#' Formats GWAS summary statistics specifically for compatibility with MiXeR analysis
#' by standardizing column names to a predefined set required by MiXeR. It builds upon
#' the default \code{format_gwas} method to ensure the output meets MiXeR's requirements.
#'
#' @param col.chr Character string or \code{NULL}. Name of the chromosome column in
#'   the input data.
#' @param col.pos Character string or \code{NULL}. Name of the position column in
#'   the input data.
#' @param col.snp Character string or \code{NULL}. Name of the SNP identifier column
#'   in the input data.
#' @param col.ref Character string or \code{NULL}. Name of the reference allele
#'   column in the input data.
#' @param col.alt Character string or \code{NULL}. Name of the alternative allele
#'   column in the input data.
#' @param col.metrics Character string or \code{NULL}. Name of the effect size column
#'   in the input data (e.g., beta, OR).
#' @param col.eaf Character string or \code{NULL}. Name of the effect allele frequency
#'   column in the input data.
#' @param col.se Character string or \code{NULL}. Name of the standard error column
#'   in the input data.
#' @param col.pval Character string or \code{NULL}. Name of the p-value column in
#'   the input data.
#' @param col.ncase Character string or \code{NULL}. Name of the case sample size
#'   column in the input data.
#' @param col.nsample Character string or \code{NULL}. Name of the total sample size
#'   column in the input data.
#' @param col.metrics.set Character string or \code{NULL}. Desired output column name
#'   for effect size. Must be either \code{"BETA"} or \code{"OR"}. If \code{NULL}, uses
#'   the value of \code{col.metrics}.
#' @param save.path Character string or \code{NULL}. File path to save the formatted
#'   summary statistics as a tab-separated file. If \code{NULL}, no file is saved.
#'   Default is \code{NULL}.
#' @param verbose Logical. If \code{TRUE}, prints progress and debugging messages
#'   using the \code{logger} package. Default is \code{TRUE}.
#'
#' @return A \code{tibble} object containing the formatted summary statistics with
#'   column names standardized for MiXeR analysis.
#'
#' @details
#' This method formats GWAS summary statistics to meet the specific requirements of
#' MiXeR analysis. It uses fixed output column names (e.g., \code{"CHROM"} for chromosome,
#' \code{"A1"} for reference allele) and validates the effect size metric to ensure it
#' is either \code{"BETA"} or \code{"OR"}. The function delegates the core formatting
#' to \code{format_gwas.default} while enforcing MiXeR-specific constraints.
#'
#' @examples
#' \dontrun{
#' # Load the data
#' sumstats.file <- system.file("extdata", "GWAS-test.txt", package = "GWASkitR")
#' sumstats.dt <- data.table::fread(sumstats.file)
#'
#' # Assign the "mixer" class to sumstats.dt
#' class(sumstats.dt) <- c("mixer", class(sumstats.dt))
#'
#' # Call format_gwas with the mixer method
#' mixer.dt <- format_gwas(
#'   sumstats = sumstats.dt,
#'   col.chr = "chr",
#'   col.pos = "pos",
#'   col.snp = "SNP",
#'   col.ref = "other_allele",
#'   col.alt = "effect_allele",
#'   col.metrics = "beta",
#'   col.eaf = "eaf",
#'   col.se = "se",
#'   col.pval = "pval",
#'   col.metrics.set = "BETA",
#'   verbose = TRUE
#' )
#'
#' # Check the result
#' head(mixer.dt)
#' }
#'
#' @importFrom logger log_info log_debug log_warn log_error
#'
#' @export
format_gwas.mixer <- function(
    sumstats,
    col.chr = NULL,
    col.pos = NULL,
    col.snp = NULL,
    col.ref = NULL,
    col.alt = NULL,
    col.metrics = NULL,
    col.eaf = NULL,
    col.se = NULL,
    col.pval = NULL,
    col.ncase = NULL,
    col.nsample = NULL,
    col.metrics.set = NULL,
    save.path = NULL,
    verbose = TRUE
) {

  # Log the start of the formatting process at INFO level
  if (verbose) {
    logger::log_info("Starting to format summary statistics for MiXeR analysis.")
  }

  # Log debugging information about input data dimensions at DEBUG level
  if (verbose) {
    logger::log_debug("Input data dimensions: {nrow(sumstats)} rows, {ncol(sumstats)} columns.")
  }

  # Define fixed output column names for MiXeR compatibility
  col.chr.set <- 'CHROM'
  col.pos.set <- 'POS'
  col.snp.set <- 'SNP'
  col.ref.set <- 'A1'
  col.alt.set <- 'A2'
  col.eaf.set <- 'FRQ'
  col.se.set <- 'SE'
  col.pval.set <- 'PVAL'
  col.ncase.set <- NULL
  col.nsample.set <- NULL

  # Validate the effect size metric (BETA or OR) if provided
  if (!is.null(col.metrics.set)) {
    col.metrics.set = toupper(col.metrics.set)
    if (grepl("^(OR|BETA)$", col.metrics.set, ignore.case = TRUE)) {
      if (verbose) {
        logger::log_info("Effect size metric '{col.metrics.set}' is valid.")
      }
    } else {
      if (verbose) {
        logger::log_error("Invalid effect size metric: '{col.metrics.set}'. Must be 'BETA' or 'OR'.")
      }
      stop("Invalid effect size metric provided. Must be either 'BETA' or 'OR'.")
    }
  } else {
    if (verbose) {
      logger::log_warn("No effect size metric specified. Proceeding with default or user-provided column mapping.")
    }
  }

  # Log the initiation of the default formatting process at INFO level
  if (verbose) {
    logger::log_info("Calling format_gwas.default to standardize column names.")
  }

  # Call the default formatting function with all parameters
  result <- format_gwas.default(
    sumstats = sumstats,
    col.chr = col.chr,
    col.pos = col.pos,
    col.snp = col.snp,
    col.ref = col.ref,
    col.alt = col.alt,
    col.metrics = col.metrics,
    col.eaf = col.eaf,
    col.se = col.se,
    col.pval = col.pval,
    col.ncase = col.ncase,
    col.nsample = col.nsample,
    col.chr.set = col.chr.set,
    col.pos.set = col.pos.set,
    col.snp.set = col.snp.set,
    col.ref.set = col.ref.set,
    col.alt.set = col.alt.set,
    col.metrics.set = col.metrics.set,
    col.eaf.set = col.eaf.set,
    col.se.set = col.se.set,
    col.pval.set = col.pval.set,
    col.ncase.set = col.ncase.set,
    col.nsample.set = col.nsample.set,
    save.path = save.path,
    verbose = verbose
  )

  # Log successful completion of formatting at INFO level
  if (verbose) { logger::log_info("Formatting completed successfully for MiXeR analysis.") }

  return(result)
}

