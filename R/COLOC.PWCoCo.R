#' Perform Pairwise Conditional and Colocalization Analysis using PWCOCO
#'
#' This function performs pairwise conditional and colocalization analysis between two sets of summary statistics
#' using the PWCOCO tool. It prepares input arguments, runs the external PWCOCO command, and processes the output
#' to return colocalization results and associated files.
#'
#' @param sum_stats1 Character. Path to the first summary statistics file.
#' @param sum_stats2 Character. Path to the second summary statistics file.
#' @param n1 Numeric. Sample size for the first dataset. Must be a positive integer.
#' @param n2 Numeric. Sample size for the second dataset. Must be a positive integer.
#' @param p_cutoff1 Numeric. P-value threshold for the first dataset (default: 5e-08).
#' @param p_cutoff2 Numeric. P-value threshold for the second dataset (default: 5e-08).
#' @param top_snp Numeric. Maximum number of top SNPs to consider (default: 1e+10).
#' @param ld_window Numeric. LD window size in base pairs for analysis (default: 1e+07).
#' @param collinear Numeric. Collinearity threshold for LD pruning (default: 0.9).
#' @param maf Numeric. Minor allele frequency threshold (default: 0.1).
#' @param freq_threshold Numeric. Frequency threshold for allele matching (default: 0.2).
#' @param init_h4 Numeric. Initial value for H4 hypothesis in colocalization (default: 80).
#' @param coloc_pp Numeric vector. Posterior probability thresholds for colocalization (default: c(1e-04, 1e-04, 1e-05)).
#' @param n1_case Numeric. Number of cases for the first dataset if binary trait (default: NULL).
#' @param n2_case Numeric. Number of cases for the second dataset if binary trait (default: NULL).
#' @param chr Character or numeric. Chromosome to analyze (default: NULL, analyzes all chromosomes).
#' @param out_cond Logical. Whether to output conditional analysis results (default: TRUE).
#' @param pve1 Numeric. Proportion of variance explained for the first dataset (default: NULL).
#' @param pve2 Numeric. Proportion of variance explained for the second dataset (default: NULL).
#' @param pve_file1 Character. Path to file with PVE for the first dataset (default: NULL).
#' @param pve_file2 Character. Path to file with PVE for the second dataset (default: NULL).
#' @param threads Numeric. Number of threads to use for PWCOCO (default: 4).
#' @param snp.index Numeric. Index of SNP column in input files (default: 1).
#' @param bfile Character. Path to PLINK binary file for LD reference (default: NULL).
#' @param pwcoco Character. Path to the PWCOCO executable. Must be provided.
#' @param save_path Character. Directory path to save output files. Must be provided.
#' @param prefix Character. Prefix for output file names. Must be provided.
#' @param verbose Logical. Whether to output verbose messages for debugging (default: TRUE).
#'
#' @return A list of class \code{pwcoco} containing the following elements:
#'   \item{coloc.dat}{A data frame with colocalization results from PWCOCO.}
#'   \item{sum_stats1}{A named list of paths to conditional analysis output files for the first dataset, indexed by rsID.}
#'   \item{sum_stats2}{A named list of paths to conditional analysis output files for the second dataset, indexed by rsID.}
#'   \item{pop}{Character. The basename of the reference bfile used.}
#'
#' @details
#' This function interfaces with the PWCOCO tool to perform pairwise conditional and colocalization analysis.
#' It validates input parameters, constructs the command-line arguments, executes the PWCOCO tool, and processes
#' the output files. Ensure that the PWCOCO executable and necessary input files are accessible, and the output
#' directory has write permissions. The function logs progress and errors if \code{verbose = TRUE}.
#'
#' @examples
#' \dontrun{
#' # Load required packages
#' library(GWASkitR)
#' library(locuscomparer)
#'
#' # Perform colocalization analysis
#' PWCoCo.res <- GWASkitR::COLOC.PWCoCo(
#'   sum_stats1 = glue::glue("{datadir}/cis-p-QTL.txt"),
#'   sum_stats2 = glue::glue("{datadir}/sumstats-coloc.txt"),
#'   p_cutoff1 = 5e-08,
#'   p_cutoff2 = 5e-08,
#'   top_snp = 1e+10,
#'   ld_window = 1e+07,
#'   collinear = 0.9,
#'   maf = 0.1,
#'   freq_threshold = 0.2,
#'   init_h4 = 100,
#'   coloc_pp = c(1e-04, 1e-04, 1e-05),
#'   n1 = 35373,
#'   n2 = 318014,
#'   n1_case = NULL,
#'   n2_case = 33043,
#'   chr = 2,
#'   out_cond = TRUE,
#'   pve1 = NULL,
#'   pve2 = NULL,
#'   pve_file1 = NULL,
#'   pve_file2 = NULL,
#'   threads = 24,
#'   snp.index = 1,
#'   bfile = "your/path/to/plink/1kg.v3/EUR",
#'   pwcoco = "your/path/to/PWCoCo",
#'   save_path = "your/path/to/output/pwcoco",
#'   prefix = "pwcoco_result",
#'   verbose = TRUE
#' )
#'
#' # Prepare data for visualization
#' coloc.dat <- coloc_result_to_plot_pwcoco(
#'   x = PWCoCo.res,
#'   rs = "rs1260326",
#'   plink = "your/path/to/plink/plink.exe",
#'   bfile = "your/path/to/plink/1kg.v3/EUR"
#' )
#'
#' # Visualize results using locuscomparer
#' locuscomparer::locuscompare(
#'   in_fn1 = coloc.dat[["in_fn1"]],
#'   in_fn2 = coloc.dat[["in_fn2"]],
#'   title1 = "eQTL",
#'   title2 = "GWAS",
#'   genome = "hg19",
#'   snp = coloc.dat[["SNP"]],
#'   population = "EUR",
#'   combine = TRUE
#' )
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom logger log_info log_error log_warn
#' @importFrom data.table fread
#' @importFrom dplyr filter arrange
#' @importFrom glue glue
#'
#' @export
COLOC.PWCoCo <- function(
    sum_stats1,
    sum_stats2,
    n1,
    n2,
    p_cutoff1 = 5e-08,
    p_cutoff2 = 5e-08,
    top_snp = 1e+10,
    ld_window = 1e+07,
    collinear = 0.9,
    maf = 0.1,
    freq_threshold = 0.2,
    init_h4 = 80,
    coloc_pp = c(1e-04, 1e-04, 1e-05),
    n1_case = NULL,
    n2_case = NULL,
    chr = NULL,
    out_cond = TRUE,
    pve1 = NULL,
    pve2 = NULL,
    pve_file1 = NULL,
    pve_file2 = NULL,
    threads = 4,
    snp.index = 1,
    bfile = NULL,
    pwcoco = NULL,
    save_path = NULL,
    prefix = NULL,
    verbose = TRUE
) {

  # Initialize logging
  start.time <- Sys.time()
  if (verbose) logger::log_info("Starting pairwise conditional and colocalization analysis with PWCOCO.")

  if (!is.numeric(n1) || length(n1) != 1 || is.na(n1)) {
    logger::log_error("Parameter 'n1' must be a non-missing numeric value.")
    stop("Parameter 'n1' must be a non-missing numeric value.")
  }
  if (n1 %% 1 != 0) {
    logger::log_warn("Parameter 'n1' is not an integer (n1 = {n1}). Rounding down to nearest integer: {as.integer(n1)}.")
    n1 <- as.integer(n1)
  }
  if (n1 <= 0) {
    logger::log_error("Parameter 'n1' must be greater than 0.")
    stop("Parameter 'n1' must be greater than 0.")
  }

  if (!is.numeric(n2) || length(n2) != 1 || is.na(n2)) {
    logger::log_error("Parameter 'n2' must be a non-missing numeric value.")
    stop("Parameter 'n2' must be a non-missing numeric value.")
  }
  if (n2 %% 1 != 0) {
    logger::log_warn("Parameter 'n2' is not an integer (n2 = {n2}). Rounding down to nearest integer: {as.integer(n2)}.")
    n2 <- as.integer(n2)
  }
  if (n2 <= 0) {
    logger::log_error("Parameter 'n2' must be greater than 0.")
    stop("Parameter 'n2' must be greater than 0.")
  }

  # Log input parameters for debugging
  if (verbose) logger::log_info("Input parameters: p_cutoff1={p_cutoff1}, p_cutoff2={p_cutoff2}, ld_window={ld_window}, threads={threads}")
  if (verbose) logger::log_info("Output save path: {save_path}, prefix: {prefix}")

  # Validate input file paths
  if (verbose) logger::log_info("Validating input file paths.")
  if (!file.exists(pwcoco)) {
    error_msg <- paste0("PWCOCO executable not found at: ", pwcoco)
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  }
  if (!file.exists(sum_stats1)) {
    error_msg <- paste0("Summary statistics file 1 not found at: ", sum_stats1)
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  }
  if (!file.exists(sum_stats2)) {
    error_msg <- paste0("Summary statistics file 2 not found at: ", sum_stats2)
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  }
  if (verbose) logger::log_info("Input file paths validated successfully.")

  # Prepare output directory and check for existing files
  out <- file.path(save_path, prefix)
  coloc.file <- paste0(out, '.coloc')
  log.file <- file.path(save_path, 'pwcoco_log.txt')
  if (verbose) logger::log_info("Setting up output directory and log file at: {save_path}.")

  if (dir.exists(save_path)) {
    if (verbose) logger::log_warn("Existing output directory found at {save_path}. Deleting it.")
    unlink(save_path, recursive = TRUE)
  }
  if (!dir.exists(save_path)) {
    if (verbose) logger::log_info("Creating output directory at {save_path}.")
    dir.create(save_path, recursive = TRUE)
  }
  if (file.exists(coloc.file)) {
    error_msg <- paste0("Output colocalization file already exists at: ", coloc.file, ". Please delete it before running.")
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  }

  # Build command-line arguments for PWCOCO
  if (verbose) logger::log_info("Building PWCOCO command-line arguments.")
  args <- c(
    "--bfile", bfile,
    "--sum_stats1", sum_stats1,
    "--sum_stats2", sum_stats2,
    "--p_cutoff1", p_cutoff1,
    "--p_cutoff2", p_cutoff2,
    "--top_snp", top_snp,
    "--ld_window", ld_window,
    "--collinear", collinear,
    "--maf", maf,
    "--freq_threshold", freq_threshold,
    "--init_h4", init_h4,
    "--coloc_pp", paste(coloc_pp, collapse = " "),
    "--n1", n1,
    "--n2", n2,
    ifelse(is.null(chr), "", paste0("--chr ", chr)),
    ifelse(out_cond == FALSE, "", "--out_cond"),
    ifelse(is.null(pve1), "", paste0("--pve1 ", pve1)),
    ifelse(is.null(pve2), "", paste0("--pve2 ", pve2)),
    ifelse(is.null(pve_file1), "", paste0("--pve_file1 ", pve_file1)),
    ifelse(is.null(pve_file2), "", paste0("--pve_file2 ", pve_file2)),
    ifelse(is.null(n1_case), "", paste0("--n1_case ", n1_case)),
    ifelse(is.null(n2_case), "", paste0("--n2_case ", n2_case)),
    "--threads", threads,
    "--out", out,
    "--log", log.file
  )
  if (verbose) logger::log_info("PWCOCO arguments prepared: {paste(args, collapse = ' ')}")

  # Execute PWCOCO command
  if (verbose) logger::log_info("Executing PWCOCO analysis.")
  if (verbose) logger::log_info("Currently performing pairwise conditional and colocalization analysis using PWCoCo ...")
  result <- tryCatch({
    if (verbose) {
      system2(command = pwcoco, args = args)
    } else {
      system2(command = pwcoco, args = args, stdout = TRUE, stderr = TRUE)
    }
  }, error = function(e) {
    if (verbose) error_msg <- paste("Error running PWCOCO:", e$message)
    logger::log_error(error_msg)
    stop(error_msg)
  })

  # Check if status is numeric before proceeding
  status <- result[length(result)]
  if (verbose) logger::log_info("System command execution completed, checking exit status.")

  if (!is.numeric(status)) {
    if (verbose) logger::log_warn("Exit status is not numeric. Unable to validate command execution status. Proceeding with caution.")
  } else {
    if (status != 0) {
      if (verbose) logger::log_error("System command '{pwcoco} {paste(args, collapse = ' ')}' failed with exit status {status}. Check if the command and its arguments are correct, and verify permissions and file existence.")
      stop(sprintf("The command failed and returned exit status %d. Please check the logs for more details.", status))
    } else {
      if (verbose) logger::log_info("System command executed successfully with exit status {status}.")
    }
  }

  if (verbose) logger::log_info("PWCOCO analysis completed.")

  # Read and process colocalization results
  if (verbose) logger::log_info("Reading colocalization results from {coloc.file}.")
  if (!file.exists(coloc.file)) {
    error_msg <- paste("Colocalization output file not found at:", coloc.file)
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  }
  coloc.dat <- tryCatch({
    data.table::fread(coloc.file) %>%
      dplyr::filter(SNP1 != 'unconditioned', SNP2 != 'unconditioned') %>%
      dplyr::arrange(-log_abf_all, -H4)
  }, error = function(e) {
    error_msg <- paste("Error reading colocalization file:", e$message)
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  })
  if (verbose) logger::log_info("Colocalization results processed with {nrow(coloc.dat)} rows.")

  # Compile output data structure
  data <- list()
  data[["coloc.dat"]] <- coloc.dat
  data[["pop"]] <- basename(bfile)

  # Collect conditional analysis output files for each dataset
  if (verbose) logger::log_info("Collecting conditional analysis output files from {save_path}.")
  files <- dir(save_path)
  sumstat1.list <- list()
  sumstat2.list <- list()
  pattern1 <- glue::glue("^{prefix}\\.{basename(sum_stats1)}\\.(rs[0-9]+)\\.cojo$")
  pattern2 <- glue::glue("^{prefix}\\.{basename(sum_stats2)}\\.(rs[0-9]+)\\.cojo$")
  for (cojo in files) {
    if (grepl(pattern1, cojo)) {
      rs <- sub(pattern1, "\\1", cojo)
      sumstat1.list[[rs]] <- file.path(save_path, cojo)
      if (verbose) logger::log_info("Found conditional file for sum_stats1, rsID={rs}: {cojo}")
    }
    if (grepl(pattern2, cojo)) {
      rs <- sub(pattern2, "\\1", cojo)
      sumstat2.list[[rs]] <- file.path(save_path, cojo)
      if (verbose) logger::log_info("Found conditional file for sum_stats2, rsID={rs}: {cojo}")
    }
  }
  data[["sum_stats1"]] <- sumstat1.list
  data[["sum_stats2"]] <- sumstat2.list
  if (verbose) logger::log_info("Collected {length(sumstat1.list)} files for sum_stats1 and {length(sumstat2.list)} files for sum_stats2.")

  # Log completion and execution time
  end.time <- Sys.time()
  execution.time <- end.time - start.time
  if (verbose) logger::log_info("Returning results.")
  if (verbose) logger::log_info("Pairwise conditional and colocalization analysis completed in {round(execution.time, 2)} seconds.")
  if (verbose) logger::log_info(execution.time)

  class(data) <- c('PWCoCo', class(data))
  return(data)
}
