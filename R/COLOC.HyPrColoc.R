#' Multi-trait Colocalization Analysis via HyPrColoc
#'
#' Conducts colocalization analysis across multiple GWAS/QTL datasets using the
#' \pkg{hyprcoloc} package. The function harmonizes SNPs, constructs beta and
#' standard error matrices, and runs colocalization with customizable thresholds.
#' It is suitable for regional multi-trait colocalization (e.g., GWAS vs QTL).
#'
#' @param sumstats A named list of data frames, each containing summary statistics
#'   for one dataset. Each data frame must include columns: \code{"SNP"}, \code{"CHR"},
#'   \code{"POS"}, \code{"beta"}, and \code{"se"}.
#' @param reg.thresh Numeric. The regional colocalization posterior probability threshold
#'   (default: \code{0.7}).
#' @param align.thresh Numeric. The alignment threshold for colocalization (default: \code{0.7}).
#' @param prior.c Numeric. The prior probability of colocalization (default: \code{0.02}).
#' @param verbose Logical. If \code{TRUE}, log progress messages (default: \code{TRUE}).
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{coloc.dat}}{Results from \code{\link[hyprcoloc]{hyprcoloc}} analysis.}
#'     \item{\code{assoc}}{Data frame: SNP info and Z-scores for each dataset.}
#'     \item{\code{betas}}{Beta matrix (SNPs x datasets).}
#'     \item{\code{ses}}{Standard error matrix (SNPs x datasets).}
#'     \item{\code{traits}}{Traits/dataset names.}
#'     \item{\code{snps}}{Common SNPs analyzed.}
#'     \item{\code{snp.sig}}{Colocalization-candidate SNP(s).}
#'   }
#'
#' @details
#' All datasets must contain the required columns and the same set of SNPs in the region.
#' Non-overlapping SNPs are excluded. Internally, the function logs key steps and
#' checks input validity.
#'
#' @section Typical Workflow:
#' \enumerate{
#'   \item Perform colocalization analysis with this function.
#'   \item Prepare data for visualization (e.g., using \code{coloc_result_to_plot()}).
#'   \item Visualize regional colocalization (e.g., \code{geni.plots::fig_region_stack()}).
#'   \item Optionally, run sensitivity analysis (\code{hyprcoloc::sensitivity.plot()}).
#' }
#'
#' @examples
#' \dontrun{
#' # Step 0. Load required packages
#' library(GWASkitR)
#' library(locuscomparer)
#'
#' # Step 1. Perform colocalization analysis
#' HyPrColoc.res <- GWASkitR::COLOC.HyPrColoc(
#'   sumstats = list(
#'     pQTLdat = GWASkitR::pQTLdat,
#'     pQTLdat2 = GWASkitR::pQTLdat,
#'     GWASdat = GWASkitR::sumstats,
#'     GWASdat2 = GWASkitR::sumstats
#'   ),
#'   reg.thresh = 0.7,
#'   align.thresh = 0.7,
#'   prior.c = 0.02
#' )
#'
#' # Step 2. Data preparation for visualization
#' coloc.dat <- coloc_result_to_plot(
#'   x = HyPrColoc.res,
#'   plink = "your/path/to/plink/plink.exe",
#'   bfile = "your/path/to/plink/1kg.v3/EUR"
#' )
#'
#' # Step 3. Visualize the colocalization analysis results
#' geni.plots::fig_region_stack(
#'   data = coloc.dat[['geni.dat']],
#'   traits = coloc.dat[['traits']],
#'   corr = coloc.dat[['ld_matrix']],
#'   build = 37,
#'   title_center = TRUE
#' )
#'
#' # Step 4. Sensitivity plot
#' hyprcoloc::sensitivity.plot(
#'   effect.est = coloc.dat[['betas']],
#'   effect.se = coloc.dat[['ses']],
#'   trait.names = coloc.dat[['traits']],
#'   snp.id = coloc.dat[['snps']],
#'   reg.thresh = 0.7,
#'   align.thresh = 0.7,
#'   prior.c = 0.02
#' )
#' }
#'
#' @author Abao Xing (\email{albertxn7@@gmail.com})
#'
#' @importFrom magrittr %>%
#' @importFrom logger log_info log_warn log_error
#' @importFrom hyprcoloc hyprcoloc
#' @export
COLOC.HyPrColoc <- function(
    sumstats,
    reg.thresh = 0.7,
    align.thresh = 0.7,
    prior.c = 0.02,
    verbose = TRUE
) {

  # Initialize logging
  start.time <- Sys.time()
  if (verbose) logger::log_info("Starting colocalization analysis with HyPrColoc.")

  # Log input parameters for debugging purposes
  if (verbose) logger::log_info("Input parameters: reg.thresh={reg.thresh}, align.thresh={align.thresh}, prior.c={prior.c}")
  if (verbose) logger::log_info("Number of datasets provided: {length(sumstats)}")

  # Identify common SNPs across all datasets
  if (verbose) logger::log_info("Identifying common SNPs across datasets.")
  snplists <- lapply(sumstats, function(x) x$SNP)
  snps <- Reduce(intersect, snplists)
  snps <- unique(snps)
  if (verbose) logger::log_info("Found {length(snps)} common SNPs across all datasets.")

  # Validate required columns in each dataset
  required.columns <- c("SNP", "CHR", "POS", "beta", "se")
  if (verbose) logger::log_info("Validating required columns in datasets.")
  for (dataset_name in names(sumstats)) {
    dataset <- sumstats[[dataset_name]]
    missing_columns <- required.columns[!required.columns %in% colnames(dataset)]
    if (length(missing_columns) > 0) {
      error_msg <- paste("Dataset", dataset_name, "is missing required columns:", paste(missing_columns, collapse = ", "))
      if (verbose) logger::log_error(error_msg)
      stop(error_msg)
    }
  }
  if (verbose) logger::log_info("All datasets contain required columns.")

  # Create beta matrix for common SNPs
  if (verbose) logger::log_info("Extracting beta values for common SNPs.")
  betas <- data.frame(row.names = snps)
  for (dataset_name in names(sumstats)) {
    dataset <- sumstats[[dataset_name]]
    # Filter dataset to include only common SNPs
    filtered_dataset <- dataset[dataset$SNP %in% snps, ]
    beta_values <- filtered_dataset$beta
    names(beta_values) <- filtered_dataset$SNP
    betas[[dataset_name]] <- beta_values[row.names(betas)]
  }
  betas <- as.matrix(betas)
  if (verbose) logger::log_info("Beta matrix created with dimensions: {dim(betas)[1]} x {dim(betas)[2]}")

  # Create standard error (se) matrix for common SNPs
  if (verbose) logger::log_info("Extracting standard error values for common SNPs.")
  ses <- data.frame(row.names = snps)
  for (dataset_name in names(sumstats)) {
    dataset <- sumstats[[dataset_name]]
    filtered_dataset <- dataset[dataset$SNP %in% snps, ]
    se_values <- filtered_dataset$se
    names(se_values) <- filtered_dataset$SNP
    ses[[dataset_name]] <- se_values[row.names(ses)]
  }
  ses <- as.matrix(ses)
  if (verbose) logger::log_info("Standard error matrix created with dimensions: {dim(ses)[1]} x {dim(ses)[2]}")

  # Prepare trait names and SNP IDs for HyPrColoc
  traits <- colnames(betas)
  rsid <- rownames(betas)
  if (verbose) logger::log_info("Traits for colocalization: {paste(traits, collapse = ', ')}")

  # Perform HyPrColoc colocalization analysis
  if (verbose) logger::log_info("Running HyPrColoc colocalization analysis.")
  coloc.dat <- tryCatch({
    hyprcoloc::hyprcoloc(
      betas, ses,
      trait.names = traits,
      snp.id = rsid,
      reg.thresh = reg.thresh,
      align.thresh = align.thresh,
      prior.c = prior.c
    )
  }, error = function(e) {
    if (verbose) logger::log_error("Error in HyPrColoc analysis: {e$message}")
    stop(e)
  })
  snp.sig <- coloc.dat$results$candidate_snp
  if (verbose) logger::log_info("Colocalization analysis completed. Candidate SNP(s): {ifelse(length(snp.sig) > 0, paste(snp.sig, collapse = ', '), 'None')}")

  # Build association data frame with marker info and Z-scores
  if (verbose) logger::log_info("Building association data frame with Z-scores.")
  assoc <- data.frame(marker = snps)
  first_dataset_filtered <- sumstats[[1]][sumstats[[1]]$SNP %in% snps, ]
  assoc$chr <- first_dataset_filtered$CHR[match(assoc$marker, first_dataset_filtered$SNP)]
  assoc$pos <- first_dataset_filtered$POS[match(assoc$marker, first_dataset_filtered$SNP)]

  # Calculate and add Z-scores for each dataset
  for (i in seq_along(sumstats)) {
    dataset <- sumstats[[i]]
    filtered_dataset <- dataset[dataset$SNP %in% snps, ]
    z_values <- filtered_dataset$beta / filtered_dataset$se
    assoc[[paste0("z_", i)]] <- z_values[match(assoc$marker, filtered_dataset$SNP)]
  }
  if (verbose) logger::log_info("Z-scores calculated for {length(sumstats)} datasets.")

  # Compile results into a list
  if (verbose) logger::log_info("Compiling results.")
  data <- list()
  data[["coloc.dat"]] <- coloc.dat
  data[["assoc"]] <- assoc
  data[["betas"]] <- betas
  data[["ses"]] <- ses
  data[["traits"]] <- traits
  data[["snps"]] <- snps
  data[["snp.sig"]] <- snp.sig

  # Log completion and execution time
  end.time <- Sys.time()
  execution.time <- end.time - start.time
  if (verbose) logger::log_info("Returning results.")
  if (verbose) logger::log_info("Colocalization analysis completed in {round(execution.time, 2)} seconds.")
  if (verbose) logger::log_info(execution.time)

  class(data) <- c('HyPrColoc', class(data))
  return(data)
}
