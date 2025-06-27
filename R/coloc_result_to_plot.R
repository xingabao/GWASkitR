#' Prepare Colocalization Results for Plotting
#'
#' This generic function processes colocalization analysis results from various methods
#' (e.g., ABF, HyPrColoc, PWCOCO, SuSiE) and prepares them for downstream plotting
#' (e.g., with locuscomparer). It extracts necessary data such as SNP sets, LD matrices,
#' and association statistics, ensuring all requirements for visualization are met.
#'
#' @param x A list containing colocalization results from a specific method (e.g., output from \code{COLOC.ABF}).
#' @param ... Additional arguments passed to the specific method implementation.
#'
#' @return A list containing processed data for plotting, typically including:
#'   \itemize{
#'     \item Input data frames for traits (e.g., SNP and p-value pairs).
#'     \item LD matrix for the intersecting SNPs.
#'     \item Merged association data for visualization.
#'     \item Lead SNP or significant SNPs for focus in plots.
#'   }
#'
#' @examples
#' \dontrun{
#'   # Assume ABF.res is the output from COLOC.ABF
#'   plot_dat <- coloc_result_to_plot(
#'     x = COLOC.ABF.res,
#'     plink = "/path/to/plink/plink.exe",
#'     bfile = "/path/to/plink/1kg.v3/EUR"
#'   )
#' }
#'
#' @importFrom logger log_info log_error
#' @export
coloc_result_to_plot <- function(x, ...) {
  UseMethod("coloc_result_to_plot")
}

#' Prepare ABF Colocalization Results for Plotting
#'
#' @description
#' Processes colocalization results from \code{COLOC.ABF} for plotting, including
#' extracting the lead SNP, validating inputs, and computing a local LD matrix.
#'
#' @param x List. Output from \code{COLOC.ABF}, containing summary, results, and input SNPs.
#' @param plink Character. Path to the PLINK executable for LD computation.
#' @param bfile Character. Path prefix for PLINK binary reference files (e.g., 1000G EUR).
#' @param verbose Logical. Whether to print detailed logging messages during execution. Default: TRUE.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{in_fn1}{Data frame for trait 1 (SNP, p-value).}
#'   \item{in_fn2}{Data frame for trait 2 (SNP, p-value).}
#'   \item{ld_matrix}{SNP-by-SNP LD matrix.}
#'   \item{geni.dat}{Merged association data for plotting.}
#'   \item{SNP}{The lead SNP for colocalization.}
#'   \item{pop}{Population identifier extracted from the basename of \code{bfile}.}
#' }
#'
#' @details
#' This function processes the output of \code{COLOC.ABF} to prepare data for visualization.
#' It extracts the lead SNP, validates input data, and computes a local linkage disequilibrium (LD)
#' matrix using PLINK. The function ensures that the association data matches the LD matrix for
#' accurate plotting.
#'
#' @examples
#' \dontrun{
#' # Assume ABF.res is the output from COLOC.ABF
#' coloc.dat <- coloc_result_to_plot.abf(
#'   x = ABF.res,
#'   plink = "/your/path/to/plink/plink.exe",
#'   bfile = "/your/path/to/plink/plink/1kg.v3/EUR"
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
#'   population = coloc.dat[["pop"]],
#'   combine = TRUE
#' )
#' }
#'
#' @seealso \code{\link{COLOC.ABF}}
#'
#' @importFrom glue glue
#' @importFrom ieugwasr ld_matrix_local
#' @importFrom dplyr filter
#' @importFrom logger log_info log_error
#' @export
coloc_result_to_plot.abf <- function(x, plink = NULL, bfile = NULL, verbose = TRUE) {
  # Validate required arguments
  if (is.null(plink) || is.null(bfile)) {
    if (verbose) logger::log_error("Both 'plink' and 'bfile' parameters must be provided for LD matrix calculation.")
    stop("Both 'plink' and 'bfile' are required.")
  }

  # Extract lead SNP, coloc results, and merged association data
  rs <- x[['snp.sig']]
  coloc.dat <- x[['coloc.dat']]
  assoc <- x[['assoc']]
  if (verbose) logger::log_info("Preparing plotting data for lead SNP: {rs}")

  # Validate lead SNP format
  if (!grepl("^rs[0-9]+$", rs)) {
    error_msg <- glue::glue("Lead SNP '{rs}' does not follow the correct format (should be 'rs' followed by numbers, e.g., rs1234175).")
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  }
  # Validate that lead SNP is present in coloc results
  if (!rs %in% coloc.dat$results$snp) {
    error_msg <- glue::glue("Lead SNP '{rs}' not found in coloc results. Please check your input.")
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  }

  # Retrieve input SNP-pval tables
  in_fn1 <- x[['in_fn1']]
  in_fn2 <- x[['in_fn2']]
  if (verbose) logger::log_info("Loaded in_fn1 ({nrow(in_fn1)} SNPs) and in_fn2 ({nrow(in_fn2)} SNPs).")

  # Intersect SNP sets from both datasets
  snps <- intersect(in_fn1$rsid, in_fn2$rsid)
  snps <- unique(snps)
  if (verbose) logger::log_info("Number of SNPs in intersection for LD matrix: {length(snps)}")

  if (length(snps) < 2) {
    error_msg <- glue::glue("Not enough overlapping SNPs ({length(snps)}) to compute LD matrix. At least two required.")
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  }

  # Compute local LD matrix
  if (verbose) logger::log_info("Computing LD matrix for {length(snps)} SNPs using ...")
  if (verbose) logger::log_info("plink='{plink}'...")
  if (verbose) logger::log_info("bfile='{bfile}'...")
  ld_matrix <- tryCatch({
    ieugwasr::ld_matrix_local(
      variants = snps,
      bfile = bfile,
      plink_bin = plink,
      with_alleles = FALSE
    )

  }, error = function(e) {
    error_msg <- glue::glue("LD matrix calculation failed: {e$message}")
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  })
  if (verbose) logger::log_info("LD matrix computed: {nrow(ld_matrix)} x {ncol(ld_matrix)}.")

  # Filter and reorder association data to match LD matrix
  geni.dat <- assoc %>% dplyr::filter(marker %in% colnames(ld_matrix))
  geni.dat <- geni.dat[match(colnames(ld_matrix), geni.dat$marker), ]

  if (nrow(geni.dat) != ncol(ld_matrix)) {
    if (verbose) logger::log_info("Mismatch: geni.dat ({nrow(geni.dat)}) rows vs LD matrix ({ncol(ld_matrix)}) columns.")
  } else {
    if (verbose) logger::log_info("geni.dat order matches LD matrix columns.")
  }

  # Assemble output
  .data <- list(
    in_fn1 = in_fn1,
    in_fn2 = in_fn2,
    ld_matrix = ld_matrix,
    geni.dat = geni.dat,
    SNP = rs,
    pop = basename(bfile)
  )

  if (verbose) logger::log_info("Colocalization plotting data ready for SNP: {rs}")
  return(.data)
}


#' Prepare SuSiE Colocalization Results for Plotting
#'
#' This function processes colocalization results obtained from \code{COLOC.SUSIE} to prepare data
#' for visualization, focusing on a specific SNP. It extracts relevant subsets of the data and
#' organizes them into a format suitable for plotting tools such as \code{locuscomparer}.
#'
#' @param x List. Output from \code{COLOC.SUSIE}, containing significant SNPs, association data,
#'   and other colocalization results.
#' @param rs Character. Specific SNP ID (e.g., \code{"rs123"}) to focus on for plotting. Must follow
#'   the format \code{"rs"} followed by numbers.
#' @param verbose Logical. If \code{TRUE}, log detailed information about the data preparation process.
#'   Default: \code{TRUE}.
#'
#' @return A list containing the following elements:
#'   \describe{
#'     \item{\code{in_fn1}}{Data frame for the first trait (exposure), with columns for SNP ID (\code{rsid})
#'       and p-value (\code{pval}).}
#'     \item{\code{in_fn2}}{Data frame for the second trait (outcome), with columns for SNP ID (\code{rsid})
#'       and p-value (\code{pval}).}
#'     \item{\code{ld_matrix}}{Matrix of linkage disequilibrium (LD) correlations between SNPs.}
#'     \item{\code{geni.dat}}{Merged association data frame for plotting, containing marker information
#'       and Z-scores for both traits.}
#'     \item{\code{SNP}}{Character. The specified SNP ID for focus in visualization.}
#'     \item{\code{pop}}{Population information from the input data, if available.}
#'   }
#'
#' @details
#' This function validates the input SNP ID format and checks if it exists among significant SNPs identified
#' in the colocalization analysis. It then organizes the relevant data subsets for visualization, ensuring
#' compatibility with downstream plotting functions like \code{locuscomparer::locuscompare}.
#'
#' @examples
#' \dontrun{
#' # Assuming SUSIE.res is the output from COLOC.SUSIE
#' coloc.dat <- coloc_result_to_plot(
#'   x = SUSIE.res,
#'   rs = "rs1260326"
#' )
#'
#' # Use the prepared data for visualization
#' locuscomparer::locuscompare(
#'   in_fn1 = coloc.dat[["in_fn1"]],
#'   in_fn2 = coloc.dat[["in_fn2"]],
#'   title1 = "eQTL",
#'   title2 = "GWAS",
#'   genome = "hg19",
#'   snp = coloc.dat[["SNP"]],
#'   population = coloc.dat[["pop"]],
#'   combine = TRUE
#' )
#' }
#'
#' @seealso \code{\link{COLOC.SUSIE}}
#'
#' @author [Abao Xing]
#'
#' @importFrom glue glue
#' @importFrom dplyr filter
#' @importFrom logger log_info log_error
#'
#' @export
coloc_result_to_plot.susie <- function(x, rs, verbose = TRUE) {

  # Extract data
  snp.sig <- x[['snp.sig']]
  in_fn1 <- x[["in_fn1"]]
  in_fn2 <- x[["in_fn2"]]
  assoc <- x[["assoc"]]
  ld_matrix <- x[["ld_matrix"]]
  if (verbose) logger::log_info("Preparing SuSiE plotting data for SNP: {rs}")

  # Validate SNP format and presence
  if (!grepl("^rs[0-9]+$", rs)) {
    error_msg <- glue::glue("SNP '{rs}' does not follow the correct format (should be 'rs' followed by numbers, e.g., rs1234175).")
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  } else if (!(rs %in% snp.sig)) {
    error_msg <- glue::glue("SNP '{rs}' not found in significant SNPs (snp.sig). Please check your input.")
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  }

  # Intersect SNPs
  snps <- intersect(in_fn1$rsid, in_fn2$rsid)
  snps <- unique(snps)
  if (verbose) logger::log_info("Number of SNPs in intersection: {length(snps)}")

  # Assemble output
  data <- list(
    in_fn1 = in_fn1,
    in_fn2 = in_fn2,
    ld_matrix = ld_matrix,
    geni.dat = assoc,
    SNP = rs,
    pop = x[["pop"]]
  )

  if (verbose) logger::log_info("SuSiE plotting data ready for SNP: {rs}")
  return(data)
}


#' Prepare PWCOCO Colocalization Results for Plotting
#'
#' @description
#' Processes colocalization results from \code{COLOC.PWCoCo} for plotting, focusing on a specific SNP. This function
#' extracts relevant data for the specified SNP, computes an LD matrix, and prepares association data for visualization
#' using tools like \code{locuscomparer}.
#'
#' @param x List of class \code{pwcoco}. Output from \code{COLOC.PWCoCo}, containing colocalization and summary statistics data.
#' @param rs Character. Specific SNP ID (e.g., 'rs123') to focus on for plotting. Must follow the 'rs' + number format.
#' @param plink Character. Path to the PLINK executable for LD computation. Must be provided.
#' @param bfile Character. Path prefix for PLINK binary reference files (e.g., 1000G EUR). Must be provided.
#' @param verbose Logical. Whether to output verbose messages for debugging (default: TRUE).
#'
#' @return A list containing the following elements for plotting:
#'   \item{in_fn1}{Data frame for trait 1 with columns \code{rsid} (SNP ID) and \code{pval} (p-value).}
#'   \item{in_fn2}{Data frame for trait 2 with columns \code{rsid} (SNP ID) and \code{pval} (p-value).}
#'   \item{ld_matrix}{Matrix. SNP-by-SNP LD (linkage disequilibrium) matrix computed using PLINK.}
#'   \item{geni.dat}{Data frame. Merged association data for plotting with columns \code{marker}, \code{chr}, \code{pos}, \code{z_1}, and \code{z_2}.}
#'   \item{SNP}{Character. The specified SNP ID for focus.}
#'   \item{pop}{Character. The basename of the reference bfile used.}
#'
#' @details
#' This function processes the output of \code{COLOC.PWCoCo} to prepare data for colocalization plotting. It validates
#' the input SNP ID, loads summary statistics for the specified SNP, computes an LD matrix using PLINK, and merges
#' association data for the two traits. The resulting data structures are formatted for compatibility with plotting
#' tools such as \code{locuscomparer::locuscompare}. Ensure that PLINK and the reference files are accessible, and
#' the SNP ID exists in the colocalization results.
#'
#' @examples
#' \dontrun{
#' # Load required packages
#' library(GWASkitR)
#' library(locuscomparer)
#'
#' # Assuming PWCoCo.res is the output from COLOC.PWCoCo
#' # Prepare data for visualization
#' coloc.dat <- coloc_result_to_plot_pwcoco(
#'   data = PWCoCo.res,
#'   rs = "rs1260326",
#'   plink = "your/path/to/plink/plink.exe",
#'   bfile = "your/path/to/plink/1kg.v3/EUR",
#'   verbose = TRUE
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
#' @seealso \code{\link{COLOC.PWCoCo}}
#'
#' @author Abao Xing
#'
#' @importFrom glue glue
#' @importFrom data.table fread
#' @importFrom ieugwasr ld_matrix_local
#' @importFrom dplyr filter left_join
#' @importFrom logger log_info log_error log_info
#'
#' @export
coloc_result_to_plot.PWCoCo <- function(x, rs, plink = NULL, bfile = NULL, verbose = TRUE) {

  # Validate required arguments
  if (is.null(plink) || is.null(bfile)) {
    if (verbose) logger::log_error("Both 'plink' and 'bfile' parameters must be provided for LD matrix calculation.")
    stop("Both 'plink' and 'bfile' are required.")
  }

  # Extract data
  coloc.dat <- x[['coloc.dat']]
  sumstat1.list <- x[["sum_stats1"]]
  sumstat2.list <- x[["sum_stats2"]]
  if (verbose) logger::log_info("Preparing PWCOCO plotting data for SNP: {rs}")

  # Validate SNP format and presence
  if (!grepl("^rs[0-9]+$", rs)) {
    error_msg <- glue::glue("SNP '{rs}' does not follow the correct format (should be 'rs' followed by numbers, e.g., rs1234175).")
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  } else if (!(rs %in% coloc.dat$SNP1 | rs %in% coloc.dat$SNP2)) {
    error_msg <- glue::glue("SNP '{rs}' not found in coloc.dat$SNP1 or coloc.dat$SNP2. Please check your input.")
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  }

  # Load summary statistics for the specified SNP
  if (verbose) logger::log_info("Loading summary statistics for SNP {rs}.")
  AA <- tryCatch({
    data.table::fread(sumstat1.list[[rs]])
  }, error = function(e) {
    error_msg <- glue::glue("Error loading sum_stats1 for SNP {rs}: {e$message}")
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  })
  BB <- tryCatch({
    data.table::fread(sumstat2.list[[rs]])
  }, error = function(e) {
    error_msg <- glue::glue("Error loading sum_stats2 for SNP {rs}: {e$message}")
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  })

  # Create input data frames
  in_fn1 <- data.frame(rsid = AA$SNP, pval = AA$p)
  in_fn2 <- data.frame(rsid = BB$SNP, pval = BB$p)
  if (verbose) logger::log_info("Created in_fn1 ({nrow(in_fn1)} SNPs) and in_fn2 ({nrow(in_fn2)} SNPs).")

  # Intersect SNPs
  snps <- intersect(in_fn1$rsid, in_fn2$rsid)
  snps <- unique(snps)
  if (verbose) logger::log_info("Number of SNPs in intersection for LD matrix: {length(snps)}")

  if (length(snps) < 2) {
    error_msg <- glue::glue("Not enough overlapping SNPs ({length(snps)}) to compute LD matrix. At least two required.")
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  }

  # Compute local LD matrix
  if (verbose) logger::log_info("Computing LD matrix for {length(snps)} SNPs using ...")
  if (verbose) logger::log_info("plink='{plink}'...")
  if (verbose) logger::log_info("bfile='{bfile}'...")
  ld_matrix <- tryCatch({
    ieugwasr::ld_matrix_local(
      variants = snps,
      bfile = bfile,
      plink_bin = plink,
      with_alleles = FALSE
    )
  }, error = function(e) {
    error_msg <- glue::glue("LD matrix calculation failed: {e$message}")
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  })
  if (verbose) logger::log_info("LD matrix computed: {nrow(ld_matrix)} x {ncol(ld_matrix)}.")

  # Prepare association data
  AAA <- AA[AA$SNP %in% colnames(ld_matrix), c('SNP', 'Chr', 'bp', 'b', 'se')]
  BBB <- BB[BB$SNP %in% colnames(ld_matrix), c('SNP', 'b', 'se')]
  AAA$z_1 <- AAA$b / AAA$se
  BBB$z_2 <- BBB$b / BBB$se
  AAA$b <- AAA$se <- BBB$b <- BBB$se <- NULL
  assoc <- dplyr::left_join(AAA, BBB, by = 'SNP')
  colnames(assoc) <- c('marker', 'chr', 'pos', 'z_1', 'z_2')

  geni.dat <- assoc %>% dplyr::filter(marker %in% colnames(ld_matrix))
  geni.dat <- geni.dat[match(colnames(ld_matrix), geni.dat$marker), ]

  if (nrow(geni.dat) != ncol(ld_matrix)) {
    if (verbose) logger::log_info("Mismatch: geni.dat ({nrow(geni.dat)}) rows vs LD matrix ({ncol(ld_matrix)}) columns.")
  } else {
    if (verbose) logger::log_info("geni.dat order matches LD matrix columns.")
  }

  # Assemble output
  data <- list(
    in_fn1 = in_fn1,
    in_fn2 = in_fn2,
    ld_matrix = ld_matrix,
    geni.dat = geni.dat,
    SNP = rs,
    pop = basename(bfile)
  )

  if (verbose) logger::log_info("PWCOCO plotting data ready for SNP: {rs}")
  return(data)
}


#' Prepare HyPrColoc Colocalization Results for Visualization
#'
#' Processes colocalization results from \code{\link{COLOC.HyPrColoc}} and computes
#' a local LD matrix for intersecting SNPs using PLINK, preparing the results for
#' downstream regional visualization or integration with plotting packages.
#'
#' @param x A list, output from \code{\link{COLOC.HyPrColoc}}, containing association data,
#'   SNPs, betas, ses, and traits.
#' @param plink Character. Path to the PLINK executable, used for LD matrix computation.
#' @param bfile Character. Path prefix for PLINK binary reference files (e.g., 1000G EUR).
#' @param verbose Logical. If \code{TRUE}, log informational messages (default: \code{TRUE}).
#'
#' @return A list with elements:
#'   \describe{
#'     \item{\code{geni.dat}}{Merged association data, ready for visualization.}
#'     \item{\code{ld_matrix}}{SNP-by-SNP local LD matrix.}
#'     \item{\code{snps}}{Vector of included SNPs.}
#'     \item{\code{betas}}{Matrix of beta values for included SNPs.}
#'     \item{\code{ses}}{Matrix of standard errors for included SNPs.}
#'     \item{\code{traits}}{Vector of trait/dataset names.}
#'     \item{\code{pop}}{Population/biobank tag, extracted from \code{bfile} name.}
#'   }
#'
#' @details
#' This function is typically used after running \code{\link{COLOC.HyPrColoc}} to
#' harmonize association and LD data for regional visualization (e.g., with \pkg{geni.plots}).
#' Both \code{plink} and \code{bfile} are required for local LD matrix calculation.
#'
#' @section Workflow Example:
#' \enumerate{
#'   \item Run \code{COLOC.HyPrColoc} to obtain multi-trait colocalization results.
#'   \item Use this function to integrate LD matrix and regional association data for visualization.
#'   \item Visualize regions using \pkg{geni.plots} or \pkg{locuscomparer}.
#'   \item Optional: Perform sensitivity analysis and trait comparisons.
#' }
#'
#' @examples
#' \dontrun{
#' # Load required packages
#' library(GWASkitR)
#' library(locuscomparer)
#'
#' # Assuming PWCoCo.res is the output from COLOC.HyPrColoc
#' # Prepare data for visualization
#' coloc.dat <- coloc_result_to_plot(
#'   x = HyPrColoc.res,
#'   plink = "your/path/to/plink/plink.exe",
#'   bfile = "your/path/to/plink/1kg.v3/EUR"
#' )
#'
#' # Visualize with geni.plots (or other locus plotting tools)
#' geni.plots::fig_region_stack(
#'   data = coloc.dat[['geni.dat']],
#'   traits = coloc.dat[['traits']],
#'   corr = coloc.dat[['ld_matrix']],
#'   build = 37,
#'   title_center = TRUE
#' )
#'
#' # Sensitivity analysis (optional)
#' hyprcoloc::sensitivity.plot(
#'   effect.est = coloc.dat[['betas']],
#'   effect.se = coloc.dat[['ses']],
#'   trait.names = coloc.dat[['traits']],
#'   snp.id = coloc.dat[['snps']]
#' )
#' }
#'
#' @author Abao Xing (\email{albertxn7@@gmail.com})
#'
#' @seealso \code{\link{COLOC.HyPrColoc}}
#'
#' @importFrom ieugwasr ld_matrix_local
#' @importFrom dplyr filter
#' @importFrom logger log_info log_warn log_error
#' @export
coloc_result_to_plot.HyPrColoc <- function(x, plink = NULL, bfile = NULL, verbose = TRUE) {

  # Validate required arguments
  if (is.null(plink) || is.null(bfile)) {
    if (verbose) logger::log_error("Both 'plink' and 'bfile' parameters must be provided for LD matrix calculation.")
    stop("Both 'plink' and 'bfile' are required.")
  }

  # Extract data from input
  assoc <- x[["assoc"]]
  snps <- x[["snps"]]
  betas <- x[["betas"]]
  ses <- x[["ses"]]
  traits <- x[["traits"]]
  if (verbose) logger::log_info("Preparing HyPrColoc plotting data for {length(snps)} SNPs across {length(traits)} traits.")

  # Compute local LD matrix
  if (verbose) logger::log_info("Computing LD matrix for {length(snps)} SNPs using PLINK.")
  if (verbose) logger::log_info("PLINK executable: '{plink}'.")
  if (verbose) logger::log_info("Reference panel: '{bfile}'.")
  ld_matrix <- tryCatch({
    ieugwasr::ld_matrix_local(
      variants = snps,
      bfile = bfile,
      plink_bin = plink,
      with_alleles = FALSE
    )
  }, error = function(e) {
    error_msg <- glue::glue("LD matrix calculation failed: {e$message}")
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  })
  if (verbose) logger::log_info("LD matrix computed: {nrow(ld_matrix)} x {ncol(ld_matrix)} matrix.")

  # Filter and reorder association data to match LD matrix
  if (verbose) logger::log_info("Filtering association data to match LD matrix SNPs.")
  geni.dat <- assoc %>% dplyr::filter(marker %in% colnames(ld_matrix))
  geni.dat <- geni.dat[match(colnames(ld_matrix), geni.dat$marker), ]

  if (nrow(geni.dat) != ncol(ld_matrix)) {
    if (verbose) logger::log_warn("Mismatch detected: geni.dat has {nrow(geni.dat)} rows, while LD matrix has {ncol(ld_matrix)} columns.")
  } else {
    if (verbose) logger::log_info("geni.dat order successfully matched to LD matrix columns.")
  }

  # Assemble output
  if (verbose) logger::log_info("Assembling final output data structure.")
  data <- list(
    geni.dat = geni.dat,
    ld_matrix = ld_matrix,
    snps = snps,
    betas = betas,
    ses = ses,
    traits = traits,
    pop = basename(bfile)
  )

  if (verbose) logger::log_info("HyPrColoc plotting data prepared successfully for population '{data$pop}'.")
  return(data)
}
