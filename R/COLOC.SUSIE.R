#' Perform Colocalization Analysis using COLOC with SuSiE
#'
#' This function performs colocalization analysis between two datasets (e.g., eQTL and GWAS data)
#' using the \code{coloc} package with SuSiE (Sum of Single Effects) for fine-mapping.
#' It processes input datasets, computes linkage disequilibrium (LD) matrices, and evaluates
#' colocalization evidence based on specified prior probabilities and thresholds.
#'
#' @param e.dat Data frame or character string. Input dataset for exposure (e.g., eQTL data) or path to the data file.
#' @param o.dat Data frame or character string. Input dataset for outcome (e.g., GWAS data) or path to the data file.
#' @param e.col.chr Character. Column name for chromosome in exposure data. Default: \code{"CHR"}.
#' @param e.col.pos Character. Column name for position in exposure data. Default: \code{"POS"}.
#' @param e.col.snp Character. Column name for SNP ID in exposure data. Default: \code{"SNP"}.
#' @param e.col.beta Character. Column name for beta effect in exposure data. Default: \code{"beta"}.
#' @param e.col.se Character. Column name for standard error in exposure data. Default: \code{"se"}.
#' @param e.col.eaf Character. Column name for effect allele frequency in exposure data. Default: \code{"eaf"}.
#' @param e.col.pval Character. Column name for p-value in exposure data. Default: \code{"pval"}.
#' @param e.col.ncase Character or numeric. Column name or value for number of cases in exposure data. Default: \code{NA}.
#' @param e.col.nsample Character or numeric. Column name or value for sample size in exposure data. Default: \code{NA}.
#' @param e.col.maf Character. Column name for minor allele frequency in exposure data. Default: \code{NA}.
#' @param e.type Character. Type of exposure data, either \code{"quant"} for quantitative or \code{"cc"} for case-control. Default: \code{"quant"}.
#' @param o.col.snp Character. Column name for SNP ID in outcome data. Default: \code{"SNP"}.
#' @param o.col.beta Character. Column name for beta effect in outcome data. Default: \code{"beta"}.
#' @param o.col.se Character. Column name for standard error in outcome data. Default: \code{"se"}.
#' @param o.col.eaf Character. Column name for effect allele frequency in outcome data. Default: \code{"eaf"}.
#' @param o.col.pval Character. Column name for p-value in outcome data. Default: \code{"pval"}.
#' @param o.col.ncase Character or numeric. Column name or value for number of cases in outcome data. Default: \code{NA}.
#' @param o.col.nsample Character or numeric. Column name or value for sample size in outcome data. Default: \code{NA}.
#' @param o.col.maf Character. Column name for minor allele frequency in outcome data. Default: \code{NA}.
#' @param o.type Character. Type of outcome data, either \code{"quant"} for quantitative or \code{"cc"} for case-control. Default: \code{"cc"}.
#' @param o.prevalence Numeric. Prevalence of outcome for case-control data. Default: \code{NA}.
#' @param p1 Numeric. Prior probability for association in exposure dataset. Default: \code{1e-04}.
#' @param p2 Numeric. Prior probability for association in outcome dataset. Default: \code{1e-04}.
#' @param p12 Numeric. Prior probability for colocalization. Default: \code{1e-05}.
#' @param threshold Numeric. Posterior probability threshold for significant colocalization. Default: \code{0.8}.
#' @param plink Character. Path to PLINK executable for LD computation. Default: \code{NULL}.
#' @param bfile Character. Path to PLINK binary file for LD reference. Default: \code{NULL}.
#' @param verbose Logical. If \code{TRUE}, log detailed information about the analysis process. Default: \code{TRUE}.
#'
#' @return A list containing the following elements:
#'   \describe{
#'     \item{\code{coloc.dat}}{Full colocalization results from \code{coloc.susie}.}
#'     \item{\code{dataset1}}{List of processed exposure data for colocalization.}
#'     \item{\code{dataset2}}{List of processed outcome data for colocalization.}
#'     \item{\code{summary}}{Summary of colocalization results.}
#'     \item{\code{results}}{Data frame of detailed colocalization results sorted by posterior probability.}
#'     \item{\code{priors}}{Prior probabilities used in the analysis.}
#'     \item{\code{assoc}}{Data frame with marker information and Z-scores for both datasets.}
#'     \item{\code{snp.sig}}{Vector of significant SNPs based on the colocalization threshold.}
#'     \item{\code{in_fn1}}{Subset of exposure data with rsID and p-values.}
#'     \item{\code{in_fn2}}{Subset of outcome data with rsID and p-values.}
#'     \item{\code{ld_matrix}}{Matrix of LD correlations between SNPs.}
#'   }
#'
#' @details
#' This function integrates the \code{coloc} package's SuSiE framework to perform fine-mapping and colocalization
#' analysis between two traits, such as eQTL and GWAS data. It requires LD information, which is computed using
#' PLINK if provided. The function handles data preprocessing, including filtering for common SNPs and preparing
#' input for the \code{coloc.susie} function. Results are returned in a comprehensive list format for downstream
#' visualization and interpretation.
#'
#' @examples
#' \dontrun{
#' # Step 0. Load packages.
#' library(GWASkitR)
#' library(locuscomparer)
#'
#' # Step 1: Perform analysis
#' SUSIE.res <- COLOC.SUSIE(
#'   e.dat = GWASkitR::pQTLdat,
#'   o.dat = GWASkitR::sumstats,
#'   e.col.chr = "CHR",
#'   e.col.pos = "POS",
#'   e.col.snp = "SNP",
#'   e.col.beta = "beta",
#'   e.col.se = "se",
#'   e.col.eaf = "eaf",
#'   e.col.pval = "pval",
#'   e.col.nsample = 35373,
#'   e.col.maf = "MAF",
#'   e.type = "quant",
#'   o.col.snp = "SNP",
#'   o.col.beta = "beta",
#'   o.col.se = "se",
#'   o.col.eaf = "eaf",
#'   o.col.pval = "pval",
#'   o.col.ncase = 33043,
#'   o.col.nsample = 318014,
#'   o.col.maf = "MAF",
#'   o.type = "cc",
#'   plink = "/your/path/to/plink/plink.exe",
#'   bfile = "/your/path/to/plink/1kg.v3/EUR"
#' )
#'
#' # Step 2: Data preparation for visualization
#' coloc.dat <- coloc_result_to_plot(
#'   x = SUSIE.res,
#'   rs = "rs1260326"
#' )
#'
#' # Step 3: Visualization with locuscompare
#' title1 <- "eQTL"
#' title2 <- "GWAS"
#' locuscomparer::locuscompare(
#'   in_fn1 = coloc.dat[["in_fn1"]],
#'   in_fn2 = coloc.dat[["in_fn2"]],
#'   title1 = title1,
#'   title2 = title2,
#'   genome = "hg19",
#'   snp = coloc.dat[["SNP"]],
#'   population = coloc.dat[["pop"]],
#'   combine = TRUE
#' )
#'
#' # Step 4: Sensitivity plot
#' coloc::sensitivity(
#'   obj = SUSIE.res[["coloc.dat"]],
#'   rule = "H4 > 0.8",
#'   row = 1,
#'   dataset1 = SUSIE.res[["dataset1"]],
#'   dataset2 = SUSIE.res[["dataset2"]],
#'   plot.manhattans = TRUE
#' )
#' }
#'
#' @author [Abao Xing]
#'
#' @references
#' [1] Liu, B., Gloudemans, M.J., Rao, A.S. et al. Abundant associations with gene
#' expression complicate GWAS follow-up. Nat Genet 51, 768–769 (2019).
#'  https://doi.org/10.1038/s41588-019-0404-0
#'
#' @importFrom magrittr %>%
#' @importFrom logger log_info log_warn log_error
#' @importFrom data.table fread
#' @importFrom dplyr select rename filter arrange left_join
#' @importFrom ieugwasr ld_matrix_local
#' @importFrom coloc coloc.susie
#'
#' @export
COLOC.SUSIE <- function(
    e.dat,
    o.dat,
    e.col.chr = 'CHR',
    e.col.pos = 'POS',
    e.col.snp = 'SNP',
    e.col.beta = 'beta',
    e.col.se = 'se',
    e.col.eaf = 'eaf',
    e.col.pval = 'pval',
    e.col.ncase = NA,
    e.col.nsample = NA,
    e.col.maf = NA,
    e.type = 'quant',
    o.col.snp = 'SNP',
    o.col.beta = 'beta',
    o.col.se = 'se',
    o.col.eaf = 'eaf',
    o.col.pval = 'pval',
    o.col.ncase = NA,
    o.col.nsample = NA,
    o.col.maf = NA,
    o.type = 'cc',
    o.prevalence = NA,
    p1 = 1e-04,
    p2 = 1e-04,
    p12 = 1e-05,
    threshold = 0.8,
    plink = NULL,
    bfile = NULL,
    verbose = TRUE
) {

  # Initialize logging
  start.time <- Sys.time()
  if (verbose) logger::log_info("Starting colocalization analysis with COLOC.SuSiE.")

  # Log input parameters for debugging
  if (verbose) logger::log_info("Input parameters: e.type={e.type}, o.type={o.type}, p1={p1}, p2={p2}, p12={p12}, threshold={threshold}")

  # Validate input data types
  if (verbose) logger::log_info("Validating input data types for exposure and outcome.")
  stopifnot(e.type %in% c("quant", "cc"))
  stopifnot(o.type %in% c("quant", "cc"))
  if (verbose) logger::log_info("Data types validated: exposure={e.type}, outcome={o.type}.")

  # Load data if provided as file paths
  if (verbose) logger::log_info("Loading input data if provided as file paths.")
  if (is.character(e.dat) && file.exists(e.dat)) {
    e.dat <- tryCatch({
      as.data.frame(data.table::fread(e.dat, showProgress = FALSE))
    }, error = function(e) {
      error_msg <- paste("Error loading exposure data file:", e$message)
      if (verbose) logger::log_error(error_msg)
      stop(error_msg)
    })
    if (verbose) logger::log_info("Exposure data loaded from file: {e.dat}")
  }
  if (is.character(o.dat) && file.exists(o.dat)) {
    o.dat <- tryCatch({
      as.data.frame(data.table::fread(o.dat, showProgress = FALSE))
    }, error = function(e) {
      error_msg <- paste("Error loading outcome data file:", e$message)
      if (verbose) logger::log_error(error_msg)
      stop(error_msg)
    })
    if (verbose) logger::log_info("Outcome data loaded from file: {o.dat}")
  }

  # Assign sample size and case numbers if provided as numeric values
  if (verbose) logger::log_info("Assigning sample size and case numbers if provided as numeric.")
  if (is.numeric(e.col.nsample)) { e.dat$N <- e.col.nsample; e.col.nsample <- 'N' }
  if (is.numeric(o.col.nsample)) { o.dat$N <- o.col.nsample; o.col.nsample <- 'N' }
  if (is.numeric(e.col.ncase)) { e.dat$ncase <- e.col.ncase; e.col.ncase <- 'ncase' }
  if (is.numeric(o.col.ncase)) { o.dat$ncase <- o.col.ncase; o.col.ncase <- 'ncase' }
  if (verbose) logger::log_info("Sample size and case assignments completed.")

  # Validate required columns in datasets
  if (verbose) logger::log_info("Validating required columns in datasets.")
  e.cols <- c(e.col.snp, e.col.chr, e.col.pos, e.col.beta, e.col.se, e.col.eaf, e.col.pval, e.col.maf, e.col.nsample, e.col.ncase)
  o.cols <- c(o.col.snp, o.col.beta, o.col.se, o.col.eaf, o.col.pval, o.col.maf, o.col.nsample, o.col.ncase)
  e.cols <- e.cols[!is.na(e.cols)]
  o.cols <- o.cols[!is.na(o.cols)]
  if (!all(e.cols %in% names(e.dat))) {
    missing_cols <- e.cols[!e.cols %in% names(e.dat)]
    error_msg <- paste("Missing required columns in exposure data:", paste(missing_cols, collapse = ", "))
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  }
  if (!all(o.cols %in% names(o.dat))) {
    missing_cols <- o.cols[!o.cols %in% names(o.dat)]
    error_msg <- paste("Missing required columns in outcome data:", paste(missing_cols, collapse = ", "))
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  }
  if (verbose) logger::log_info("All required columns validated in both datasets.")

  # Calculate sample sizes and case proportions
  if (verbose) logger::log_info("Calculating sample sizes and case proportions.")
  e.sample <- as.integer(mean(e.dat[[e.col.nsample]], na.rm = TRUE))
  o.sample <- as.integer(mean(o.dat[[o.col.nsample]], na.rm = TRUE))
  if (e.type == 'cc') { e.ncase <- as.integer(mean(e.dat[[e.col.ncase]], na.rm = TRUE)) }
  if (o.type == 'cc') { o.ncase <- as.integer(mean(o.dat[[o.col.ncase]], na.rm = TRUE)) }
  if (verbose) logger::log_info("Sample sizes: exposure={e.sample}, outcome={o.sample}")
  if (e.type == 'cc') logger::log_info("Exposure cases: {e.ncase}")
  if (o.type == 'cc') logger::log_info("Outcome cases: {o.ncase}")

  # Identify common SNPs and filter datasets
  if (verbose) logger::log_info("Identifying common SNPs between datasets.")
  snps <- intersect(e.dat[[e.col.snp]], o.dat[[o.col.snp]])
  snps <- unique(snps)
  if (verbose) logger::log_info("Found {length(snps)} common SNPs.")
  if (length(snps) == 0) {
    error_msg <- "No common SNPs found between exposure and outcome datasets."
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  }
  e.dat <- e.dat[e.dat[[e.col.snp]] %in% snps, ] %>% na.omit() %>% as.data.frame()
  o.dat <- o.dat[o.dat[[o.col.snp]] %in% snps, ] %>% na.omit() %>% as.data.frame()
  e.dat <- e.dat[order(e.dat[[e.col.snp]]), ]
  o.dat <- o.dat[order(o.dat[[o.col.snp]]), ]
  if (verbose) logger::log_info("Filtered datasets to common SNPs; exposure rows={nrow(e.dat)}, outcome rows={nrow(o.dat)}")

  # Prepare subsets for p-values
  if (verbose) logger::log_info("Preparing subsets with rsID and p-values.")
  e.fn <- dplyr::select(e.dat, rsid = dplyr::all_of(e.col.snp), pval = dplyr::all_of(e.col.pval))
  o.fn <- dplyr::select(o.dat, rsid = dplyr::all_of(o.col.snp), pval = dplyr::all_of(o.col.pval))
  if (verbose) logger::log_info("Subsets for p-values created.")

  # Compute or assign MAF if not provided
  if (verbose) logger::log_info("Computing or assigning MAF for datasets.")
  if (e.col.maf %in% names(e.dat)) {
    e.dat <- e.dat %>% dplyr::rename(MAF = dplyr::all_of(e.col.maf))
  } else {
    e.dat$MAF <- pmin(e.dat[[e.col.eaf]], 1 - e.dat[[e.col.eaf]])
  }
  if (o.col.maf %in% names(o.dat)) {
    o.dat <- o.dat %>% dplyr::rename(MAF = dplyr::all_of(o.col.maf))
  } else {
    o.dat$MAF <- pmin(o.dat[[o.col.eaf]], 1 - o.dat[[o.col.eaf]])
  }
  if (sum(is.na(e.dat[['MAF']])) > 0) { e.dat[['MAF']][is.na(e.dat[['MAF']])] <- 0.5 }
  if (sum(is.na(o.dat[['MAF']])) > 0) { o.dat[['MAF']][is.na(o.dat[['MAF']])] <- 0.5 }
  if (verbose) logger::log_info("MAF values assigned or computed.")

  # Compute LD matrix using PLINK
  if (verbose) logger::log_info("Computing LD matrix for common SNPs.")
  ld_matrix <- tryCatch({
    ieugwasr::ld_matrix_local(
      variants = snps,
      bfile = bfile,
      plink_bin = plink,
      with_alleles = FALSE
    )
  }, error = function(e) {
    error_msg <- paste("Error computing LD matrix:", e$message)
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  })
  if (verbose) logger::log_info("LD matrix computed with dimensions: {dim(ld_matrix)[1]} x {dim(ld_matrix)[2]}")

  # Reorder datasets based on LD matrix row names
  if (verbose) logger::log_info("Reordering datasets based on LD matrix row names.")
  rownames(e.dat) <- e.dat[[e.col.snp]]
  e.dat <- e.dat[rownames(ld_matrix), ]
  rownames(o.dat) <- o.dat[[o.col.snp]]
  o.dat <- o.dat[rownames(ld_matrix), ]
  if (verbose) logger::log_info("Datasets reordered to match LD matrix.")

  # Prepare data lists for coloc.susie
  if (verbose) logger::log_info("Preparing data lists for coloc.susie analysis.")
  e.list <- list(
    snp = e.dat[[e.col.snp]],
    position = e.dat[[e.col.pos]],
    beta = e.dat[[e.col.beta]],
    varbeta = e.dat[[e.col.se]]^2,
    pvalues = e.dat[[e.col.pval]],
    MAF = e.dat[['MAF']],
    type = e.type,
    N = e.sample
  )
  if (e.type == "cc") { e.list$s <- as.numeric(e.ncase / e.sample) }

  o.list <- list(
    snp = o.dat[[o.col.snp]],
    position = e.dat[[e.col.pos]],  # Note: using exposure position as in original code
    beta = o.dat[[o.col.beta]],
    varbeta = o.dat[[o.col.se]]^2,
    pvalues = o.dat[[o.col.pval]],
    MAF = o.dat[['MAF']],
    type = o.type,
    N = o.sample
  )
  if (o.type == "cc") { o.list$s <- as.numeric(o.ncase / o.sample) }

  e.list$LD <- ld_matrix
  o.list$LD <- ld_matrix
  if (verbose) logger::log_info("Data lists prepared for colocalization.")

  # Perform colocalization analysis with coloc.susie
  if (verbose) logger::log_info("Running coloc.susie for colocalization analysis.")
  coloc.dat <- tryCatch({
    suppressMessages({ coloc::coloc.susie(e.list, o.list, p1 = p1, p2 = p2, p12 = p12) })
  }, error = function(e) {
    error_msg <- paste("Error in coloc.susie analysis:", e$message)
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  })
  if (verbose) logger::log_info("Colocalization analysis completed.")

  # Process results and identify significant SNPs
  if (verbose) logger::log_info("Processing colocalization results.")
  results <- coloc.dat$results %>% dplyr::arrange(-SNP.PP.H4.abf)
  snp.sig <- results %>% dplyr::filter(SNP.PP.H4.abf >= threshold)
  snp.sig <- snp.sig$snp
  if (verbose) logger::log_info("Identified {length(snp.sig)} significant SNPs with PP.H4 >= {threshold}.")

  # Compute Z-scores and build association data frame
  if (verbose) logger::log_info("Computing Z-scores and building association data frame.")
  e.dat$z_1 <- e.dat[[e.col.beta]] / e.dat[[e.col.se]]
  e.dat.c <- dplyr::select(e.dat, marker = dplyr::all_of(e.col.snp), chr = dplyr::all_of(e.col.chr), pos = dplyr::all_of(e.col.pos), z_1)
  o.dat$z_2 <- o.dat[[o.col.beta]] / o.dat[[o.col.se]]
  o.dat.c <- dplyr::select(o.dat, marker = dplyr::all_of(o.col.snp), z_2)
  assoc <- dplyr::left_join(e.dat.c, o.dat.c, by = "marker")
  if (verbose) logger::log_info("Association data frame created with Z-scores.")

  # Compile output data structure
  if (verbose) logger::log_info("Compiling output results.")

  data <- list()
  data[["coloc.dat"]] <- coloc.dat
  data[["dataset1"]] <- e.list
  data[["dataset2"]] <- o.list
  data[["summary"]] <- coloc.dat$summary
  data[["results"]] <- results
  data[["priors"]] <- coloc.dat[["priors"]]
  data[["assoc"]] <- assoc
  data[["snp.sig"]] <- snp.sig
  data[["in_fn1"]] <- e.fn
  data[["in_fn2"]] <- o.fn
  data[["ld_matrix"]] <- ld_matrix
  data[["pop"]] <- basename(bfile)

  # Log completion and execution time
  end.time <- Sys.time()
  execution.time <- end.time - start.time
  if (verbose) logger::log_info("Returning results.")
  if (verbose) logger::log_info("COLOC.SuSiE analysis completed in {round(execution.time, 2)} seconds.")
  if (verbose) logger::log_info(execution.time)

  class(data) <- c('susie', class(data))
  return(data)
}
