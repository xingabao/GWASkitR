#' Perform Colocalization Analysis using Approximate Bayes Factor (ABF)
#'
#' @description
#' This function performs colocalization analysis between two traits using the
#' Approximate Bayes Factor (ABF) method as implemented in the \code{coloc} package.
#' It evaluates the evidence for shared causal variants between two datasets (e.g., eQTL and GWAS data).
#'
#' @param e.dat Data frame or file path for the first dataset (e.g., eQTL data).
#' @param o.dat Data frame or file path for the second dataset (e.g., GWAS data).
#' @param e.col.chr Character. Column name for chromosome in \code{e.dat}. Default: "CHR".
#' @param e.col.pos Character. Column name for position in \code{e.dat}. Default: "POS".
#' @param e.col.snp Character. Column name for SNP identifier in \code{e.dat}. Default: "SNP".
#' @param e.col.beta Character. Column name for effect size (beta) in \code{e.dat}. Default: "beta".
#' @param e.col.se Character. Column name for standard error in \code{e.dat}. Default: "se".
#' @param e.col.eaf Character. Column name for effect allele frequency in \code{e.dat}. Default: "eaf".
#' @param e.col.pval Character. Column name for p-value in \code{e.dat}. Default: "pval".
#' @param e.col.nsample Character or Numeric. Column name or numeric value for sample size in \code{e.dat}. Default: "N".
#' @param e.col.maf Character. Column name for minor allele frequency in \code{e.dat}. Default: \code{NA}.
#' @param e.type Character. Type of trait for \code{e.dat}, either "quant" (quantitative) or "cc" (case-control). Default: "quant".
#' @param e.prevalence Numeric. Disease prevalence for \code{e.dat} if \code{e.type = "cc"}. Default: \code{NA}.
#' @param o.col.snp Character. Column name for SNP identifier in \code{o.dat}. Default: "SNP".
#' @param o.col.beta Character. Column name for effect size (beta) in \code{o.dat}. Default: "beta".
#' @param o.col.se Character. Column name for standard error in \code{o.dat}. Default: "se".
#' @param o.col.eaf Character. Column name for effect allele frequency in \code{o.dat}. Default: "eaf".
#' @param o.col.pval Character. Column name for p-value in \code{o.dat}. Default: "pval".
#' @param o.col.nsample Character or Numeric. Column name or numeric value for sample size in \code{o.dat}. Default: "N".
#' @param o.col.maf Character. Column name for minor allele frequency in \code{o.dat}. Default: \code{NA}.
#' @param o.type Character. Type of trait for \code{o.dat}, either "quant" (quantitative) or "cc" (case-control). Default: "cc".
#' @param o.prevalence Numeric. Disease prevalence for \code{o.dat} if \code{o.type = "cc"}. Default: \code{NA}.
#' @param p1 Numeric. Prior probability for a SNP to be associated with trait 1. Default: 1e-04.
#' @param p2 Numeric. Prior probability for a SNP to be associated with trait 2. Default: 1e-04.
#' @param p12 Numeric. Prior probability for a SNP to be associated with both traits. Default: 1e-05.
#' @param verbose Logical. Whether to print detailed logging messages during execution. Default: TRUE.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{coloc.dat}{Full colocalization results from \code{coloc::coloc.abf}.}
#'   \item{summary}{Summary of posterior probabilities for each hypothesis (H0 to H4).}
#'   \item{results}{Detailed results per SNP, including posterior probabilities.}
#'   \item{priors}{Prior probabilities used in the analysis.}
#'   \item{assoc}{Merged dataset with z-scores for both traits.}
#'   \item{snp.sig}{The top SNP with the highest posterior probability for colocalization (H4).}
#'   \item{in_fn1}{Input data for trait 1 (SNP and p-value).}
#'   \item{in_fn2}{Input data for trait 2 (SNP and p-value).}
#' }
#'
#' @details
#' This function integrates data from two traits and performs colocalization analysis
#' using the \code{coloc::coloc.abf} function. It supports both quantitative and
#' case-control traits and allows customization of column names for input datasets.
#' Missing MAF values are imputed from effect allele frequencies if not provided.
#'
#' @examples
#' \dontrun{
#'
#' # Step 0. Load packages.
#' library(GWASkitR)
#' library(locuscomparer)
#'
#' # Step 1. Perform colocalization analysis.
#' ABF.res <- COLOC.ABF(
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
#'   e.prevalence = NA,
#'   o.col.snp = "SNP",
#'   o.col.beta = "beta",
#'   o.col.se = "se",
#'   o.col.eaf = "eaf",
#'   o.col.pval = "pval",
#'   o.col.nsample = 318014,
#'   o.col.maf = "MAF",
#'   o.type = "cc",
#'   o.prevalence = NA
#' )
#'
#' # Step 2. Organize the results for visualization.
#' coloc.dat <- coloc_result_to_plot(
#'   x = ABF.res,
#'   plink = "/your/path/to/plink/plink.exe",
#'   bfile = "/your/path/to/plink/plink/1kg.v3/EUR"
#' )
#'
#' # Step 3. Visualize the colocalization analysis results.
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
#' # Step 4. Sensitivity plot
#' coloc::sensitivity(
#'   obj = ABF.res[["coloc.dat"]],
#'   rule = "H4 > 0.8",
#'   plot.manhattans = TRUE
#' )
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom data.table fread
#' @importFrom dplyr select rename arrange left_join all_of
#' @importFrom coloc coloc.abf
#' @importFrom logger log_info log_warn log_error
#' @export
COLOC.ABF <- function(
    e.dat,
    o.dat,
    e.col.chr = 'CHR',
    e.col.pos = 'POS',
    e.col.snp = 'SNP',
    e.col.beta = 'beta',
    e.col.se = 'se',
    e.col.eaf = 'eaf',
    e.col.pval = 'pval',
    e.col.nsample = 'N',
    e.col.maf = NA,
    e.type = 'quant',
    e.prevalence = NA,
    o.col.snp = 'SNP',
    o.col.beta = 'beta',
    o.col.se = 'se',
    o.col.eaf = 'eaf',
    o.col.pval = 'pval',
    o.col.nsample = 'N',
    o.col.maf = NA,
    o.type = 'cc',
    o.prevalence = NA,
    p1 = 1e-04,
    p2 = 1e-04,
    p12 = 1e-05,
    verbose = TRUE
) {

  # Initialize logging
  start.time <- Sys.time()
  if (verbose) logger::log_info("Starting COLOC.ABF colocalization analysis...")

  # Validate trait types
  if (!e.type %in% c("quant", "cc")) {
    if (verbose) logger::log_error("e.type must be 'quant' or 'cc', got '{e.type}'")
    stop("e.type must be 'quant' or 'cc'")
  }
  if (!o.type %in% c("quant", "cc")) {
    if (verbose) logger::log_error("o.type must be 'quant' or 'cc', got '{o.type}'")
    stop("o.type must be 'quant' or 'cc'")
  }
  if (verbose) logger::log_info("Trait types validated: e.type='{e.type}', o.type='{o.type}'")

  # Read input data if provided as file paths
  if (is.character(e.dat) && file.exists(e.dat)) {
    if (verbose) logger::log_info("Reading eQTL data from file: {e.dat}")
    e.dat <- as.data.frame(data.table::fread(e.dat, showProgress = TRUE))
  }
  if (is.character(o.dat) && file.exists(o.dat)) {
    if (verbose) logger::log_info("Reading GWAS data from file: {o.dat}")
    o.dat <- as.data.frame(data.table::fread(o.dat, showProgress = TRUE))
  }

  # Assign fixed sample size if provided as numeric
  if (is.numeric(e.col.nsample)) {
    if (verbose) logger::log_info("Assigning fixed sample size for eQTL: {e.col.nsample}")
    e.dat$N <- e.col.nsample; e.col.nsample <- 'N'
  }
  if (is.numeric(o.col.nsample)) {
    if (verbose) logger::log_info("Assigning fixed sample size for GWAS: {o.col.nsample}")
    o.dat$N <- o.col.nsample; o.col.nsample <- 'N'
  }

  # Define required columns for each dataset
  e.cols <- c(e.col.snp, e.col.beta, e.col.se, e.col.eaf, e.col.pval, e.col.nsample, e.col.maf)
  o.cols <- c(o.col.snp, o.col.beta, o.col.se, o.col.eaf, o.col.pval, o.col.nsample, o.col.maf)
  e.cols <- e.cols[!is.na(e.cols)]
  o.cols <- o.cols[!is.na(o.cols)]

  # Validate that required columns exist
  missing_e <- setdiff(e.cols, names(e.dat))
  missing_o <- setdiff(o.cols, names(o.dat))
  if (length(missing_e) > 0) {
    if (verbose) logger::log_error("Missing required columns in eQTL data: {paste(missing_e, collapse=', ')}")
    stop("Missing columns in e.dat: ", paste(missing_e, collapse=', '))
  }
  if (length(missing_o) > 0) {
    if (verbose) logger::log_error("Missing required columns in GWAS data: {paste(missing_o, collapse=', ')}")
    stop("Missing columns in o.dat: ", paste(missing_o, collapse=', '))
  }

  if (verbose) logger::log_info("Required columns present in both datasets.")

  # SNP intersection and filtering
  snps <- intersect(e.dat[[e.col.snp]], o.dat[[o.col.snp]])
  snps <- unique(snps)
  if (verbose) logger::log_info("Number of SNPs in common: {length(snps)}")

  if (length(snps) == 0) {
    if (verbose) logger::log_error("No overlapping SNPs found between datasets!")
    stop("No overlapping SNPs.")
  }

  # Filter datasets to common SNPs and remove missing values
  e.dat <- e.dat[e.dat[[e.col.snp]] %in% snps, ] %>% na.omit()
  o.dat <- o.dat[o.dat[[o.col.snp]] %in% snps, ] %>% na.omit()
  e.dat <- e.dat[order(e.dat[[e.col.snp]]), ]
  o.dat <- o.dat[order(o.dat[[o.col.snp]]), ]
  if (verbose) logger::log_info("Filtered datasets to {nrow(e.dat)} (eQTL) and {nrow(o.dat)} (GWAS) SNPs after intersecting and omitting NAs.")

  # Extract input SNP and p-value tables for visualization or plotting
  e.fn <- dplyr::select(e.dat, rsid = dplyr::all_of(e.col.snp), pval = dplyr::all_of(e.col.pval))
  o.fn <- dplyr::select(o.dat, rsid = dplyr::all_of(o.col.snp), pval = dplyr::all_of(o.col.pval))
  if (verbose) logger::log_info("Extracted SNP and p-value tables for visualization.")

  # Compute or assign MAF if not provided
  if (!is.na(e.col.maf) && e.col.maf %in% names(e.dat)) {
    if (verbose) logger::log_info("Using provided MAF column for eQTL: {e.col.maf}")
    e.dat <- e.dat %>% dplyr::rename(MAF = all_of(e.col.maf))
  } else {
    if (verbose) logger::log_warn("MAF column missing for eQTL data; will impute MAF from EAF.")
    e.dat$MAF <- pmin(e.dat[[e.col.eaf]], 1 - e.dat[[e.col.eaf]])
  }
  if (!is.na(o.col.maf) && o.col.maf %in% names(o.dat)) {
    if (verbose) logger::log_info("Using provided MAF column for GWAS: {o.col.maf}")
    o.dat <- o.dat %>% dplyr::rename(MAF = all_of(o.col.maf))
  } else {
    if (verbose) logger::log_warn("MAF column missing for GWAS data; will impute MAF from EAF.")
    o.dat$MAF <- pmin(o.dat[[o.col.eaf]], 1 - o.dat[[o.col.eaf]])
  }

  # Replace NA MAF values with 0.5 as a default
  if (sum(is.na(e.dat[['MAF']])) > 0) {
    if (verbose) logger::log_warn("Some eQTL MAF values missing, replacing with 0.5")
    e.dat[['MAF']][is.na(e.dat[['MAF']])] <- 0.5
  }
  if (sum(is.na(o.dat[['MAF']])) > 0) {
    if (verbose) logger::log_warn("Some GWAS MAF values missing, replacing with 0.5")
    o.dat[['MAF']][is.na(o.dat[['MAF']])] <- 0.5
  }

  # Compute variance of beta (varbeta) as the square of standard error
  e.dat$varbeta <- e.dat[[e.col.se]]^2
  o.dat$varbeta <- o.dat[[o.col.se]]^2
  if (verbose) logger::log_info("Computed variance of beta (varbeta) for both datasets.")

  # Prepare lists for coloc.abf input
  e.list <- list(position = e.dat[[e.col.pos]])
  o.list <- list(position = e.dat[[e.col.pos]])
  elements <- c('pvalues', 'N', 'MAF', 'beta', 'varbeta', 'snp', 'type', 's')
  e.elements <- c(e.col.pval, e.col.nsample, 'MAF', e.col.beta, 'varbeta', e.col.snp, e.type, e.prevalence)
  o.elements <- c(o.col.pval, o.col.nsample, 'MAF', o.col.beta, 'varbeta', o.col.snp, o.type, o.prevalence)

  # Populate lists for coloc.abf input
  for (i in seq_along(elements)) {
    element <- elements[i]
    e.element <- e.elements[i]
    o.element <- o.elements[i]
    # eQTL
    if (!is.na(e.element)) {
      if (element %in% c("type", "s")) {
        e.list[[element]] <- e.element
      } else {
        e.list[[element]] <- e.dat[[e.element]]
      }
    }
    # GWAS
    if (!is.na(o.element)) {
      if (element %in% c("type", "s")) {
        o.list[[element]] <- o.element
      } else {
        o.list[[element]] <- o.dat[[o.element]]
      }
    }
  }

  if (verbose) logger::log_info("Prepared coloc.abf input lists.")

  # Run coloc.abf
  if (verbose) logger::log_info("Running coloc.abf with p1={p1}, p2={p2}, p12={p12}...")
  coloc.dat <- tryCatch({
    if (verbose) {
      suppressWarnings(
        coloc::coloc.abf(
          dataset1 = e.list,
          dataset2 = o.list,
          p1 = p1,
          p2 = p2,
          p12 = p12
        )
      )
    } else {
      suppressMessages(suppressWarnings(
        coloc::coloc.abf(
          dataset1 = e.list,
          dataset2 = o.list,
          p1 = p1,
          p2 = p2,
          p12 = p12
        )
      ))
    }

  }, error = function(e) {
    if (verbose) logger::log_error("coloc.abf failed: {e$message}")
    stop(e)
  })

  if (verbose) logger::log_info("coloc.abf run complete. Processing results...")

  # Sort results by posterior probability of H4 (colocalization)
  results <- coloc.dat[["results"]] %>% dplyr::arrange(-SNP.PP.H4)
  if (verbose) logger::log_info("Top SNP by H4: {results$snp[1]}, PP.H4: {results$SNP.PP.H4[1]}")

  # Compute z-scores for association data
  e.dat$z_1 <- e.dat[[e.col.beta]]/e.dat[[e.col.se]]
  o.dat$z_2 <- o.dat[[o.col.beta]]/o.dat[[o.col.se]]
  e.dat.c <- dplyr::select(e.dat, marker = dplyr::all_of(e.col.snp), chr = dplyr::all_of(e.col.chr), pos = dplyr::all_of(e.col.pos), z_1)
  o.dat.c <- dplyr::select(o.dat, marker = dplyr::all_of(o.col.snp), z_2)
  assoc <- dplyr::left_join(e.dat.c, o.dat.c, by = "marker")
  snp.sig <- results$snp[1]
  if (verbose) logger::log_info("Computed z-scores and merged association data for both traits.")

  # Summarize posterior probabilities
  summary <- data.frame(
    topSNP = results$snp[1],
    PP.H0.abf = coloc.dat[["summary"]][["PP.H0.abf"]],
    PP.H1.abf = coloc.dat[["summary"]][["PP.H1.abf"]],
    PP.H2.abf = coloc.dat[["summary"]][["PP.H2.abf"]],
    PP.H3.abf = coloc.dat[["summary"]][["PP.H3.abf"]],
    PP.H4.abf = coloc.dat[["summary"]][["PP.H4.abf"]]
  )

  if (verbose) logger::log_info("Posterior probability for colocalization (H4) is {summary$PP.H4.abf}")
  if (summary$PP.H4.abf > 0.8) {
    if (verbose) logger::log_info("Strong evidence for colocalization (H4 > 0.8).")
  } else if (summary$PP.H4.abf > 0.5) {
    if (verbose) logger::log_warn("Moderate evidence for colocalization (H4 > 0.5).")
  } else {
    if (verbose) logger::log_info("Weak or no evidence for colocalization (H4 <= 0.5).")
  }

  # Output assembly
  data <- list(
    coloc.dat = coloc.dat,
    summary = summary,
    results = coloc.dat[["results"]],
    priors = coloc.dat[["priors"]],
    assoc = assoc,
    snp.sig = snp.sig,
    in_fn1 = e.fn,
    in_fn2 = o.fn
  )

  # Log completion and execution time
  end.time <- Sys.time()
  execution.time <- end.time - start.time
  if (verbose) logger::log_info("Returning results.")
  if (verbose) logger::log_info("COLOC.ABF analysis completed in {round(execution.time, 2)} seconds.")
  if (verbose) logger::log_info(execution.time)

  class(data) <- c("abf", class(data))
  return(data)
}
