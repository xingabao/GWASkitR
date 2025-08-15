# remotes::install_github("alexploner/cfdr.pleio", build_vignettes = TRUE)
#' @export
cond.conjFDR <- function(
    sumstats.1,
    sumstats.2,
    col.snp.1 = 'SNP',
    col.beta.1 = 'beta',
    col.pval.1 = 'pval',
    col.snp.2 = 'SNP',
    col.beta.2 = 'beta',
    col.pval.2 = 'pval',
    n_iter = 50,
    seed = 1234175,
    trait1_name = "trait1_name",
    trait2_name = "trait2_name",
    refdata_location,
    temp_dir = tempdir(),
    save.path = NULL,
    verbose = TRUE
) {

  start.time = Sys.time()
  # exe_py <- system.file("bin/rc/cond_conjFDR.Rc", package = "MRsoft")
  # compiler::loadcmp(exe_py, env = environment())
  # a <- suppressWarnings(tryCatch(if (!"cfdr.pleio" %in% installed.packages()[, "Package"]) {
  #   cli::cli_alert_info("正在安装所需的 cfdr.pleio 包...")
  #   QTLMRget:::install_soft("cfdr.pleio")
  #   message("")
  # }, error = function(e) "error"))
  # cli::cli_alert_info("正在读取数据中，请稍后...")
  # verbose = TRUE
  # sumstats.1 = "C:/Users/Administrator/Downloads/27089181-GCST003769-EFO_0007006-build37.f.tsv.gz"
  # sumstats.2 = "C:/Users/Administrator/Downloads/ukb-b-12493.vcf.gz"
  # sumstats.1 = "E:/GWAS/GWAS/SNP_gwas_mc_merge_nogc.tbl.uniq.gz"
  # sumstats.2 = "E:/GWAS/GWAS/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz"
  # col.snp.1 = 'SNP'
  # col.beta.1 = 'beta'
  # col.pval.1 = 'pval'
  # #
  # col.snp.2 = 'SNP'
  # col.beta.2 = 'beta'
  # col.pval.2 = 'pval'
  #
  # col.snp.1 = 'variant_id'
  # col.beta.1 = 'beta'
  # col.pval.1 = 'p_value'
  #
  # col.snp.2 = 'ID'
  # col.beta.2 = 'ES'
  # col.pval.2 = 'PVAL'

  # col.snp.1 = 'SNP'
  # col.beta.1 = 'b'
  # col.pval.1 = 'p'
  #
  # col.snp.2 = 'MarkerName'
  # col.beta.2 = 'b'
  # col.pval.2 = 'p'


  # n_iter = 50
  # seed = 154226

  # trait1_name = "BMI"
  # trait2_name = "Height"

  # refdata_location = "E:/GWAS/REF_DIR"
  # temp_dir = tempdir()
  # save.path = NULL

  sumstats.o.1 <- sumstats.1
  sumstats.o.2 <- sumstats.2

  # If input is a file path, read the file
  if (is.character(sumstats.1) && file.exists(sumstats.1)) {
    if (endsWith(sumstats.1, 'vcf') || endsWith(sumstats.1, 'vcf.gz')) {
      sumstats.1 <- GWASkitR::read.vcf(sumstats.1, nrows = 1000, verbose = verbose)
    } else {
      sumstats.1 <- as.data.frame(data.table::fread(sumstats.1, showProgress = FALSE, nrows = 1000))
    }
  } else {
    if (verbose) logger::log_error("The file '{sumstats.o.2}' does not exist.")
    stop("Input summary statistics file does not exist: ", sumstats.o.2)
  }

  if (is.character(sumstats.2) && file.exists(sumstats.2)) {
    if (endsWith(sumstats.2, 'vcf') || endsWith(sumstats.2, 'vcf.gz')) {
      sumstats.2 <- GWASkitR::read.vcf(sumstats.2, nrows = 1000, verbose = verbose)
    } else {
      sumstats.2 <- as.data.frame(data.table::fread(sumstats.2, showProgress = FALSE, nrows = 1000))
    }
  } else {
    if (verbose) logger::log_error("The file '{sumstats.o.2}' does not exist.")
    stop("Input summary statistics file does not exist: ", sumstats.o.2)
  }

  # Collect the required columns
  .cols.1 <- c(col.snp.1, col.beta.1, col.pval.1)
  .cols.2 <- c(col.snp.2, col.beta.2, col.pval.2)

  # Remove NULLs from .cols.1/.cols.2 to avoid checking for missing NULL columns
  .cols.1 <- .cols.1[!is.null(.cols.1)]
  .cols.2 <- .cols.2[!is.null(.cols.2)]

  # Validate that required columns exist in the data
  missing.cols.1 <- setdiff(.cols.1, names(sumstats.1))
  if (length(missing.cols.1) > 0) {
    if (verbose) logger::log_warn("These are the column names in the input data: {paste(names(sumstats.1), collapse = ', ')}")
    if (verbose) logger::log_info("Preview of the first few rows of input data:")
    print(utils::head(sumstats.1))
    if (verbose) logger::log_error("Missing required columns: {paste(missing.cols.1, collapse = ', ')}")
    stop("Missing columns: ", paste(missing.cols.1, collapse = ', '))
  } else if (length(.cols.1) == 0) {
    if (verbose) logger::log_warn("No columns exist in the data frame, please check.")
    if (verbose) logger::log_warn("These are the column names in the input data: {paste(names(sumstats.1), collapse = ', ')}")
    if (verbose) logger::log_info("Preview of the first few rows of input data:")
    print(utils::head(sumstats.1))
  }
  missing.cols.2 <- setdiff(.cols.2, names(sumstats.2))
  if (length(missing.cols.2) > 0) {
    if (verbose) logger::log_warn("These are the column names in the input data: {paste(names(sumstats.2), collapse = ', ')}")
    if (verbose) logger::log_info("Preview of the first few rows of input data:")
    print(utils::head(sumstats.2))
    if (verbose) logger::log_error("Missing required columns: {paste(missing.cols.2, collapse = ', ')}")
    stop("Missing columns: ", paste(missing.cols.2, collapse = ', '))
  } else if (length(.cols.2) == 0) {
    if (verbose) logger::log_warn("No columns exist in the data frame, please check.")
    if (verbose) logger::log_warn("These are the column names in the input data: {paste(names(sumstats.2), collapse = ', ')}")
    if (verbose) logger::log_info("Preview of the first few rows of input data:")
    print(utils::head(sumstats.2))
  }

  # If input is a file path, read the file
  sumstats.1 <- sumstats.o.1
  if (is.character(sumstats.1) && file.exists(sumstats.1)) {
    if (verbose) logger::log_info("Detected local file: '{sumstats.1}'. Reading as table.")
    if (endsWith(sumstats.1, 'vcf') || endsWith(sumstats.1, 'vcf.gz')) {
      sumstats.1 <- GWASkitR::read.vcf(sumstats.1, verbose = verbose)
    } else {
      sumstats.1 <- as.data.frame(data.table::fread(sumstats.1, showProgress = verbose))
    }
    if (verbose) logger::log_info("Local file loaded successfully.")
  } else {
    if (verbose) logger::log_info("Input provided as data.frame or data.table. Proceeding.")
  }
  sumstats.2 <- sumstats.o.2
  if (is.character(sumstats.2) && file.exists(sumstats.2)) {
    if (verbose) logger::log_info("Detected local file: '{sumstats.2}'. Reading as table.")
    if (endsWith(sumstats.2, 'vcf') || endsWith(sumstats.2, 'vcf.gz')) {
      sumstats.2 <- GWASkitR::read.vcf(sumstats.2, verbose = verbose)
    } else {
      sumstats.2 <- as.data.frame(data.table::fread(sumstats.2, showProgress = verbose))
    }
    if (verbose) logger::log_info("Local file loaded successfully.")
  } else {
    if (verbose) logger::log_info("Input provided as data.frame or data.table. Proceeding.")
  }

  if (verbose) logger::log_info("Selecting required columns: {paste(.cols.1, collapse=', ')} from sumstats1")
  sumstats.dt.1 <- sumstats.1 %>% dplyr::select(dplyr::all_of(.cols.1))
  if (verbose) logger::log_info("Selecting required columns: {paste(.cols.2, collapse=', ')} from sumstats2")
  sumstats.dt.2 <- sumstats.2 %>% dplyr::select(dplyr::all_of(.cols.2))

  trait1 <- sumstats.dt.1
  trait2 <- sumstats.dt.2
  colnames(trait1) <- c("SNP", "BETA", "PVAL")
  colnames(trait2) <- c("SNP", "BETA", "PVAL")

  dat <- cfdr.pleio::cfdr_pleio$new()

  dat$init_data(
    trait1 = trait1, trait2 = trait2,
    trait_names = c(trait1_name, trait2_name),
    refdat = cfdr.pleio::refdata_location(refdata_location),
    local_refdat_path = temp_dir,
    verbose = verbose
  )

  dat$initialize_pruning_index(n_iter = n_iter, seed = seed, verbose = verbose)

  # 正在计算FDR，请稍后...
  suppressWarnings(dat$calculate_cond_fdr(fdr_trait = 1, verbose = TRUE))
  suppressWarnings(dat$calculate_cond_fdr(fdr_trait = 2, verbose = TRUE))

  data_res <- dat$get_trait_results()

  data_res$PVAL1 <- p_value <- 10^(-data_res$LOG10PVAL1)
  data_res$PVAL2 <- p_value <- 10^(-data_res$LOG10PVAL2)

  data_res <- dplyr::arrange(data_res, conj_fdr)

  # Optionally, save the formatted data to file
  if (!is.null(save.path)) {
    save.sumstats(data_res, save.path)
  } else {
    if (verbose) logger::log_info("No save.path provided. Data will not be written to disk.")
  }

  # Log completion and execution time
  if (verbose) execution.time <- Sys.time() - start.time
  if (verbose) logger::log_info("GWAS summary statistics formatting completed in {round(execution.time, 2)} minutes")

  gc()

  return(data_res %>% dplyr::as_tibble())
}
