#' Perform Genomic Coordinate Liftover Between Builds
#'
#' This function converts genomic coordinates from one genome build to another (e.g., hg38 to hg19)
#' using a chain file for liftover. It supports conversion of summary statistics data and can return
#' results as a data frame or GRanges object.
#'
#' @param sumstats Data frame, input summary statistics data with chromosomal positions.
#' @param build.from Character, source genome build (e.g., 'hg38'; default: 'hg38').
#' @param build.to Character, target genome build (e.g., 'hg19'; default: 'hg19').
#' @param chain.file Character, path to a custom chain file for liftover (default: NULL, auto-download).
#' @param chain.source Character, source of chain file if not provided ('ensembl' or 'ucsc'; default: 'ensembl').
#' @param imputation.ind Logical, whether to include an indicator column for imputation (default: FALSE).
#' @param chrom.col Character, column name for chromosome in input data (default: 'CHR').
#' @param start.col Character, column name for start position in input data (default: 'POS').
#' @param as.granges Logical, whether to return results as a GRanges object (default: FALSE).
#' @param GRanges.style Character, style for GRanges object ('NCBI' or 'UCSC'; default: 'NCBI').
#' @param verbose Logical, whether to print progress messages (default: TRUE).
#' @param save.path Character, path to save the lifted-over summary statistics (default: NULL, no saving).
#'
#' @return A data frame or GRanges object (based on \code{as.granges}) containing the lifted-over
#'   summary statistics with updated genomic coordinates.
#'
#' @examples
#' \dontrun{
#'   sumstats <- data.frame(CHR = c(1, 2), POS = c(1000, 2000), SNP = c("rs1", "rs2"))
#'   result <- liftover(sumstats, build.from = "hg38", build.to = "hg19", save.path = "output/path.csv")
#' }
#'
#' @importFrom logger log_info log_debug log_warn log_error
#' @importFrom data.table as.data.table set setnames
#' @importFrom rtracklayer import.chain liftOver
#' @importFrom GenomeInfoDb mapGenomeBuilds
#' @importFrom GenomicRanges mcols
#'
#' @export
liftover <- function(
    sumstats,
    build.from = "hg38",
    build.to = "hg19",
    chain.file = NULL,
    chain.source = "ensembl",
    imputation.ind = FALSE,
    chrom.col = "CHR",
    start.col = "POS",
    as.granges = FALSE,
    GRanges.style = "NCBI",
    verbose = TRUE,
    save.path = NULL
) {

  # Initialize logging
  if (verbose) logger::log_info("Starting genomic coordinate liftover process.")
  start.time <- Sys.time()
  if (verbose) logger::log_debug("Input parameters: build.from={build.from}, build.to={build.to}, chain.source={chain.source}")

  # Define a variable named chain.dir to store the file directory path
  # This path points to the 'extdata' subdirectory of the R package "GWASkitR", used for storing external data files
  # system.file() function retrieves the full path to the specified directory within the R package
  chain.dir = file.path(system.file(package = "GWASkitR"), 'extdata')

  # Check if the directory specified by chain.dir exists
  # If the directory does not exist (!dir.exists(paths = chain.dir) returns TRUE), create the directory
  # dir.create() function is used to create the directory, with recursive = TRUE to create parent directories if needed
  if (!dir.exists(paths = chain.dir)) {
    # Log the attempt to create the directory
    logger::log_info("Attempting to create directory: {chain.dir}")
    dir.create(path = chain.dir, recursive = TRUE)
    # Log the successful creation of the directory
    logger::log_info("Directory created successfully: {chain.dir}")
  } else {
    # Log that the directory already exists
    logger::log_info("Directory already exists: {chain.dir}")
  }

  # Inner function for liftover logic
  liftover.inner <- function(
    sumstats, build.from, build.to, chain.file = NULL, chain.source = NULL,
    imputation.ind = TRUE, chrom.col = "CHR", start.col = "POS", end.col = start.col,
    as.granges = FALSE, GRanges.style = GRanges.style, verbose = TRUE
  ) {
    # Log start of inner function
    if (verbose) logger::log_info("Processing liftover for input data with {nrow(sumstats)} rows.")

    # Determine if start and end columns are the same
    if (start.col == end.col) {
      sestatus <- TRUE
    } else {
      sestatus <- FALSE
    }
    if (verbose) logger::log_debug("Start and end column status: same={sestatus}")

    # Validate start column existence
    if (!(start.col %in% colnames(sumstats))) {
      error_msg <- paste("Column", start.col, "does not exist in sumstats.")
      if (verbose) logger::log_error(error_msg)
      stop(error_msg)
    } else {
      # Rename start column to 'BP' if not already named as such
      if (start.col != "BP") {
        colnames(sumstats)[colnames(sumstats) == start.col] <- "BP"
        if (verbose) logger::log_info("Renamed column {start.col} to 'BP'.")
        start.col <- 'BP'
        if (sestatus) { end.col <- start.col }
      }
    }

    # Validate chain source
    chain.source <- tolower(chain.source)
    chain.error.msg <- "Chain file source must be either 'Ensembl' or 'UCSC'."
    if (length(chain.source) > 1 || !chain.source %in% c("ucsc", "ensembl")) {
      if (verbose) logger::log_error(chain.error.msg)
      stop(chain.error.msg)
    }
    if (verbose) logger::log_debug("Chain source validated: {chain.source}")

    # Determine source and target genome builds
    if (!is.null(build.from)) {
      source.ref <- if (!is.null(build.from)) {
        GenomeInfoDb::mapGenomeBuilds(genome = build.from)$ucscID[1]
      } else {
        build.from
      }
      target.ref <- if (!is.null(build.to)) {
        GenomeInfoDb::mapGenomeBuilds(genome = build.to)$ucscID[1]
      } else {
        build.to
      }
      null <- c("source.reference", "target.reference")[c(is.null(source.ref), is.null(target.ref))]
      if (length(null) > 0) {
        msg <- paste0("Could not recognize genome build of: ", paste(null, collapse = ", "), ". Will infer from data.")
        if (verbose) logger::log_warn(msg)
        if (verbose) message(msg)
      }

      # Perform liftover if source and target builds differ
      if ((!is.null(source.ref) & !is.null(target.ref)) && (source.ref != target.ref)) {
        if (verbose) logger::log_info("Performing liftover from {source.ref} to {target.ref}.")

        # Check if liftover is between supported builds (hg38 and hg19)
        if (any(!c(source.ref, target.ref) %in% c("hg38", "hg19")) || source.ref == target.ref) {
          error_msg <- "Can only perform liftover between hg38 and hg19."
          if (verbose) logger::log_error(error_msg)
          stop(error_msg)
        }

        # Convert data to GRanges object
        if (verbose) logger::log_info("Converting summary statistics to GRanges object.")
        gr <- tryCatch({
          suppressMessages({
            MungeSumstats:::to_granges(
              sumstats_dt = sumstats,
              style = GRanges.style,
              seqnames.field = chrom.col,
              start.field = start.col,
              end.field = end.col
            )
          })
        }, error = function(e) {
          error_msg <- paste("Error converting to GRanges:", e$message)
          if (verbose) logger::log_error(error_msg)
          stop(error_msg)
        })
        if (verbose) logger::log_info("Converting summary statistics to GenomicRanges.")

        # Load or download chain file for liftover
        if (verbose) logger::log_info("Loading chain file for liftover.")
        if (!is.null(chain.file)) {
          chain <- tryCatch({
            rtracklayer::import.chain(chain.file)
          }, error = function(e) {
            error_msg <- paste("Error loading chain file:", e$message)
            if (verbose) logger::log_error(error_msg)
            stop(error_msg)
          })
          if (verbose) logger::log_info("Custom chain file loaded from {chain.file}.")
        } else {
          chain <- tryCatch({
            suppressMessages({
              MungeSumstats:::get_chain_file(
                from = source.ref,
                to = target.ref,
                chain_source = chain.source,
                verbose = verbose,
                save_dir = chain.dir
              )
            })
          }, error = function(e) {
            error_msg <- paste("Error downloading chain file:", e$message)
            if (verbose) logger::log_error(error_msg)
            stop(error_msg)
          })
          if (verbose) logger::log_info("Chain file downloaded from {chain.source}.")
          if (verbose) logger::log_info("Using chain file from {chain.dir}.")
        }

        # Perform liftover
        if (verbose) logger::log_info("Performing liftover on GRanges object.")
        gr.lifted <- tryCatch({
          unlist(rtracklayer::liftOver(x = gr, chain = chain))
        }, error = function(e) {
          error_msg <- paste("Error during liftover:", e$message)
          if (verbose) logger::log_error(error_msg)
          stop(error_msg)
        })
        gr.lifted <- MungeSumstats:::granges_style(gr = gr.lifted, style = GRanges.style)
        if (verbose) logger::log_info("Liftover completed with {length(gr.lifted)} coordinates lifted.")

        # Return results in requested format
        if (as.granges) {
          if (imputation.ind) {
            GenomicRanges::mcols(gr.lifted)[["IMPUTATION_gen_build"]] <- TRUE
          }
          if (verbose) logger::log_info("Returning results as GRanges object.")
          return(gr.lifted)
        } else {
          # Convert GRanges back to data.table
          sumstats.dt <- data.table::as.data.table(gr.lifted)
          void.cols <- c("width", "strand")
          void.cols <- void.cols[void.cols %in% names(sumstats.dt)]
          for (col in void.cols) {
            data.table::set(sumstats.dt, j = col, value = NULL)
          }
          data.table::setnames(sumstats.dt, "seqnames", "CHR")
          if (start.col == end.col) {
            if ("end" %in% names(sumstats.dt)) {
              data.table::set(sumstats.dt, j = "end", value = NULL)
            }
            data.table::setnames(sumstats.dt, old = "start", new = start.col)
          } else {
            if (all(c("start", "end") %in% names(sumstats.dt))) {
              data.table::setnames(sumstats.dt, c("start", "end"), c(start.col, end.col))
            }
          }
          suppressMessages({ sumstats.dt <- MungeSumstats:::check_col_order(sumstats_dt = sumstats.dt, path = NULL)$sumstats_dt })
          if (verbose) logger::log_info("Reordering so first three column headers are SNP, CHR and BP in this order.")
          if (imputation.ind) {
            if (!"IMPUTATION_gen_build" %in% names(sumstats.dt)) {
              data.table::set(sumstats.dt, j = "IMPUTATION_gen_build", value = rep(TRUE, nrow(sumstats.dt)))
            } else {
              data.table::set(sumstats.dt, j = "IMPUTATION_gen_build", value = TRUE)
            }
          }
          if (verbose) logger::log_info("Converted lifted coordinates back to data frame.")
          return(sumstats.dt)
        }
      }
    }
    if (verbose) logger::log_warn("No liftover performed; source and target builds are the same or invalid.")
    return(data.table::as.data.table(sumstats))
  }

  # Call inner liftover function
  if (verbose) logger::log_info("Executing inner liftover function.")
  sumstats.dat <- tryCatch({
    liftover.inner(
      sumstats = sumstats,
      build.from = build.from,
      build.to = build.to,
      chain.file = chain.file,
      chain.source = chain.source,
      imputation.ind = imputation.ind,
      chrom.col = chrom.col,
      start.col = start.col,
      end.col = start.col,
      as.granges = as.granges,
      GRanges.style = GRanges.style,
      verbose = verbose
    )
  }, error = function(e) {
    error_msg <- paste("Error in liftover process:", e$message)
    if (verbose) logger::log_error(error_msg)
    stop(error_msg)
  })

  # Rename 'BP' back to 'POS' or original name
  data.table::setnames(sumstats.dat, "BP", "POS")
  if (verbose) logger::log_info("Renamed 'BP' column back to 'POS' in output data.")

  # Save results if save.path is provided
  if (!is.null(save.path)) {
    if (verbose) logger::log_info("Saving lifted-over summary statistics to {save.path}.")
    tryCatch({
      save.sumstats(sumstats.dat, save.path)
    }, error = function(e) {
      error_msg <- paste("Error saving summary statistics:", e$message)
      if (verbose) logger::log_error(error_msg)
      stop(error_msg)
    })
  } else {
    if (verbose) logger::log_info("No save path provided; results not saved to file.")
  }

  # Log completion and execution time
  end.time <- Sys.time()
  execution.time <- end.time - start.time
  if (verbose) logger::log_info("Liftover process completed in {round(execution.time, 2)} seconds.")
  if (verbose) logger::log_info(execution.time)

  return(sumstats.dat)
}
