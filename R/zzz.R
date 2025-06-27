.onLoad <- function(libname, pkgname) {

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("Installing BiocManager for Bioconductor dependencies...")
    install.packages("BiocManager")
  }

  if (!requireNamespace("MungeSumstats", quietly = TRUE)) {
    message("Installing MungeSumstats from Bioconductor...")
    BiocManager::install("MungeSumstats", update = FALSE, ask = FALSE)
  }

  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    message("Installing GenomicRanges from Bioconductor...")
    BiocManager::install("GenomicRanges", update = FALSE, ask = FALSE)
  }

  if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) {
    message("Installing GenomeInfoDb from Bioconductor...")
    BiocManager::install("GenomeInfoDb", update = FALSE, ask = FALSE)
  }

  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    message("Installing rtracklayer from Bioconductor...")
    BiocManager::install("rtracklayer", update = FALSE, ask = FALSE)
  }
}


