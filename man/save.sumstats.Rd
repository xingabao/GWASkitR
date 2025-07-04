% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/save.sumstats.R
\name{save.sumstats}
\alias{save.sumstats}
\title{Save summary statistics data to a text or TSV file}
\usage{
save.sumstats(sumstats.dat, save.path)
}
\arguments{
\item{sumstats.dat}{A data frame or data table containing summary statistics to be saved.}

\item{save.path}{File path where the summary statistics will be saved. Must end with `.txt` or `.tsv`.}
}
\value{
No return value. The function is called for its side effect of writing to disk and logging.
}
\description{
This function saves a data frame or data table (such as summary statistics) to a `.txt` or `.tsv` file.
It checks the file extension, creates the target directory if necessary, writes the table in tab-delimited format,
and logs each important step using the `logger` package.
}
\examples{
\dontrun{
df <- data.frame(a = 1:3, b = c("x", "y", "z"))
save.sumstats(df, "output/sumstats.txt")
}

}
