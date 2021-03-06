% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{read_sfs}
\alias{read_sfs}
\title{Read 'flattened' site frequency spectrum (SFS) from a file}
\usage{
read_sfs(ff, dims, dtype = double, bootstraps = TRUE, repolarize = FALSE,
  ...)
}
\arguments{
\item{ff}{filename}

\item{dims}{sample sizes (number of chromosomes) along each dimension of SFS}

\item{dtype}{what sort of frequences to expect; best let the procedure auto-detect}

\item{bootstraps}{logical; if \code{TRUE}, read bootstrap replicates (one per line) if they are present}

\item{repolarize}{logical; if \code{TRUE}, swap ancestral and derived states (see Details)}

\item{...}{ignored}
}
\value{
a site frequency spectrum (SFS): a k-dimensional array representing the joint frequencies of derived alleles in each of k populations
}
\description{
Read 'flattened' site frequency spectrum (SFS) from a file
}
\details{
The SFS is expected to be provided in a text file as space-separated numbers (integers or floating-point), in row-major order.  
An \emph{unfolded} SFS (ie. polarized against the derived allele) is expected, although this may be relaxed in future
  
NB: Recent versions of \code{ANGSD} apparently get the ancestral and derived alleles backwards.  Use \code{repolarize = TRUE} to correct
this issue at runtime.
}
