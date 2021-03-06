% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HDGC_VAR_multiple_pairs_I0.R
\name{HDGC_VAR_multiple_pairs_I0}
\alias{HDGC_VAR_multiple_pairs_I0}
\title{Test multiple pairs Granger causality in High Dimensional Stationary VARs}
\usage{
HDGC_VAR_multiple_pairs_I0(
  data,
  GCpairs = NULL,
  GCto = NULL,
  GCfrom = NULL,
  p = 1,
  bound = 0.5 * nrow(data),
  parallel = FALSE,
  n_cores = NULL
)
}
\arguments{
\item{data}{the data matrix or object that can be coerced to a matrix.}

\item{GCpairs}{it should contain a nested list. The outer list is all the pairs to be considered. See \code{ \link{HDGC_VAR_multiple_I0}}.
The inner list contains the GCto and GCfrom vectors needed for \code{\link{HDGC_VAR_I0}}.}

\item{GCto}{all combination variables Granger caused}

\item{GCfrom}{all combination variables Granger causing}

\item{p}{lag length of the VAR}

\item{bound}{lower bound on tuning parameter lambda}

\item{parallel}{TRUE for parallel computing}

\item{n_cores}{nr of cores to use in parallel computing, default is all but one}
}
\value{
Granger causality matrix and Lasso selections are printed to the console
}
\description{
A wrapper around \code{\link{HDGC_VAR_multiple_I0}}. If GCpairs is used, the function is the same as \code{\link{HDGC_VAR_multiple_I0}}.
Alternatively, if we want to test all combinations between variables in GCto and GCfrom, these can be given directly. See Example.
}
\examples{
\dontrun{GCto = list(c("Var 1", "Var 2")); GCfrom = list(c("Var 3", "Var 4", "Var 5"))}
\dontrun{HDGC_VAR_multiple_pairs_I0(sample_dataset_I0,GCto,GCfrom,p=2)}
}
\references{
Hecq, A., Margaritella, L., Smeekes, S., "Granger Causality Testing in High-Dimensional VARs: a Post-Double-Selection Procedure." arXiv preprint arXiv:1902.10991 (2019).
}
