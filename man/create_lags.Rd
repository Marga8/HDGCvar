% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_lags.R
\name{create_lags}
\alias{create_lags}
\title{Lags creation}
\usage{
create_lags(y, p = 1, include.original = TRUE, trim = TRUE)
}
\arguments{
\item{y}{vector or matrix to lag}

\item{p}{lag length}

\item{include.original}{logical, if TRUE the original series are left inside the output matrix.}

\item{trim}{logical, if TRUE the initial NAs due to the lag gets trimmed}
}
\value{
the lagged vector or matrix
}
\description{
Creates a lag matrix of order p. Can include or exclude original series, and trim the NAs in the
start of the sample.
}
\keyword{internal}
