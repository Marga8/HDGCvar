#' @title  Unit Root Sample Data for Examples
#'
#' @description  A dataset of 30 non-stationary time series of sample T=200 obtined by reverting a 30 stationary time series obtained through a non-sparse VAR(1) DGP.
#' This means effectively that the series are obtained from a VAR(2) in levels. All series have a unit root.
#'
#'
#' @format A matrix of 30 columns and 200 rows where each column correspond to a single non-stationary time series.
#' \describe{
#' \item{T_}{The sample size}
#' \item{g}{The number of covariates}
#' }
"sample_dataset_I1"

