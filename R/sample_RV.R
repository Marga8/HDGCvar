#' @title  Dataset of simulated Realized Volatilities via HAR(1,5,22)
#'
#' @description  A dataset of 30 Realized Volatility time series of sample T=200 obtined by simulating 30 random instances from a
#' Heterogeneous Autoregressive (HAR) model with daily, weekly and monthly lags. The simulations are obtained using \code{\link[HARModel]{HARSimulate}}.
#'
#' @format A matrix of 30 columns and 200 rows where each column correspond to a single Realized Volatility time series.
#' \describe{
#' \item{T_}{The sample size}
#' \item{g}{The number of covariates}
#' \item{len}{Time series length}
#' \item{periods}{daily, weekly and monthly lags}
#' \item{coef}{Coefficients: constant, daily lag, weekly lag, monthly lag}
#' \item{errorTermSD}{standard deviation of the error term}
#' }
#' @importFrom HARModel HARSimulate
"sample_RV"

