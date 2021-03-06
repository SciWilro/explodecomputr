% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/genetics_plots.R
\name{est_lambda}
\alias{est_lambda}
\title{Estimate the inflation factor for a distribution of P-values (from GenABEL)}
\usage{
est_lambda(data, plot = FALSE, proportion = 1, method = "regression",
  filter = TRUE, df = 1, ...)
}
\arguments{
\item{data}{A vector of reals. If all are <=1, it is assumed that this is a vector of P-values, else it is treated as a vector of chi-squares}

\item{plot}{Whether the plot should be shown or not (default).}

\item{proportion}{The proportion of lowest P (or chi^2) values to be used when estimating the inflation factor lambda. Default = 1}

\item{method}{"regression" (default) or "median"}

\item{filter}{if the test statistics with 0-value of chi^2 should be excluded prior to estimation of lambda. Default = TRUE}

\item{df}{Number of degrees of freedom. Default = 1}

\item{...}{Arguments passed to plot function}
}
\value{
A list with elements estimate and se
}
\description{
Estimate the inflation factor for a distribution of P-values or 1df chi-square test. The major use of this procedure is the Genomic Control, but can also be used to visualise the distribution of P-values coming from other tests. Methods implemented include 'median' (median(chi2)/0.455...) and regression (of observed onto expected)
}

