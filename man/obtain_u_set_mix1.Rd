% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optim.R
\name{obtain_u_set_mix1}
\alias{obtain_u_set_mix1}
\title{Obtain set of thresholds with high posterior density for the TZP-power-law mixture model}
\usage{
obtain_u_set_mix1(
  df,
  positive = FALSE,
  u_max = 2000L,
  log_diff_max = 11,
  alpha1_init = 0.01,
  theta1_init = exp(-1),
  alpha2_init = 2,
  a_psiu = 0.1,
  b_psiu = 0.9,
  m_alpha1 = 0,
  s_alpha1 = 10,
  a_theta1 = 1,
  b_theta1 = 1,
  m_alpha2 = 0,
  s_alpha2 = 10,
  x_max = 1e+05
)
}
\arguments{
\item{df}{A data frame with at least two columns, x & count}

\item{positive}{Boolean, is alpha1 positive (TRUE) or unbounded (FALSE, default)?}

\item{u_max}{Positive integer for the maximum threshold}

\item{log_diff_max}{Positive real number, the value such that thresholds with profile posterior density not less than the maximum posterior density - \code{log_diff_max} will be kept}

\item{alpha1_init}{Scalar, initial value of alpha1}

\item{theta1_init}{Scalar, initial value of theta1}

\item{alpha2_init}{Scalar, initial value of alpha2}

\item{a_psiu, b_psiu, m_alpha1, s_alpha1, a_theta1, b_theta1, m_alpha2, s_alpha2}{Scalars, hyperparameters of the priors for the parameters}

\item{x_max}{Scalar (default 100000), positive integer limit for computing the normalising constant}
}
\value{
A list: \code{u_set} is the vector of thresholds with high posterior density, \code{init} is the data frame with the maximum profile posterior density and associated parameter values, \code{profile} is the data frame with all thresholds with high posterior density and associated parameter values, \code{scalars} is the data frame with all arguments (except df)
}
\description{
\code{obtain_u_set_mix1} computes the profile posterior density of the threshold u, and subsets the thresholds (and other parameter values) with high profile values i.e. within a certain value from the maximum posterior density. The set of u can then be used for \code{\link{mcmc_mix1}}.
}
\seealso{
\code{\link{mcmc_mix1_wrapper}} that wraps \code{obtain_u_set_mix1} and \code{\link{mcmc_mix1}}, \code{\link{obtain_u_set_mix2}} for the equivalent function for the 2-component mixture model
}
