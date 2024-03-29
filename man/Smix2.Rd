% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{Smix2}
\alias{Smix2}
\title{Survival function of 2-component discrete extreme value mixture distribution}
\usage{
Smix2(x, u, alpha, theta, shape, sigma, phiu)
}
\arguments{
\item{x}{Vector of positive integers}

\item{u}{Positive integer representing the threshold}

\item{alpha}{Real number, first parameter of the Zipf-polylog component}

\item{theta}{Real number in (0, 1], second parameter of the Zipf-polylog component}

\item{shape}{Real number, shape parameter of the generalised Pareto component}

\item{sigma}{Real number, scale parameter of the generalised Pareto component}

\item{phiu}{Real number in (0, 1), exceedance rate of the threshold u}
}
\value{
A numeric vector of the same length as x
}
\description{
\code{Smix2} returns the survival function at x for the 2-component discrete extreme value mixture distribution. The components below and above the threshold u are the (truncated) Zipf-polylog(alpha,theta) and the generalised Pareto(shape, sigma) distributions, respectively.
}
\seealso{
\code{\link{dmix2}} for the corresponding probability mass function, \code{\link{Spol}} and \code{\link{Smix3}} for the survival functions of the Zipf-polylog and 3-component discrete extreme value mixture distributions, respectively.
}
