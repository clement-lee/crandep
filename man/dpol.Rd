% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{dpol}
\alias{dpol}
\title{Probability mass function (PMF) of Zipf-polylog distribution}
\usage{
dpol(x, alpha, theta, x_max = 100000L)
}
\arguments{
\item{x}{Vector of positive integers}

\item{alpha}{Real number greater than 1}

\item{theta}{Real number in (0, 1]}

\item{x_max}{Scalar (default 100000), positive integer limit for computing the normalising constant}
}
\value{
A numeric vector of the same length as x
}
\description{
\code{dpol} returns the PMF at x for the Zipf-polylog distribution with parameters (alpha, theta). The distribution is reduced to the discrete power law when theta = 1.
}
\details{
The PMF is proportional to x^(-alpha) * theta^x. It is normalised in order to be a proper PMF.
}
\examples{
dpol(c(1,2,3,4,5), 1.2, 0.5)
}
\seealso{
\code{\link{Spol}} for the corresponding survival function, \code{\link{dmix2}} and \code{\link{dmix3}} for the PMFs of the 2-component and 3-component discrete extreme value mixture distributions, respectively.
}
