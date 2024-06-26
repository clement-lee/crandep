# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Probability mass function (PMF) of Zipf-polylog distribution
#'
#' \code{dpol} returns the PMF at x for the Zipf-polylog distribution with parameters (alpha, theta). The distribution is reduced to the discrete power law when theta = 1.
#'
#' The PMF is proportional to x^(-alpha) * theta^x. It is normalised in order to be a proper PMF.
#' @param x Vector of positive integers
#' @param alpha Real number greater than 1
#' @param theta Real number in (0, 1]
#' @param x_max Scalar (default 100000), positive integer limit for computing the normalising constant
#' @return A numeric vector of the same length as x
#' @examples
#' dpol(c(1,2,3,4,5), 1.2, 0.5)
#' @seealso \code{\link{Spol}} for the corresponding survival function, \code{\link{dmix2}} and \code{\link{dmix3}} for the PMFs of the 2-component and 3-component discrete extreme value mixture distributions, respectively.
#' @export
dpol <- function(x, alpha, theta, x_max = 100000L) {
    .Call(`_crandep_dpol`, x, alpha, theta, x_max)
}

#' Survival function of Zipf-polylog distribution
#'
#' \code{Spol} returns the survival function at x for the Zipf-polylog distribution with parameters (alpha, theta). The distribution is reduced to the discrete power law when theta = 1.
#'
#' @param x Vector of positive integers
#' @param alpha Real number greater than 1
#' @param theta Real number in (0, 1]
#' @param x_max Scalar (default 100000), positive integer limit for computing the normalising constant
#' @return A numeric vector of the same length as x
#' @examples
#' Spol(c(1,2,3,4,5), 1.2, 0.5)
#' @seealso \code{\link{dpol}} for the corresponding probability mass function, \code{\link{Smix2}} and \code{\link{Smix3}} for the survival functions of the 2-component and 3-component discrete extreme value mixture distributions, respectively.
#' @export
Spol <- function(x, alpha, theta, x_max = 100000L) {
    .Call(`_crandep_Spol`, x, alpha, theta, x_max)
}

llik_pol <- function(par, x, count, powerlaw, x_max) {
    .Call(`_crandep_llik_pol`, par, x, count, powerlaw, x_max)
}

lpost_pol <- function(x, count, alpha, theta, a_alpha, b_alpha, a_theta, b_theta, powerlaw, x_max, llik, invt = 1.0) {
    .Call(`_crandep_lpost_pol`, x, count, alpha, theta, a_alpha, b_alpha, a_theta, b_theta, powerlaw, x_max, llik, invt)
}

#' Markov chain Monte Carlo for Zipf-polylog distribution
#'
#' \code{mcmc_pol} returns the samples from the posterior of alpha and theta, for fitting the Zipf-polylog distribution to the data x. The samples are obtained using Markov chain Monte Carlo (MCMC). In the MCMC, a Metropolis-Hastings algorithm is used.
#' @param x Vector of the unique values (positive integers) of the data
#' @param count Vector of the same length as x that contains the counts of each unique value in the full data, which is essentially rep(x, count)
#' @param alpha Real number greater than 1, initial value of the parameter
#' @param theta Real number in (0, 1], initial value of the parameter
#' @param a_alpha Real number, mean of the prior normal distribution for alpha
#' @param b_alpha Positive real number, standard deviation of the prior normal distribution for alpha
#' @param a_theta Positive real number, first parameter of the prior beta distribution for theta; ignored if pr_power = 1.0
#' @param b_theta Positive real number, second parameter of the prior beta distribution for theta; ignored if pr_power = 1.0
#' @param a_pseudo Positive real number, first parameter of the pseudoprior beta distribution for theta in model selection; ignored if pr_power = 1.0
#' @param b_pseudo Positive real number, second parameter of the pseudoprior beta distribution for theta in model selection; ignored if pr_power = 1.0
#' @param pr_power Real number in [0, 1], prior probability of the discrete power law
#' @param iter Positive integer representing the length of the MCMC output
#' @param thin Positive integer representing the thinning in the MCMC
#' @param burn Non-negative integer representing the burn-in of the MCMC
#' @param freq Positive integer representing the frequency of the sampled values being printed
#' @param invt Vector of the inverse temperatures for Metropolis-coupled MCMC
#' @param x_max Scalar, positive integer limit for computing the normalising constant
#' @param mc3_or_marg Boolean, is invt for parallel tempering / Metropolis-coupled MCMC (TRUE, default) or marginal likelihood via power posterior (FALSE)?
#' @return A list: $pars is a data frame of iter rows of the MCMC samples, $fitted is a data frame of length(x) rows with the fitted values, amongst other quantities related to the MCMC
#' @seealso \code{\link{mcmc_mix2}} and \code{\link{mcmc_mix3}} for MCMC for the 2-component and 3-component discrete extreme value mixture distributions, respectively.
#' @export
mcmc_pol <- function(x, count, alpha, theta, a_alpha, b_alpha, a_theta, b_theta, a_pseudo, b_pseudo, pr_power, iter, thin, burn, freq, invt, mc3_or_marg, x_max) {
    .Call(`_crandep_mcmc_pol`, x, count, alpha, theta, a_alpha, b_alpha, a_theta, b_theta, a_pseudo, b_pseudo, pr_power, iter, thin, burn, freq, invt, mc3_or_marg, x_max)
}

llik_bulk <- function(par, x, count, v, u, phil, powerlaw, positive) {
    .Call(`_crandep_llik_bulk`, par, x, count, v, u, phil, powerlaw, positive)
}

lpost_bulk <- function(par, x, count, v, u, phil, a_alpha, b_alpha, a_theta, b_theta, powerlaw, positive) {
    .Call(`_crandep_lpost_bulk`, par, x, count, v, u, phil, a_alpha, b_alpha, a_theta, b_theta, powerlaw, positive)
}

llik_igpd <- function(par, x, count, u, phiu) {
    .Call(`_crandep_llik_igpd`, par, x, count, u, phiu)
}

lpost_igpd <- function(par, x, count, u, m_shape, s_shape, a_sigma, b_sigma, phiu) {
    .Call(`_crandep_lpost_igpd`, par, x, count, u, m_shape, s_shape, a_sigma, b_sigma, phiu)
}

lpost_mix1 <- function(x, count, u, alpha1, theta1, alpha2, a_psiu, b_psiu, a_alpha1, b_alpha1, a_theta1, b_theta1, a_alpha2, b_alpha2, positive, x_max, llik, invt = 1.0) {
    .Call(`_crandep_lpost_mix1`, x, count, u, alpha1, theta1, alpha2, a_psiu, b_psiu, a_alpha1, b_alpha1, a_theta1, b_theta1, a_alpha2, b_alpha2, positive, x_max, llik, invt)
}

#' Markov chain Monte Carlo for TZP-power-law mixture
#'
#' \code{mcmc_mix1} returns the posterior samples of the parameters, for fitting the TZP-power-law mixture distribution. The samples are obtained using Markov chain Monte Carlo (MCMC).
#'
#' In the MCMC, a componentwise Metropolis-Hastings algorithm is used. The threshold u is treated as a parameter and therefore sampled. The hyperparameters are used in the following priors: u is such that the implied unique exceedance probability psiu ~ Uniform(a_psi, b_psi); alpha1 ~ Normal(mean = a_alpha1, sd = b_alpha1); theta1 ~ Beta(a_theta1, b_theta1); alpha2 ~ Normal(mean = a_alpha2, sd = b_alpha2)
#' @param x Vector of the unique values (positive integers) of the data
#' @param count Vector of the same length as x that contains the counts of each unique value in the full data, which is essentially rep(x, count)
#' @param u_set Positive integer vector of the values u will be sampled from
#' @param u Positive integer, initial value of the threshold
#' @param alpha1 Real number, initial value of the parameter
#' @param theta1 Real number in (0, 1], initial value of the parameter
#' @param alpha2 Real number greater than 1, initial value of the parameter
#' @param a_psiu,b_psiu,a_alpha1,b_alpha1,a_theta1,b_theta1,a_alpha2,b_alpha2 Scalars, real numbers representing the hyperparameters of the prior distributions for the respective parameters. See details for the specification of the priors.
#' @param positive Boolean, is alpha positive (TRUE) or unbounded (FALSE)?
#' @param iter Positive integer representing the length of the MCMC output
#' @param thin Positive integer representing the thinning in the MCMC
#' @param burn Non-negative integer representing the burn-in of the MCMC
#' @param freq Positive integer representing the frequency of the sampled values being printed
#' @param invt Vector of the inverse temperatures for Metropolis-coupled MCMC
#' @param x_max Scalar, positive integer limit for computing the normalising constant
#' @param mc3_or_marg Boolean, is invt for parallel tempering / Metropolis-coupled MCMC (TRUE, default) or marginal likelihood via power posterior (FALSE)?
#' @return A list: $pars is a data frame of iter rows of the MCMC samples, $fitted is a data frame of length(x) rows with the fitted values, amongst other quantities related to the MCMC
#' @seealso \code{\link{mcmc_pol}}, \code{\link{mcmc_mix2}} and \code{\link{mcmc_mix3}} for MCMC for the Zipf-polylog, and 2-component and 3-component discrete extreme value mixture distributions, respectively.
#' @export
mcmc_mix1 <- function(x, count, u_set, u, alpha1, theta1, alpha2, a_psiu, b_psiu, a_alpha1, b_alpha1, a_theta1, b_theta1, a_alpha2, b_alpha2, positive, iter, thin, burn, freq, invt, mc3_or_marg, x_max) {
    .Call(`_crandep_mcmc_mix1`, x, count, u_set, u, alpha1, theta1, alpha2, a_psiu, b_psiu, a_alpha1, b_alpha1, a_theta1, b_theta1, a_alpha2, b_alpha2, positive, iter, thin, burn, freq, invt, mc3_or_marg, x_max)
}

#' Probability mass function (PMF) of 2-component discrete extreme value mixture distribution
#'
#' \code{dmix2} returns the PMF at x for the 2-component discrete extreme value mixture distribution. The components below and above the threshold u are the (truncated) Zipf-polylog(alpha,theta) and the generalised Pareto(shape, sigma) distributions, respectively.
#' @param x Vector of positive integers
#' @param u Positive integer representing the threshold
#' @param alpha Real number, first parameter of the Zipf-polylog component
#' @param theta Real number in (0, 1], second parameter of the Zipf-polylog component
#' @param shape Real number, shape parameter of the generalised Pareto component
#' @param sigma Real number, scale parameter of the generalised Pareto component
#' @param phiu Real number in (0, 1), exceedance rate of the threshold u
#' @return A numeric vector of the same length as x
#' @seealso \code{\link{Smix2}} for the corresponding survival function, \code{\link{dpol}} and \code{\link{dmix3}} for the PMFs of the Zipf-polylog and 3-component discrete extreme value mixture distributions, respectively.
#' @export
dmix2 <- function(x, u, alpha, theta, shape, sigma, phiu) {
    .Call(`_crandep_dmix2`, x, u, alpha, theta, shape, sigma, phiu)
}

#' Survival function of 2-component discrete extreme value mixture distribution
#'
#' \code{Smix2} returns the survival function at x for the 2-component discrete extreme value mixture distribution. The components below and above the threshold u are the (truncated) Zipf-polylog(alpha,theta) and the generalised Pareto(shape, sigma) distributions, respectively.
#' @param x Vector of positive integers
#' @param u Positive integer representing the threshold
#' @param alpha Real number, first parameter of the Zipf-polylog component
#' @param theta Real number in (0, 1], second parameter of the Zipf-polylog component
#' @param shape Real number, shape parameter of the generalised Pareto component
#' @param sigma Real number, scale parameter of the generalised Pareto component
#' @param phiu Real number in (0, 1), exceedance rate of the threshold u
#' @return A numeric vector of the same length as x
#' @seealso \code{\link{dmix2}} for the corresponding probability mass function, \code{\link{Spol}} and \code{\link{Smix3}} for the survival functions of the Zipf-polylog and 3-component discrete extreme value mixture distributions, respectively.
#' @export
Smix2 <- function(x, u, alpha, theta, shape, sigma, phiu) {
    .Call(`_crandep_Smix2`, x, u, alpha, theta, shape, sigma, phiu)
}

lpost_mix2 <- function(x, count, u, alpha, theta, shape, sigma, a_psiu, b_psiu, a_alpha, b_alpha, a_theta, b_theta, m_shape, s_shape, a_sigma, b_sigma, powerlaw, positive, llik, invt = 1.0, constrained = FALSE) {
    .Call(`_crandep_lpost_mix2`, x, count, u, alpha, theta, shape, sigma, a_psiu, b_psiu, a_alpha, b_alpha, a_theta, b_theta, m_shape, s_shape, a_sigma, b_sigma, powerlaw, positive, llik, invt, constrained)
}

#' Markov chain Monte Carlo for 2-component discrete extreme value mixture distribution
#'
#' \code{mcmc_mix2} returns the posterior samples of the parameters, for fitting the 2-component discrete extreme value mixture distribution. The samples are obtained using Markov chain Monte Carlo (MCMC).
#'
#' In the MCMC, a componentwise Metropolis-Hastings algorithm is used. The threshold u is treated as a parameter and therefore sampled. The hyperparameters are used in the following priors: u is such that the implied unique exceedance probability psiu ~ Uniform(a_psi, b_psi); alpha ~ Normal(mean = a_alpha, sd = b_alpha); theta ~ Beta(a_theta, b_theta); shape ~ Normal(mean = m_shape, sd = s_shape); sigma ~ Gamma(a_sigma, scale = b_sigma). If pr_power = 1.0, the discrete power law (below u) is assumed, and the samples of theta will be all 1.0. If pr_power is in (0.0, 1.0), model selection between the polylog distribution and the discrete power law will be performed within the MCMC.
#' @param x Vector of the unique values (positive integers) of the data
#' @param count Vector of the same length as x that contains the counts of each unique value in the full data, which is essentially rep(x, count)
#' @param u_set Positive integer vector of the values u will be sampled from
#' @param u Positive integer, initial value of the threshold
#' @param alpha Real number greater than 1, initial value of the parameter
#' @param theta Real number in (0, 1], initial value of the parameter
#' @param shape Real number, initial value of the parameter
#' @param sigma Positive real number, initial value of the parameter
#' @param a_psiu,b_psiu,a_alpha,b_alpha,a_theta,b_theta,m_shape,s_shape,a_sigma,b_sigma Scalars, real numbers representing the hyperparameters of the prior distributions for the respective parameters. See details for the specification of the priors.
#' @param positive Boolean, is alpha positive (TRUE) or unbounded (FALSE)? Ignored if constrained is TRUE
#' @param a_pseudo Positive real number, first parameter of the pseudoprior beta distribution for theta in model selection; ignored if pr_power = 1.0
#' @param b_pseudo Positive real number, second parameter of the pseudoprior beta distribution for theta in model selection; ignored if pr_power = 1.0
#' @param pr_power Real number in [0, 1], prior probability of the discrete power law (below u). Overridden if constrained is TRUE
#' @param iter Positive integer representing the length of the MCMC output
#' @param thin Positive integer representing the thinning in the MCMC
#' @param burn Non-negative integer representing the burn-in of the MCMC
#' @param freq Positive integer representing the frequency of the sampled values being printed
#' @param invt Vector of the inverse temperatures for Metropolis-coupled MCMC
#' @param mc3_or_marg Boolean, is invt for parallel tempering / Metropolis-coupled MCMC (TRUE, default) or marginal likelihood via power posterior (FALSE)?
#' @param constrained Boolean, are alpha & shape constrained such that 1/shape+1 > alpha > 1 with the powerlaw assumed in the body & "continuity" at the threshold u (TRUE), or is there no constraint between alpha & shape, with the former governed by positive, and no powerlaw and continuity enforced (FALSE, default)?
#' @return A list: $pars is a data frame of iter rows of the MCMC samples, $fitted is a data frame of length(x) rows with the fitted values, amongst other quantities related to the MCMC
#' @seealso \code{\link{mcmc_pol}} and \code{\link{mcmc_mix3}} for MCMC for the Zipf-polylog and 3-component discrete extreme value mixture distributions, respectively.
#' @export
mcmc_mix2 <- function(x, count, u_set, u, alpha, theta, shape, sigma, a_psiu, b_psiu, a_alpha, b_alpha, a_theta, b_theta, m_shape, s_shape, a_sigma, b_sigma, positive, a_pseudo, b_pseudo, pr_power, iter, thin, burn, freq, invt, mc3_or_marg = TRUE, constrained = FALSE) {
    .Call(`_crandep_mcmc_mix2`, x, count, u_set, u, alpha, theta, shape, sigma, a_psiu, b_psiu, a_alpha, b_alpha, a_theta, b_theta, m_shape, s_shape, a_sigma, b_sigma, positive, a_pseudo, b_pseudo, pr_power, iter, thin, burn, freq, invt, mc3_or_marg, constrained)
}

#' Probability mass function (PMF) of 3-component discrete extreme value mixture distribution
#'
#' \code{dmix3} returns the PMF at x for the 3-component discrete extreme value mixture distribution. The component below v is the (truncated) Zipf-polylog(alpha1,theta1) distribution, between v & u the (truncated) Zipf-polylog(alpha2,theta2) distribution, and above u the generalised Pareto(shape, sigma) distribution.
#' @param x Vector of positive integers
#' @param v Positive integer representing the lower threshold
#' @param u Positive integer representing the upper threshold
#' @param alpha1 Real number, first parameter of the Zipf-polylog component below v
#' @param theta1 Real number in (0, 1], second parameter of the Zipf-polylog component below v
#' @param alpha2 Real number, first parameter of the Zipf-polylog component between v & u
#' @param theta2 Real number in (0, 1], second parameter of the Zipf-polylog component between v & u
#' @param shape Real number, shape parameter of the generalised Pareto component
#' @param sigma Real number, scale parameter of the generalised Pareto component
#' @param phi1 Real number in (0, 1), proportion of values below v
#' @param phi2 Real number in (0, 1), proportion of values between v & u
#' @param phiu Real number in (0, 1), exceedance rate of the threshold u
#' @return A numeric vector of the same length as x
#' @seealso \code{\link{Smix3}} for the corresponding survival function, \code{\link{dpol}} and \code{\link{dmix2}} for the PMFs of the Zipf-polylog and 2-component discrete extreme value mixture distributions, respectively.
#' @export
dmix3 <- function(x, v, u, alpha1, theta1, alpha2, theta2, shape, sigma, phi1, phi2, phiu) {
    .Call(`_crandep_dmix3`, x, v, u, alpha1, theta1, alpha2, theta2, shape, sigma, phi1, phi2, phiu)
}

#' Survival function of 3-component discrete extreme value mixture distribution
#'
#' \code{Smix3} returns the survival function at x for the 3-component discrete extreme value mixture distribution. The component below v is the (truncated) Zipf-polylog(alpha1,theta1) distribution, between v & u the (truncated) Zipf-polylog(alpha2,theta2) distribution, and above u the generalised Pareto(shape, sigma) distribution.
#' @param x Vector of positive integers
#' @param v Positive integer representing the lower threshold
#' @param u Positive integer representing the upper threshold
#' @param alpha1 Real number, first parameter of the Zipf-polylog component below v
#' @param theta1 Real number in (0, 1], second parameter of the Zipf-polylog component below v
#' @param alpha2 Real number, first parameter of the Zipf-polylog component between v & u
#' @param theta2 Real number in (0, 1], second parameter of the Zipf-polylog component between v & u
#' @param shape Real number, shape parameter of the generalised Pareto component
#' @param sigma Real number, scale parameter of the generalised Pareto component
#' @param phi1 Real number in (0, 1), proportion of values below v
#' @param phi2 Real number in (0, 1), proportion of values between v & u
#' @param phiu Real number in (0, 1), exceedance rate of the threshold u
#' @return A numeric vector of the same length as x
#' @seealso \code{\link{dmix3}} for the corresponding probability mass function, \code{\link{Spol}} and \code{\link{Smix2}} for the survival functions of the Zipf-polylog and 2-component discrete extreme value mixture distributions, respectively.
#' @export
Smix3 <- function(x, v, u, alpha1, theta1, alpha2, theta2, shape, sigma, phi1, phi2, phiu) {
    .Call(`_crandep_Smix3`, x, v, u, alpha1, theta1, alpha2, theta2, shape, sigma, phi1, phi2, phiu)
}

lpost_mix3 <- function(x, count, v, u, alpha1, theta1, alpha2, theta2, shape, sigma, a_psi1, a_psi2, a_psiu, b_psiu, a_alpha1, b_alpha1, a_theta1, b_theta1, a_alpha2, b_alpha2, a_theta2, b_theta2, m_shape, s_shape, a_sigma, b_sigma, powerlaw1, powerlaw2, positive1, positive2, llik, invt = 1.0) {
    .Call(`_crandep_lpost_mix3`, x, count, v, u, alpha1, theta1, alpha2, theta2, shape, sigma, a_psi1, a_psi2, a_psiu, b_psiu, a_alpha1, b_alpha1, a_theta1, b_theta1, a_alpha2, b_alpha2, a_theta2, b_theta2, m_shape, s_shape, a_sigma, b_sigma, powerlaw1, powerlaw2, positive1, positive2, llik, invt)
}

#' Markov chain Monte Carlo for 3-component discrete extreme value mixture distribution
#'
#' \code{mcmc_mix3} returns the posterior samples of the parameters, for fitting the 3-component discrete extreme value mixture distribution. The samples are obtained using Markov chain Monte Carlo (MCMC).
#'
#' In the MCMC, a componentwise Metropolis-Hastings algorithm is used. The thresholds v and u are treated as parameters and therefore sampled. The hyperparameters are used in the following priors: psi1 / (1.0 - psiu) ~ Beta(a_psi1, a_psi2); u is such that the implied unique exceedance probability psiu ~ Uniform(a_psi, b_psi); alpha1 ~ Normal(mean = a_alpha1, sd = b_alpha1); theta1 ~ Beta(a_theta1, b_theta1); alpha2 ~ Normal(mean = a_alpha2, sd = b_alpha2); theta2 ~ Beta(a_theta2, b_theta2); shape ~ Normal(mean = m_shape, sd = s_shape); sigma ~ Gamma(a_sigma, scale = b_sigma). If pr_power2 = 1.0, the discrete power law (between v and u) is assumed, and the samples of theta2 will be all 1.0. If pr_power2 is in (0.0, 1.0), model selection between the polylog distribution and the discrete power law will be performed within the MCMC.
#' @param x Vector of the unique values (positive integers) of the data
#' @param count Vector of the same length as x that contains the counts of each unique value in the full data, which is essentially rep(x, count)
#' @param v_set Positive integer vector of the values v will be sampled from
#' @param u_set Positive integer vector of the values u will be sampled from
#' @param v Positive integer, initial value of the lower threshold
#' @param u Positive integer, initial value of the upper threshold
#' @param alpha1 Real number greater than 1, initial value of the parameter
#' @param theta1 Real number in (0, 1], initial value of the parameter
#' @param alpha2 Real number greater than 1, initial value of the parameter
#' @param theta2 Real number in (0, 1], initial value of the parameter
#' @param shape Real number, initial value of the parameter
#' @param sigma Positive real number, initial value of the parameter
#' @param a_psi1,a_psi2,a_psiu,b_psiu,a_alpha1,b_alpha1,a_theta1,b_theta1,a_alpha2,b_alpha2,a_theta2,b_theta2,m_shape,s_shape,a_sigma,b_sigma Scalars, real numbers representing the hyperparameters of the prior distributions for the respective parameters. See details for the specification of the priors.
#' @param powerlaw1 Boolean, is the discrete power law assumed for below v?
#' @param positive1 Boolean, is alpha1 positive (TRUE) or unbounded (FALSE)?
#' @param positive2 Boolean, is alpha2 positive (TRUE) or unbounded (FALSE)?
#' @param a_pseudo Positive real number, first parameter of the pseudoprior beta distribution for theta2 in model selection; ignored if pr_power2 = 1.0
#' @param b_pseudo Positive real number, second parameter of the pseudoprior beta distribution for theta2 in model selection; ignored if pr_power2 = 1.0
#' @param pr_power2 Real number in [0, 1], prior probability of the discrete power law (between v and u)
#' @param iter Positive integer representing the length of the MCMC output
#' @param thin Positive integer representing the thinning in the MCMC
#' @param burn Non-negative integer representing the burn-in of the MCMC
#' @param freq Positive integer representing the frequency of the sampled values being printed
#' @param invt Vector of the inverse temperatures for Metropolis-coupled MCMC
#' @param mc3_or_marg Boolean, is invt for parallel tempering / Metropolis-coupled MCMC (TRUE, default) or marginal likelihood via power posterior (FALSE)?
#' @return A list: $pars is a data frame of iter rows of the MCMC samples, $fitted is a data frame of length(x) rows with the fitted values, amongst other quantities related to the MCMC
#' @seealso \code{\link{mcmc_pol}} and \code{\link{mcmc_mix2}} for MCMC for the Zipf-polylog and 2-component discrete extreme value mixture distributions, respectively.
#' @export
mcmc_mix3 <- function(x, count, v_set, u_set, v, u, alpha1, theta1, alpha2, theta2, shape, sigma, a_psi1, a_psi2, a_psiu, b_psiu, a_alpha1, b_alpha1, a_theta1, b_theta1, a_alpha2, b_alpha2, a_theta2, b_theta2, m_shape, s_shape, a_sigma, b_sigma, powerlaw1, positive1, positive2, a_pseudo, b_pseudo, pr_power2, iter, thin, burn, freq, invt, mc3_or_marg = TRUE) {
    .Call(`_crandep_mcmc_mix3`, x, count, v_set, u_set, v, u, alpha1, theta1, alpha2, theta2, shape, sigma, a_psi1, a_psi2, a_psiu, b_psiu, a_alpha1, b_alpha1, a_theta1, b_theta1, a_alpha2, b_alpha2, a_theta2, b_theta2, m_shape, s_shape, a_sigma, b_sigma, powerlaw1, positive1, positive2, a_pseudo, b_pseudo, pr_power2, iter, thin, burn, freq, invt, mc3_or_marg)
}

