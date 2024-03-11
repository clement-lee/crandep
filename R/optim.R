#' Unnormalised log-posterior density of discrete power law
#'
#' @param alpha Real number greater than 1
#' @param df A data frame with at least two columns, x & count
#' @param m_alpha Real number, mean of the prior normal distribution for alpha
#' @param s_alpha Positive real number, standard deviation of the prior normal distribution for alpha
#' @importFrom gsl hzeta
#' @importFrom stats dnorm
#' @return A real number
#' @keywords internal
lpost_pow <- function(alpha, df, m_alpha, s_alpha) {
  y <- rep(df$x, df$count)
  llik <- ifelse(alpha <= 1.0, -Inf, -alpha * sum(log(y)) - length(y) * log(gsl::hzeta(alpha, min(df$x))))
  llik + dnorm(alpha, m_alpha, s_alpha, log = TRUE)
}

#' Marginal log-likelihood and posterior density of discrete power law via numerical integration
#'
#' @param df A data frame with at least two columns, x & count
#' @param lower Real number greater than 1, lower limit for numerical integration
#' @param upper Real number greater than lower, upper limit for numerical integration
#' @param m_alpha Real number, mean of the prior normal distribution for alpha
#' @param s_alpha Positive real number, standard deviation of the prior normal distribution for alpha
#' @param by Positive real number, the width of subintervals between lower and upper, for numerical integration and posterior density evaluation
#' @importFrom pracma integral
#' @return A list: \code{log.marginal} is the marginal log-likelihood, \code{posterior} is a data frame of non-zero posterior densities
#' @export
marg_pow <- function(df, lower, upper, m_alpha, s_alpha, by = 0.001) {
  if (lower <= 1.0) {
    stop("marg_pow: lower bound of alpha has to be greater than 1.0.")
  } else if (lower >= upper) {
    stop("marg_pow: lower bound of alpha has to be smaller than the upper bound.")
  }
  optim0 <-
    optim(
      1.5,
      lpost_pow,
      df = df,
      m_alpha = m_alpha,
      s_alpha = s_alpha,
      control = list(fnscale = -1, maxit = 50000),
      method = "Brent",
      lower = lower,
      upper = upper
    )
  foo <- function(alpha, offset) {
    exp(lpost_pow(alpha, df, m_alpha, s_alpha) - offset)
  }
  l1 <- pracma::integral(foo, lower, upper, no_intervals = (upper - lower) / by, offset = optim0$value)
  lmarg1 <- log(l1) + optim0$value
  alphas <- seq(lower, upper, by = by)
  post <- data.frame(alpha = alphas, density = foo(alphas, lmarg1))
  list(
    log.marginal = lmarg1,
    posterior = post[post$density > .Machine$double.eps, ]
  )
}

#' Wrapper of lpost_pol, assuming power law (theta = 1.0)
#'
#' @param alpha A scalar, positive
#' @param ... Other arguments passed to lpost_pol
#' @return A scalar of the log-posterior density
#' @keywords internal
lpost_pol_wrapper <- function(alpha, x, count, ...) {
  llik_temp <- -Inf # won't get updated
  lpost_pol(x, count, alpha, 1.0, llik = llik_temp, invt = 1.0, ...)
}

#' Wrapper of mcmc_pol
#'
#' @param df A data frame with at least two columns, x & count
#' @param seed Integer for \code{set.seed}
#' @param alpha_init Real number greater than 1, initial value of the parameter
#' @param theta_init Real number in (0, 1], initial value of the parameter
#' @param m_alpha Real number, mean of the prior normal distribution for alpha
#' @param s_alpha Positive real number, standard deviation of the prior normal distribution for alpha
#' @param a_theta Positive real number, first parameter of the prior beta distribution for theta; ignored if pr_power = 1.0
#' @param b_theta Positive real number, second parameter of the prior beta distribution for theta; ignored if pr_power = 1.0
#' @param a_pseudo Positive real number, first parameter of the pseudoprior beta distribution for theta in model selection; ignored if pr_power = 1.0
#' @param b_pseudo Positive real number, second parameter of the pseudoprior beta distribution for theta in model selection; ignored if pr_power = 1.0
#' @param pr_power Real number in [0, 1], prior probability of the discrete power law
#' @param iter Positive integer representing the length of the MCMC output
#' @param thin Positive integer representing the thinning in the MCMC
#' @param burn Non-negative integer representing the burn-in of the MCMC
#' @param freq Positive integer representing the frequency of the sampled values being printed
#' @param mc3 Boolean, is Metropolis-coupled MCMC to be used?
#' @param invts Vector of the inverse temperatures for Metropolis-coupled MCMC; ignored if mc3 = FALSE
#' @param xmax Scalar (default 100000), positive integer limit for computing the normalising constant
#' @return A list returned by \code{mcmc_pol}
#' @export
mcmc_pol_wrapper <- function(df, seed,
                             alpha_init = 1.5,
                             theta_init = 0.5,
                             m_alpha = 0.00,
                             s_alpha = 10.0,
                             a_theta = 1.00,
                             b_theta = 1.00,
                             a_pseudo = 10.0,
                             b_pseudo = 1.0,
                             pr_power = 0.5,
                             iter = 2e+4L,
                             thin = 2e+1L,
                             burn = 1e+5L,
                             freq = 1e+3L,
                             mc3 = FALSE,
                             invts = 0.001 ^ ((0:8)/8),
                             xmax = 100000) {
  set.seed(seed)
  t0 <- system.time({
    mcmc0 <-
      mcmc_pol(
        x = df$x,
        count = df$count,
        alpha = alpha_init,
        theta = theta_init,
        a_alpha = m_alpha,
        b_alpha = s_alpha,
        a_theta = a_theta,
        b_theta = b_theta,
        a_pseudo = a_pseudo,
        b_pseudo = b_pseudo,
        pr_power = pr_power,
        iter = iter,
        thin = thin * ifelse(mc3, 2L, 1L),
        burn = burn,
        freq = freq,
        invt = if (mc3) invts else 1.0,
        xmax = xmax
      )
  })
  mcmc0$scalars$seed <- as.integer(seed)
  mcmc0$scalars$seconds <- as.numeric(t0["elapsed"])
  mcmc0
}

#' Obtain set of thresholds with high posterior density for the TZP-power-law mixture model
#'
#' \code{obtain_u_set_mix1} computes the profile posterior density of the threshold u, and subsets the thresholds (and other parameter values) with high profile values i.e. within a certain value from the maximum posterior density. The set of u can then be used for \code{\link{mcmc_mix1}}.
#' @param df A data frame with at least two columns, x & count
#' @param positive Boolean, is alpha1 positive (TRUE) or unbounded (FALSE, default)?
#' @param log_diff_max Positive real number, the value such that thresholds with profile posterior density not less than the maximum posterior density - \code{log_diff_max} will be kept
#' @param u_max Positive integer for the maximum threshold
#' @param alpha1_init Scalar, initial value of alpha1
#' @param theta1_init Scalar, initial value of theta1
#' @param alpha2_init Scalar, initial value of alpha2
#' @param a_psiu,b_psiu,m_alpha1,s_alpha1,a_theta1,b_theta1,m_alpha2,s_alpha2 Scalars, hyperparameters of the priors for the parameters
#' @param xmax Scalar (default 100000), positive integer limit for computing the normalising constant
#' @importFrom dplyr bind_rows
#' @importFrom stats optim dunif
#' @return A list: \code{u_set} is the vector of thresholds with high posterior density, \code{init} is the data frame with the maximum profile posterior density and associated parameter values, \code{profile} is the data frame with all thresholds with high posterior density and associated parameter values, \code{scalars} is the data frame with all arguments (except df)
#' @seealso \code{\link{mcmc_mix1_wrapper}} that wraps \code{obtain_u_set_mix1} and \code{\link{mcmc_mix1}}, \code{\link{obtain_u_set_mix2}} for the equivalent function for the 2-component mixture model
#' @export
obtain_u_set_mix1 <- function(df,
                              positive = FALSE,
                              u_max = 2000L,
                              log_diff_max = 11.0,
                              alpha1_init = 0.01,
                              theta1_init = exp(-1.0),
                              alpha2_init = 2.0,
                              a_psiu = 0.1,
                              b_psiu = 0.9,
                              m_alpha1 = 0.00,
                              s_alpha1 = 10.0,
                              a_theta1 = 1.00,
                              b_theta1 = 1.00,
                              m_alpha2 = 0.00,
                              s_alpha2 = 10.0,
                              xmax = 100000) {
  x <- df$x
  if (any(x <= 0L)) {
    stop("df$x has to be positive integers.")
  }
  count <- df$count
  y <- rep(x, count) # full data
  u_max <- min(max(x[x != max(x)]) - 1L, u_max)
  l1 <- list()
  i <- 0L # counter
  u <- min(x) + 1L # smallest threshold for profile
  print(paste0("Threshold = ", u))
  df0 <-
    data.frame(
      u = u,
      alpha1 = alpha1_init,
      theta1 = theta1_init,
      alpha2 = alpha2_init,
      ll = as.numeric(NA),
      lp = as.numeric(NA),
      l_check = as.numeric(NA)
    )
  obj_bulk <- obj_pol <- NULL
  both <- !inherits(obj_bulk, "try-error") && !inherits(obj_pol, "try-error")
  ## loop
  while (u < u_max && both) {
    i <- i + 1L
    u <- u + 1L
    phiu <- mean(y > u)
    psiu <- mean(x > u)
    obj_bulk <-
      try(
        optim(
          c(df0$alpha1, df0$theta1),
          fn = lpost_bulk,
          x = x,
          count = count,
          v = min(x) - 1,
          u = u,
          a_alpha = m_alpha1,
          b_alpha = s_alpha1,
          a_theta = a_theta1,
          b_theta = b_theta1,
          phil = 1.0 - phiu,
          powerlaw = FALSE,
          positive = positive,
          control = list(fnscale = -1, maxit = 50000)
        ),
        silent = TRUE
      )
    obj_pol <-
      try(
        optim(
          2.0,
          fn = lpost_pol_wrapper,
          x = x[x > u],
          count = count[x > u],
          a_alpha = m_alpha2,
          b_alpha = s_alpha2,
          a_theta = 1.0, # shouldn't matter
          b_theta = 1.0, # shouldn't matter
          powerlaw = TRUE,
          xmax = xmax,
          control = list(fnscale = -1, maxit = 50000),
          method = "Brent",
          lower = 1.00001,
          upper = 500.0
        ),
        silent = TRUE
      )
    both <- !inherits(obj_bulk, "try-error") && !inherits(obj_pol, "try-error")
    if (both) {
      la <- obj_bulk$value + obj_pol$value
      par_bulk <- obj_bulk$par
      par_pol <- c(obj_pol$par, 1.0)
      llik_temp <- -Inf # won't get updated
      lb <-
        lpost_mix1(
          x, count, u,
          par_bulk[1], par_bulk[2], par_pol[1],
          a_psiu, b_psiu,
          m_alpha1, s_alpha1,
          a_theta1, b_theta1,
          m_alpha2, s_alpha2,
          positive, xmax,
          llik = llik_temp, invt = 1.0
        )
      lc <-
        llik_bulk(
          par = par_bulk,
          x = x,
          count = count,
          v = min(x) - 1L,
          u = u,
          phil = 1.0 - phiu,
          powerlaw = FALSE,
          positive = positive
        ) +
        llik_pol(
          par = par_pol,
          x = x,
          count = count,
          powerlaw = TRUE,
          xmax = xmax
        ) +
        sum(y > u) * log(phiu)
      lp <- la +
        sum(y > u) * log(phiu) +
        dunif(psiu, a_psiu, b_psiu, log = TRUE)
      l_check <- lb
      ll <- lc
      df0 <-
        data.frame(
          u = u,
          alpha1 = par_bulk[1],
          theta1 = par_bulk[2],
          alpha2 = par_pol[1],
          ll = ll,
          lp = lp,
          l_check = l_check
        )
      df0$phiu <- phiu
      df0$psiu <- psiu
      l1[[i]] <- df0
    } else {
      print(paste0("Threshold = ", u))
    }
    if (u == u_max) {
      print(paste0("Threshold = ", u, " (max)"))
    }
  }
  cat("\n")
  df1 <- dplyr::bind_rows(l1)
  df1 <- df1[df1$lp > -Inf, ]
  df1$ldiff <- max(df1$lp) - df1$lp
  df2 <- data.frame(
    positive = positive,
    u_max = u_max,
    log_diff_max = log_diff_max,
    alpha1_init = alpha1_init,
    theta1_init = theta1_init,
    alpha2_init = alpha2_init,
    a_psiu = a_psiu,
    b_psiu = b_psiu,
    m_alpha1 = m_alpha1,
    s_alpha1 = s_alpha1,
    a_theta1 = a_theta1,
    b_theta1 = b_theta1,
    m_alpha2 = m_alpha2,
    s_alpha2 = s_alpha2
  )
  list(
    u_set = df1[df1$ldiff <= log_diff_max, ]$u,
    init = df1[df1$ldiff == 0.0, ],
    profile = df1,
    scalars = df2
  )
}

#' Wrapper of mcmc_mix1
#'
#' @param df A data frame with at least two columns, x & count
#' @param seed Integer for \code{set.seed}
#' @param a_psiu,b_psiu,m_alpha1,s_alpha1,a_theta1,b_theta1,m_alpha2,s_alpha2 Scalars, real numbers representing the hyperparameters of the prior distributions for the respective parameters. See details for the specification of the priors.
#' @param positive Boolean, is alpha1 positive (TRUE) or unbounded (FALSE)?
#' @param iter Positive integer representing the length of the MCMC output
#' @param thin Positive integer representing the thinning in the MCMC
#' @param burn Non-negative integer representing the burn-in of the MCMC
#' @param freq Positive integer representing the frequency of the sampled values being printed
#' @param mc3 Boolean, is Metropolis-coupled MCMC to be used?
#' @param invts Vector of the inverse temperatures for Metropolis-coupled MCMC; ignored if mc3 = FALSE
#' @param xmax Scalar (default 100000), positive integer limit for computing the normalising constant
#' @return A list returned by \code{mcmc_mix1}
#' @export
mcmc_mix1_wrapper <- function(df, seed,
                              a_psiu = 0.1,
                              b_psiu = 0.9,
                              m_alpha1 = 0.00,
                              s_alpha1 = 10.0,
                              a_theta1 = 1.00,
                              b_theta1 = 1.00,
                              m_alpha2 = 0.00,
                              s_alpha2 = 10.0,
                              positive = FALSE,
                              iter = 2e+4L,
                              thin = 1L,
                              burn = 1e+4L,
                              freq = 1e+2L,
                              mc3 = FALSE,
                              invts = 0.001 ^ ((0:8)/8),
                              xmax = 100000) {
  print("Obtaining profile")
  obj0 <-
    obtain_u_set_mix1(
      df, positive,
      a_psiu = a_psiu, b_psiu = b_psiu,
      m_alpha1 = m_alpha1, s_alpha1 = s_alpha1,
      a_theta1 = a_theta1, b_theta1 = b_theta1,
      m_alpha2 = m_alpha2, s_alpha2 = s_alpha2,
      xmax = xmax
    )
  print("Running MCMC")
  set.seed(seed)
  t0 <- system.time({
    mcmc0 <-
      mcmc_mix1(
        x = df$x,
        count = df$count,
        u_set = obj0$u_set,
        u = obj0$init$u,
        alpha1 = obj0$init$alpha1,
        theta1 = obj0$init$theta1,
        alpha2 = obj0$init$alpha2,
        a_psiu = a_psiu,
        b_psiu = b_psiu,
        a_alpha1 = m_alpha1,
        b_alpha1 = s_alpha1,
        a_theta1 = a_theta1,
        b_theta1 = b_theta1,
        a_alpha2 = m_alpha2,
        b_alpha2 = s_alpha2,
        positive = positive,
        iter = iter,
        thin = thin,
        burn = burn,
        freq = freq,
        invt = if (mc3) invts else 1.0,
        xmax = xmax
      )
  })
  mcmc0$scalars$seed <- as.integer(seed)
  mcmc0$scalars$seconds <- as.numeric(t0["elapsed"])
  mcmc0
}

#' Wrapper of lpost_bulk, assuming power law (theta = 1.0)
#'
#' @param alpha A scalar, positive
#' @param ... Other arguments passed to lpost_bulk
#' @return A scalar of the log-posterior density
#' @keywords internal
lpost_bulk_wrapper <- function(alpha, ...) {
  lpost_bulk(c(alpha, 1.0), ...)
}

#' Obtain set of thresholds with high posterior density for the 2-component mixture model
#'
#' \code{obtain_u_set_mix2} computes the profile posterior density of the threshold u, and subsets the thresholds (and other parameter values) with high profile values i.e. within a certain value from the maximum posterior density. The set of u can then be used for \code{\link{mcmc_mix2}}.
#' @param df A data frame with at least two columns, x & count
#' @param powerlaw Boolean, is the power law (TRUE) or polylogarithm (FALSE, default) assumed?
#' @param positive Boolean, is alpha positive (TRUE) or unbounded (FALSE, default)?
#' @param log_diff_max Positive real number, the value such that thresholds with profile posterior density not less than the maximum posterior density - \code{log_diff_max} will be kept
#' @param u_max Positive integer for the maximum threshold
#' @param alpha_init Scalar, initial value of alpha
#' @param theta_init Scalar, initial value of theta
#' @param shape_init Scalar, initial value of shape parameter
#' @param sigma_init Scalar, initial value of sigma
#' @param a_psiu,b_psiu,m_alpha,s_alpha,a_theta,b_theta,m_shape,s_shape,a_sigma,b_sigma Scalars, hyperparameters of the priors for the parameters
#' @importFrom dplyr bind_rows
#' @importFrom stats optim dunif
#' @return A list: \code{u_set} is the vector of thresholds with high posterior density, \code{init} is the data frame with the maximum profile posterior density and associated parameter values, \code{profile} is the data frame with all thresholds with high posterior density and associated parameter values, \code{scalars} is the data frame with all arguments (except df)
#' @seealso \code{\link{mcmc_mix2_wrapper}} that wraps \code{obtain_u_set_mix2} and \code{\link{mcmc_mix2}}, \code{\link{obtain_u_set_mix1}} for the equivalent function for the TZP-power-law mixture model
#' @export
obtain_u_set_mix2 <- function(df,
                              powerlaw = FALSE,
                              positive = FALSE,
                              u_max = 2000L,
                              log_diff_max = 11.0, 
                              alpha_init = 0.01,
                              theta_init = exp(-1.0),
                              shape_init = 0.1,
                              sigma_init = 1.0,
                              a_psiu = 0.001,
                              b_psiu = 0.9,
                              m_alpha = 0.00,
                              s_alpha = 10.0,
                              a_theta = 1.00,
                              b_theta = 1.00,
                              m_shape = 0.00,
                              s_shape = 10.0,
                              a_sigma = 1.00,
                              b_sigma = 0.01) {
  x <- df$x
  if (any(x <= 0L)) {
    stop("df$x has to be positive integers.")
  }
  count <- df$count
  y <- rep(x, count) # full data
  u_max <- min(max(x[x != max(x)]) - 1L, u_max)
  l1 <- list()
  i <- 0L # counter
  u <- min(x) + 1L # smallest threshold for profile
  print(paste0("Threshold = ", u))
  df0 <-
    data.frame(
      u = u,
      alpha = alpha_init,
      theta = theta_init,
      shape = shape_init,
      sigma = sigma_init,
      ll = as.numeric(NA),
      lp = as.numeric(NA),
      l_check = as.numeric(NA)
    )
  obj_bulk <- obj_igpd <- NULL
  both <- !inherits(obj_bulk, "try-error") && !inherits(obj_igpd, "try-error")
  ## loop
  while (u < u_max && both) {
    i <- i + 1L
    u <- u + 1L
    phiu <- mean(y > u)
    psiu <- mean(x > u)
    if (powerlaw) {
      obj_bulk <-
        try(
          optim(
            2.0,
            fn = lpost_bulk_wrapper,
            x = x,
            count = count,
            v = min(x) - 1,
            u = u,
            a_alpha = m_alpha,
            b_alpha = s_alpha,
            a_theta = a_theta,
            b_theta = b_theta,
            phil = 1.0 - phiu,
            powerlaw = TRUE,
            positive = positive, # overridden by TRUE powerlaw
            control = list(fnscale = -1, maxit = 50000),
            method = "Brent",
            lower = 1.00001,
            upper = 500.0
          ),
          silent = TRUE
        )
    } else {
      obj_bulk <-
        try(
          optim(
            c(df0$alpha, df0$theta),
            fn = lpost_bulk,
            x = x,
            count = count,
            v = min(x) - 1,
            u = u,
            a_alpha = m_alpha,
            b_alpha = s_alpha,
            a_theta = a_theta,
            b_theta = b_theta,
            phil = 1.0 - phiu,
            powerlaw = FALSE,
            positive = positive,
            control = list(fnscale = -1, maxit = 50000)
          ),
          silent = TRUE
        )
    }
    obj_igpd <-
      try(
        optim(
          c(df0$shape, df0$sigma),
          fn = lpost_igpd,
          x = x,
          count = count,
          u = u,
          m_shape = m_shape,
          s_shape = s_shape,
          a_sigma = a_sigma,
          b_sigma = b_sigma,
          phiu = phiu,
          control = list(fnscale = -1, maxit = 50000)
        ),
        silent = TRUE
      )
    both <- !inherits(obj_bulk, "try-error") && !inherits(obj_igpd, "try-error")
    if (both) {
      la <- obj_bulk$value + obj_igpd$value
      par_bulk <- obj_bulk$par
      if (powerlaw) {
        par_bulk <- c(par_bulk, 1.0)
      }
      par_igpd <- obj_igpd$par
      llik_temp <- -Inf # won't get updated
      lb <-
        lpost_mix2(
          x, count, u,
          par_bulk[1], par_bulk[2],
          par_igpd[1], par_igpd[2],
          a_psiu, b_psiu,
          m_alpha, s_alpha,
          a_theta, b_theta,
          m_shape, s_shape,
          a_sigma, b_sigma,
          powerlaw = powerlaw,
          positive = positive,
          llik = llik_temp,
          invt = 1.0
        )
      lc <-
        llik_bulk(
          par = par_bulk,
          x = x,
          count = count,
          v = min(x) - 1L,
          u = u,
          phil = 1.0 - phiu,
          powerlaw = powerlaw,
          positive = positive
        ) +
        llik_igpd(
          par = par_igpd,
          x = x,
          count = count,
          u = u,
          phiu = phiu
        )
      lp <- la + dunif(psiu, a_psiu, b_psiu, log = TRUE)
      l_check <- lb
      ll <- lc
      df0 <-
        data.frame(
          u = u,
          alpha = par_bulk[1],
          theta = par_bulk[2],
          shape = par_igpd[1],
          sigma = par_igpd[2],
          ll = ll,
          ll_check = llik_temp,
          lp = lp,
          l_check = l_check
        )
      df0$phiu <- phiu
      df0$psiu <- psiu
      l1[[i]] <- df0
    } else {
      print(paste0("Threshold = ", u))
    }
    if (u == u_max) {
      print(paste0("Threshold = ", u, " (max)"))
    }
  }
  cat("\n")
  df1 <- dplyr::bind_rows(l1)
  df1 <- df1[df1$lp > -Inf, ]
  df1$ldiff <- max(df1$lp) - df1$lp
  df2 <- data.frame(
    powerlaw = powerlaw,
    positive = positive,
    u_max = u_max,
    log_diff_max = log_diff_max,
    alpha_init = alpha_init,
    theta_init = theta_init,
    shape_init = shape_init,
    sigma_init = sigma_init,
    a_psiu = a_psiu,
    b_psiu = b_psiu,
    m_alpha = m_alpha,
    s_alpha = s_alpha,
    a_theta = a_theta,
    b_theta = b_theta,
    m_shape = m_shape,
    s_shape = s_shape,
    a_sigma = a_sigma,
    b_sigma = b_sigma
  )
  list(
    u_set = df1[df1$ldiff <= log_diff_max, ]$u,
    init = df1[df1$ldiff == 0.0, ],
    profile = df1,
    scalars = df2
  )
}

#' Wrapper of mcmc_mix2
#'
#' @param df A data frame with at least two columns, x & count
#' @param seed Integer for \code{set.seed}
#' @param a_psiu,b_psiu,m_alpha,s_alpha,a_theta,b_theta,m_shape,s_shape,a_sigma,b_sigma Scalars, real numbers representing the hyperparameters of the prior distributions for the respective parameters. See details for the specification of the priors.
#' @param a_pseudo Positive real number, first parameter of the pseudoprior beta distribution for theta in model selection; ignored if pr_power = 1.0
#' @param b_pseudo Positive real number, second parameter of the pseudoprior beta distribution for theta in model selection; ignored if pr_power = 1.0
#' @param pr_power Real number in [0, 1], prior probability of the discrete power law (below u)
#' @param positive Boolean, is alpha positive (TRUE) or unbounded (FALSE)?
#' @param iter Positive integer representing the length of the MCMC output
#' @param thin Positive integer representing the thinning in the MCMC
#' @param burn Non-negative integer representing the burn-in of the MCMC
#' @param freq Positive integer representing the frequency of the sampled values being printed
#' @param mc3 Boolean, is Metropolis-coupled MCMC to be used?
#' @param invts Vector of the inverse temperatures for Metropolis-coupled MCMC; ignored if mc3 = FALSE
#' @return A list returned by \code{mcmc_mix2}
#' @export
mcmc_mix2_wrapper <- function(df, seed,
                              a_psiu = 0.001,
                              b_psiu = 0.9,
                              m_alpha = 0.00,
                              s_alpha = 10.0,
                              a_theta = 1.00,
                              b_theta = 1.00,
                              m_shape = 0.00,
                              s_shape = 10.0,
                              a_sigma = 1.00,
                              b_sigma = 0.01,
                              a_pseudo = 10.0,
                              b_pseudo = 1.0,
                              pr_power = 0.5,
                              positive = FALSE,
                              iter = 2e+4L,
                              thin = 2e+1L,
                              burn = 1e+5L,
                              freq = 1e+3L,
                              mc3 = FALSE,
                              invts = 0.001 ^ ((0:8)/8)) {
  print("Obtaining profile")
  obj0 <-
    obtain_u_set_mix2(
      df, powerlaw = (pr_power == 1.0), positive = positive,
      a_psiu = a_psiu, b_psiu = b_psiu,
      m_alpha = m_alpha, s_alpha = s_alpha,
      a_theta = a_theta, b_theta = b_theta,
      m_shape = m_shape, s_shape = s_shape,
      a_sigma = a_sigma, b_sigma = b_sigma
    )
  print("Running MCMC")
  set.seed(seed)
  t0 <- system.time({
    mcmc0 <-
      mcmc_mix2(
        x = df$x,
        count = df$count,
        u_set = obj0$u_set,
        u = obj0$init$u,
        alpha = obj0$init$alpha,
        theta = obj0$init$theta,
        shape = obj0$init$shape,
        sigma = obj0$init$sigma,
        a_psiu = a_psiu,
        b_psiu = b_psiu,
        a_alpha = m_alpha,
        b_alpha = s_alpha,
        a_theta = a_theta,
        b_theta = b_theta,
        m_shape = m_shape,
        s_shape = s_shape,
        a_sigma = a_sigma,
        b_sigma = b_sigma,
        positive = positive,
        a_pseudo = a_pseudo,
        b_pseudo = b_pseudo,
        pr_power = pr_power,
        iter = iter,
        thin = thin * ifelse(mc3, 2L, 1L),
        burn = burn,
        freq = freq,
        invt = if (mc3) invts else 1.0
      )
  })
  mcmc0$scalars$seed <- as.integer(seed)
  mcmc0$scalars$seconds <- as.numeric(t0["elapsed"])
  mcmc0
}

#' Obtain set of thresholds with high posterior density for the 3-component mixture model
#'
#' \code{obtain_u_set_mix3} computes the profile posterior density of the thresholds v & u, and subsets the thresholds (and other parameter values) with high profile values i.e. within a certain value from the maximum posterior density. The sets of v & u can then be used for \code{\link{mcmc_mix3}}.
#' @param df A data frame with at least two columns, degree & count
#' @param powerlaw1 Boolean, is the power law (TRUE) or polylogarithm (FALSE, default) assumed for the left tail?
#' @param powerlaw2 Boolean, is the power law (TRUE) or polylogarithm (FALSE, default) assumed for the middle bulk?
#' @param positive1 Boolean, is alpha positive (TRUE) or unbounded (FALSE, default) for the left tail?
#' @param positive2 Boolean, is alpha positive (TRUE) or unbounded (FALSE, default) for the middle bulk?
#' @param log_diff_max Positive real number, the value such that thresholds with profile posterior density not less than the maximum posterior density - \code{log_diff_max} will be kept
#' @param v_max Positive integer for the maximum lower threshold
#' @param u_max Positive integer for the maximum upper threshold
#' @param alpha_init Scalar, initial value of alpha
#' @param theta_init Scalar, initial value of theta
#' @param shape_init Scalar, initial value of shape parameter
#' @param sigma_init Scalar, initial value of sigma
#' @param a_psi1,a_psi2,a_psiu,b_psiu,m_alpha,s_alpha,a_theta,b_theta,m_shape,s_shape,a_sigma,b_sigma Scalars, hyperparameters of the priors for the parameters
#' @importFrom dplyr bind_rows arrange desc
#' @importFrom stats optim dbeta dunif
#' @return A list: \code{v_set} is the vector of lower thresholds with high posterior density, \code{u_set} is the vector of upper thresholds with high posterior density, \code{init} is the data frame with the maximum profile posterior density and associated parameter values, \code{profile} is the data frame with all thresholds with high posterior density and associated parameter values, \code{scalars} is the data frame with all arguments (except df)
#' @seealso \code{\link{mcmc_mix3_wrapper}} that wraps \code{obtain_u_set_mix3} and \code{\link{mcmc_mix3}}
#' @export
obtain_u_set_mix3 <- function(df,
                              powerlaw1 = FALSE,
                              powerlaw2 = FALSE,
                              positive1 = FALSE,
                              positive2 = TRUE,
                              log_diff_max = 11.0,
                              v_max = 100L,
                              u_max = 2000L,
                              alpha_init = 0.01,
                              theta_init = exp(-1.0),
                              shape_init = 1.0,
                              sigma_init = 1.0,
                              a_psi1 = 1.0,
                              a_psi2 = 1.0,
                              a_psiu = 0.001,
                              b_psiu = 0.9,
                              m_alpha = 0.00,
                              s_alpha = 10.0,
                              a_theta = 1.00,
                              b_theta = 1.00,
                              m_shape = 0.00,
                              s_shape = 10.0,
                              a_sigma = 1.00,
                              b_sigma = 0.01) {
  x <- df$x
  if (any(x <= 0L)) {
    stop("df$x has to be positive integers.")
  }
  count <- df$count
  y <- rep(x, count) # full data
  u_max <- min(max(x[x != max(x)]) - 1L, u_max)
  v_max <- min(u_max - 1L, v_max)
  v_seq <- seq(min(x) + 1L, v_max, by = 1L)
  j <- 0L
  l1 <- list()
  for (i in seq_along(v_seq)) {
    v <- v_seq[i]
    phi1 <- mean(y <= v)
    psi1 <- mean(x <= v)
    u <- v + 2L
    print(paste0("v = ", v))
    df0 <-
      data.frame(
        v = v,
        u = u,
        alpha1 = alpha_init,
        theta1 = theta_init,
        alpha2 = alpha_init,
        theta2 = theta_init,
        shape = shape_init,
        sigma = sigma_init,
        ll = as.numeric(NA),
        lp = as.numeric(NA),
        l_check = as.numeric(NA)
      )
    obj_pol1 <- obj_pol2 <- obj_igpd <- NULL
    three <-
      !inherits(obj_pol1, "try-error") &&
      !inherits(obj_pol2, "try-error") &&
      !inherits(obj_igpd, "try-error")
    while (u < u_max && three) {
      u <- u + 1L
      phi2 <- mean(y > v & y <= u)
      psi2 <- mean(x > v & x <= u)
      phiu <- mean(y > u)
      psiu <- mean(x > u)
      pxis <- c(phi2, psi2, phiu, psiu)
      if (any(pxis <= 0.0 | pxis >= 1.0)) {
        print(paste0("u = ", u, ": boundary issue, no optimisation")) # no data between v & u
      } else {
        j <- j + 1L
        if (powerlaw1) {
          obj_pol1 <-
            try(
              optim(
                2.0,
                fn = lpost_bulk_wrapper,
                x = x,
                count = count,
                v = min(x) - 1,
                u = v,
                a_alpha = m_alpha,
                b_alpha = s_alpha,
                a_theta = a_theta,
                b_theta = b_theta,
                phil = phi1,
                powerlaw = TRUE,
                positive = positive1, # overridden by TRUE powerlaw
                control = list(fnscale = -1, maxit = 50000),
                method = "Brent",
                lower = 1.00001,
                upper = 500.0
              ),
              silent = TRUE
            )
        } else {
          obj_pol1 <-
            try(
              optim(
                c(df0$alpha1, df0$theta1),
                fn = lpost_bulk,
                x = x,
                count = count,
                v = min(x) - 1,
                u = v,
                a_alpha = m_alpha,
                b_alpha = s_alpha,
                a_theta = a_theta,
                b_theta = b_theta,
                phil = phi1,
                powerlaw = FALSE,
                positive = positive1,
                control = list(fnscale = -1, maxit = 50000)
              ),
            silent = TRUE
          )
        }
        if (powerlaw2) {
          obj_pol2 <-
            try(
              optim(
                2.0,
                fn = lpost_bulk_wrapper,
                x = x,
                count = count,
                v = v,
                u = u,
                a_alpha = m_alpha,
                b_alpha = s_alpha,
                a_theta = a_theta,
                b_theta = b_theta,
                phil = phi2,
                powerlaw = TRUE,
                positive = positive2, # overridden by TRUE powerlaw
                control = list(fnscale = -1, maxit = 50000),
                method = "Brent",
                lower = 1.00001,
                upper = 500.0
              ),
              silent = TRUE
            )
        } else {
          obj_pol2 <-
            try(
              optim(
                c(df0$alpha2, df0$theta2),
                fn = lpost_bulk,
                x = x,
                count = count,
                v = v,
                u = u,
                a_alpha = m_alpha,
                b_alpha = s_alpha,
                a_theta = a_theta,
                b_theta = b_theta,
                phil = phi2,
                powerlaw = FALSE,
                positive = positive2,
                control = list(fnscale = -1, maxit = 50000)
              ),
              silent = TRUE
            )
        }
        obj_igpd <-
          try(
            optim(
              c(df0$shape, df0$sigma),
              fn = lpost_igpd,
              x = x,
              count = count,
              u = u,
              m_shape = m_shape,
              s_shape = s_shape,
              a_sigma = a_sigma,
              b_sigma = b_sigma,
              phiu = phiu,
              control = list(fnscale = -1, maxit = 50000)
            ),
            silent = TRUE
          )
        three <-
          !inherits(obj_pol1, "try-error") &&
          !inherits(obj_pol2, "try-error") &&
          !inherits(obj_igpd, "try-error")
        if (three) {
          la <- obj_pol1$value + obj_pol2$value + obj_igpd$value
          par_pol1 <- obj_pol1$par
          if (powerlaw1) {
            par_pol1 <- c(par_pol1, 1.0)
          }
          par_pol2 <- obj_pol2$par
          if (powerlaw2) {
            par_pol2 <- c(par_pol2, 1.0)
          }
          par_igpd <- obj_igpd$par
          llik_temp <- -Inf # won't get updated
          lb <-
            lpost_mix3(
              x, count, v, u,
              par_pol1[1], par_pol1[2],
              par_pol2[1], par_pol2[2],
              par_igpd[1], par_igpd[2],
              a_psi1, a_psi2, a_psiu, b_psiu,
              m_alpha, s_alpha,
              a_theta, b_theta,
              m_alpha, s_alpha,
              a_theta, b_theta,
              m_shape, s_shape,
              a_sigma, b_sigma,
              powerlaw1 = powerlaw1,
              powerlaw2 = powerlaw2,
              positive1 = positive1,
              positive2 = positive2,
              llik = llik_temp,
              invt = 1.0
            )
          lc <-
            llik_bulk(
              par = par_pol1,
              x = x,
              count = count,
              v = min(x) - 1,
              u = v,
              phil = phi1,
              powerlaw = powerlaw1,
              positive = positive1
            ) +
            llik_bulk(
              par = par_pol2,
              x = x,
              count = count,
              v = v,
              u = u,
              phil = phi2,
              powerlaw = powerlaw2,
              positive = positive2
            ) +
            llik_igpd(
              par = par_igpd,
              x = x,
              count = count,
              u = u,
              phiu = phiu
            )
          lp <-
            la +
            dbeta(psi1 / (1.0 - psiu), a_psi1, a_psi2, log = TRUE) +
            dunif(psiu, a_psiu, b_psiu, log = TRUE)
          l_check <- lb
          ll <- lc
          df0 <-
            data.frame(
              v = v,
              u = u,
              alpha1 = par_pol1[1],
              theta1 = par_pol1[2],
              alpha2 = par_pol2[1],
              theta2 = par_pol2[2],
              shape = par_igpd[1],
              sigma = par_igpd[2],
              ll = ll,
              lp = lp,
              l_check = l_check
            )
          df0$phi1 <- phi1
          df0$phi2 <- phi2
          df0$phiu <- phiu
          df0$psi1 <- psi1
          df0$psi2 <- psi2
          df0$psiu <- psiu
          l1[[j]] <- df0
        }
      }
    }
  }
  df1 <- dplyr::bind_rows(l1)
  df1 <- df1[df1$lp > -Inf, ]
  df1$ldiff <- max(df1$lp) - df1$lp
  df2 <- data.frame(
    powerlaw1 = powerlaw1,
    powerlaw2 = powerlaw2,
    positive1 = positive1,
    positive2 = positive2,
    log_diff_max = log_diff_max,
    v_max = v_max,
    u_max = u_max,
    alpha_init = alpha_init,
    theta_init = theta_init,
    shape_init = shape_init,
    sigma_init = sigma_init,
    a_psi1 = a_psi1,
    a_psi2 = a_psi2,
    a_psiu = a_psiu,
    b_psiu = b_psiu,
    m_alpha = m_alpha,
    s_alpha = s_alpha,
    a_theta = a_theta,
    b_theta = b_theta,
    m_shape = m_shape,
    s_shape = s_shape,
    a_sigma = a_sigma,
    b_sigma = b_sigma
  )
  list(
    v_set = df1[df1$ldiff <= log_diff_max, ]$v,
    u_set = df1[df1$ldiff <= log_diff_max, ]$u,
    init = df1[df1$lp == max(df1$lp), ],
    profile = dplyr::arrange(df1, dplyr::desc(lp)),
    scalars = df2
  )
}

#' Wrapper of mcmc_mix3
#'
#' @param df A data frame with at least two columns, x & count
#' @param seed Integer for \code{set.seed}
#' @param a_psi1,a_psi2,a_psiu,b_psiu,m_alpha,s_alpha,a_theta,b_theta,m_shape,s_shape,a_sigma,b_sigma Scalars, real numbers representing the hyperparameters of the prior distributions for the respective parameters. See details for the specification of the priors.
#' @param a_pseudo Positive real number, first parameter of the pseudoprior beta distribution for theta2 in model selection; ignored if pr_power2 = 1.0
#' @param b_pseudo Positive real number, second parameter of the pseudoprior beta distribution for theta2 in model selection; ignored if pr_power2 = 1.0
#' @param pr_power2 Real number in [0, 1], prior probability of the discrete power law (between v and u)
#' @param powerlaw1 Boolean, is the discrete power law assumed for below v?
#' @param positive1 Boolean, is alpha1 positive (TRUE) or unbounded (FALSE)?
#' @param positive2 Boolean, is alpha2 positive (TRUE) or unbounded (FALSE)?
#' @param iter Positive integer representing the length of the MCMC output
#' @param thin Positive integer representing the thinning in the MCMC
#' @param burn Non-negative integer representing the burn-in of the MCMC
#' @param freq Positive integer representing the frequency of the sampled values being printed
#' @param mc3 Boolean, is Metropolis-coupled MCMC to be used?
#' @param invts Vector of the inverse temperatures for Metropolis-coupled MCMC; ignored if mc3 = FALSE
#' @return A list returned by \code{mcmc_mix3}
#' @export
mcmc_mix3_wrapper <- function(df, seed,
                              a_psi1 = 1.0,
                              a_psi2 = 1.0,
                              a_psiu = 0.001,
                              b_psiu = 0.9,
                              m_alpha = 0.00,
                              s_alpha = 10.0,
                              a_theta = 1.00,
                              b_theta = 1.00,
                              m_shape = 0.00,
                              s_shape = 10.0,
                              a_sigma = 1.00,
                              b_sigma = 0.01,
                              a_pseudo = 10.0,
                              b_pseudo = 1.0,
                              pr_power2 = 0.5,
                              powerlaw1 = FALSE,
                              positive1 = FALSE,
                              positive2 = TRUE,
                              iter = 2e+4L,
                              thin = 2e+1L,
                              burn = 1e+5L,
                              freq = 1e+3L,
                              mc3 = FALSE,
                              invts = 0.001 ^ ((0:8)/8)) {
  print("Obtaining profile")
  obj0 <-
    obtain_u_set_mix3(
      df,
      powerlaw1 = powerlaw1,
      powerlaw2 = (pr_power2 == 1.0),
      positive1 = positive1,
      positive2 = positive2,
      a_psi1 = a_psi1, a_psi2 = a_psi2,
      a_psiu = a_psiu, b_psiu = b_psiu,
      m_alpha = m_alpha, s_alpha = s_alpha,
      a_theta = a_theta, b_theta = b_theta,
      m_shape = m_shape, s_shape = s_shape,
      a_sigma = a_sigma, b_sigma = b_sigma
    )
  print("Running MCMC")
  set.seed(seed)
  t0 <- system.time({
    mcmc0 <-
      mcmc_mix3(
        x = df$x,
        count = df$count,
        v_set = obj0$v_set,
        u_set = obj0$u_set,
        v = obj0$init$v,
        u = obj0$init$u,
        alpha1 = obj0$init$alpha1,
        theta1 = obj0$init$theta1,
        alpha2 = obj0$init$alpha2,
        theta2 = obj0$init$theta2,
        shape = obj0$init$shape,
        sigma = obj0$init$sigma,
        a_psi1 = a_psi1,
        a_psi2 = a_psi2,
        a_psiu = a_psiu,
        b_psiu = b_psiu,
        a_alpha1 = m_alpha,
        b_alpha1 = s_alpha,
        a_theta1 = a_theta,
        b_theta1 = b_theta,
        a_alpha2 = m_alpha,
        b_alpha2 = s_alpha,
        a_theta2 = a_theta,
        b_theta2 = b_theta,
        m_shape = m_shape,
        s_shape = s_shape,
        a_sigma = a_sigma,
        b_sigma = b_sigma,
        powerlaw1 = powerlaw1,
        positive1 = positive1,
        positive2 = positive2,
        a_pseudo = a_pseudo,
        b_pseudo = b_pseudo,
        pr_power2 = pr_power2,
        iter = iter,
        thin = thin * ifelse(mc3, 2L, 1L),
        burn = burn,
        freq = freq,
        invt = if (mc3) invts else 1.0
      )
  })
  mcmc0$scalars$seed <- as.integer(seed)
  mcmc0$scalars$seconds <- as.numeric(t0["elapsed"])
  mcmc0
}
