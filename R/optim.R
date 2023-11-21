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
#' @param invts Vector of the inverse temperatures for Metropolis-coupled MCMC
#' @param name Boolean; if the column \code{name} exists, are its unique values printed?
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
                             name = TRUE) {
  if (name && !is.null(df$name)) {
    print(paste0("Data set: ", unique(as.character(df$name))))
  }
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
        invt = if (mc3) invts else 1.0
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
#' @param name Boolean; if the column \code{name} exists, are its unique values printed?
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
#' @seealso \code{\link{mcmc_mix2_wrapper}} that wraps \code{obtain_u_set_mix2} and \code{\link{mcmc_mix2}}
#' @export
obtain_u_set_mix2 <- function(df,
                              powerlaw = FALSE,
                              positive = FALSE,
                              name = TRUE,
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
  if (name && !is.null(df$name)) {
    print(paste0("Data set: ", unique(as.character(df$name))))
  }
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
      l.check = as.numeric(NA)
    )
  obj.bulk <- obj.igpd <- NULL
  both <- !inherits(obj.bulk, "try-error") && !inherits(obj.igpd, "try-error")
  ## loop
  while (u < u_max && both) {
    i <- i + 1L
    u <- u + 1L
    phiu <- mean(y > u)
    psiu <- mean(x > u)
    if (powerlaw) {
      obj.bulk <-
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
      obj.bulk <-
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
    obj.igpd <-
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
    both <- !inherits(obj.bulk, "try-error") && !inherits(obj.igpd, "try-error")
    if (both) {
      la <- obj.bulk$value + obj.igpd$value
      par.bulk <- obj.bulk$par
      if (powerlaw) {
        par.bulk <- c(par.bulk, 1.0)
      }
      par.igpd <- obj.igpd$par
      lb <-
        lpost_mix2(
          x, count, u,
          par.bulk[1], par.bulk[2],
          par.igpd[1], par.igpd[2],
          a_psiu, b_psiu,
          m_alpha, s_alpha,
          a_theta, b_theta,
          m_shape, s_shape,
          a_sigma, b_sigma,
          powerlaw = powerlaw,
          positive = positive
        )
      lc <-
        llik_bulk(
          par = par.bulk,
          x = x,
          count = count,
          v = min(x) - 1L,
          u = u,
          phil = 1.0 - phiu,
          powerlaw = powerlaw,
          positive = positive
        ) +
        llik_igpd(
          par = par.igpd,
          x = x,
          count = count,
          u = u,
          phiu = phiu
        )
      lp <- la + dunif(psiu, a_psiu, b_psiu, log = TRUE)
      l.check <- lb
      ll <- lc
      df0 <-
        data.frame(
          u = u,
          alpha = par.bulk[1],
          theta = par.bulk[2],
          shape = par.igpd[1],
          sigma = par.igpd[2],
          ll = ll,
          lp = lp,
          l.check = l.check
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
    name = name,
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
#' @param name Boolean; if the column \code{name} exists, are its unique values printed?
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
                              invts = 0.001 ^ ((0:8)/8),
                              name = TRUE) {
  obj0 <-
    obtain_u_set_mix2(df, powerlaw = (pr_power == 1.0), positive = positive, name = name)
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
#' @param name Boolean; if the column \code{name} exists, are its unique values printed?
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
                              name = TRUE,
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
  if (name && !is.null(df$name)) {
    print(paste0("Data set: ", unique(as.character(df$name))))
  }
  x <- df$x
  if (any(x <= 0L)) {
    stop("df$x has to be positive integers.")
  }
  count <- df$count
  y <- rep(x, count) # full data
  u_max <- min(max(x[x != max(x)]) - 1L, u_max)
  v_max <- min(u_max - 1L, v_max)
  v.seq <- seq(min(x) + 1L, v_max, by = 1L)
  j <- 0L
  l1 <- list()
  for (i in seq_along(v.seq)) {
    v <- v.seq[i]
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
        l.check = as.numeric(NA)
      )
    obj.pol1 <- obj.pol2 <- obj.igpd <- NULL
    three <-
      !inherits(obj.pol1, "try-error") &&
      !inherits(obj.pol2, "try-error") &&
      !inherits(obj.igpd, "try-error")
    while (u < u_max && three) {
      j <- j + 1L
      u <- u + 1L
      phi2 <- mean(y > v & y <= u)
      psi2 <- mean(x > v & x <= u)
      phiu <- mean(y > u)
      psiu <- mean(x > u)
      if (powerlaw1) {
        obj.pol1 <-
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
        obj.pol1 <-
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
        obj.pol2 <-
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
        obj.pol2 <-
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
      obj.igpd <-
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
        !inherits(obj.pol1, "try-error") &&
        !inherits(obj.pol2, "try-error") &&
        !inherits(obj.igpd, "try-error")
      if (three) {
        la <- obj.pol1$value + obj.pol2$value + obj.igpd$value
        par.pol1 <- obj.pol1$par
        if (powerlaw1) {
          par.pol1 <- c(par.pol1, 1.0)
        }
        par.pol2 <- obj.pol2$par
        if (powerlaw2) {
          par.pol2 <- c(par.pol2, 1.0)
        }
        par.igpd <- obj.igpd$par
        lb <-
          lpost_mix3(
            x, count, v, u,
            par.pol1[1], par.pol1[2],
            par.pol2[1], par.pol2[2],
            par.igpd[1], par.igpd[2],
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
            positive2 = positive2
          )
        lc <-
          llik_bulk(
            par = par.pol1,
            x = x,
            count = count,
            v = min(x) - 1,
            u = v,
            phil = phi1,
            powerlaw = powerlaw1,
            positive = positive1
          ) +
          llik_bulk(
            par = par.pol2,
            x = x,
            count = count,
            v = v,
            u = u,
            phil = phi2,
            powerlaw = powerlaw2,
            positive = positive2
          ) +
          llik_igpd(
            par = par.igpd,
            x = x,
            count = count,
            u = u,
            phiu = phiu
          )
        lp <-
          la +
          dbeta(psi1 / (1.0 - psiu), a_psi1, a_psi2, log = TRUE) +
          dunif(psiu, a_psiu, b_psiu, log = TRUE)
        l.check <- lb
        ll <- lc
        df0 <-
          data.frame(
            v = v,
            u = u,
            alpha1 = par.pol1[1],
            theta1 = par.pol1[2],
            alpha2 = par.pol2[1],
            theta2 = par.pol2[2],
            shape = par.igpd[1],
            sigma = par.igpd[2],
            ll = ll,
            lp = lp,
            l.check = l.check
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
  df1 <- dplyr::bind_rows(l1)
  df1 <- df1[df1$lp > -Inf, ]
  df1$ldiff <- max(df1$lp) - df1$lp
  df2 <- data.frame(
    powerlaw1 = powerlaw1,
    powerlaw2 = powerlaw2,
    positive1 = positive1,
    positive2 = positive2,
    name = name,
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
#' @param name Boolean; if the column \code{name} exists, are its unique values printed?
#' @return A list returned by \code{mcmc_mix3}
#' @export
mcmc_mix3_wrapper <- function(df, seed,
                              a_psi1 = 1.0,
                              a_psi2 = 1.0,
                              a_psiu = 0.001,
                              b_psiu = 0.9,
                              m_alpha = 0.00,
                              s_alpha = 0.00,
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
                              invts = 0.001 ^ ((0:8)/8),
                              name = TRUE) {
  obj0 <-
    obtain_u_set_mix3(
      df,
      powerlaw1 = powerlaw1,
      powerlaw2 = (pr_power2 == 1.0),
      positive1 = positive1,
      positive2 = positive2,
      name = name
    )
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
