## ---- dpol (Zipf-polylog PMF) ------------------------------------------------

test_that("dpol returns numeric vector of correct length", {
  x <- 1:5
  p <- dpol(x, alpha = 1.5, theta = 0.8)
  expect_type(p, "double")
  expect_length(p, length(x))
})

test_that("dpol values are non-negative", {
  p <- dpol(1:20, alpha = 2.0, theta = 0.9)
  expect_true(all(p >= 0))
})

test_that("dpol sums approximately to 1 over a large support", {
  # With x_max = 10000 and well-decaying parameters the total mass should be ~1
  p <- dpol(1:10000, alpha = 2.5, theta = 0.5, x_max = 10000L)
  expect_equal(sum(p), 1, tolerance = 1e-3)
})

test_that("dpol reduces to discrete power law when theta = 1", {
  x <- 1:10
  p_pol   <- dpol(x, alpha = 2.0, theta = 1.0)
  # Direct power law: p(x) proportional to x^{-alpha}
  raw  <- x^(-2.0)
  p_pl <- raw / sum(raw[1:10000])   # normalised over 1:x_max (default 100000)
  # Both should give same relative weights
  expect_equal(p_pol / p_pol[1], (x^(-2.0)) / 1^(-2.0), tolerance = 1e-4)
})

test_that("dpol is monotone decreasing for theta = 1 (power law)", {
  p <- dpol(1:20, alpha = 2.0, theta = 1.0)
  expect_true(all(diff(p) < 0))
})

## ---- Spol (Zipf-polylog survival function) ----------------------------------

test_that("Spol returns numeric vector of correct length", {
  x <- 1:5
  s <- Spol(x, alpha = 1.5, theta = 0.8)
  expect_type(s, "double")
  expect_length(s, length(x))
})

test_that("Spol values are in [0, 1]", {
  s <- Spol(1:50, alpha = 2.0, theta = 0.9)
  expect_true(all(s >= 0 & s <= 1))
})

test_that("Spol is non-increasing", {
  s <- Spol(1:30, alpha = 2.0, theta = 0.8)
  expect_true(all(diff(s) <= 0))
})

test_that("Spol and dpol are consistent: S(x) = sum_{k > x} p(k)", {
  alpha <- 2.0; theta <- 0.8
  x <- 1:10
  p <- dpol(1:10000, alpha, theta)
  # S(x) = P(X > x) = sum of p(k) for k = x+1, ..., 10000
  s_manual <- sapply(x, function(xi) sum(p[(xi + 1):10000]))
  s_func   <- Spol(x, alpha, theta)
  expect_equal(s_func, s_manual, tolerance = 1e-4)
})

## ---- dmix2 / Smix2 ----------------------------------------------------------

test_that("dmix2 returns numeric vector of correct length", {
  x   <- 1:10
  val <- dmix2(x, u = 5L, alpha = 2.0, theta = 0.8, shape = 0.1, sigma = 2.0, phiu = 0.1)
  expect_type(val, "double")
  expect_length(val, length(x))
})

test_that("dmix2 values are non-negative", {
  val <- dmix2(1:20, u = 10L, alpha = 2.0, theta = 0.9, shape = 0.2, sigma = 3.0, phiu = 0.05)
  expect_true(all(val >= 0))
})

test_that("Smix2 returns numeric vector of correct length", {
  x   <- 1:10
  val <- Smix2(x, u = 5L, alpha = 2.0, theta = 0.8, shape = 0.1, sigma = 2.0, phiu = 0.1)
  expect_type(val, "double")
  expect_length(val, length(x))
})

test_that("Smix2 values are in [0, 1]", {
  val <- Smix2(1:30, u = 15L, alpha = 2.0, theta = 0.9, shape = 0.2, sigma = 3.0, phiu = 0.05)
  expect_true(all(val >= 0 & val <= 1))
})

test_that("Smix2 is non-increasing", {
  val <- Smix2(1:30, u = 15L, alpha = 2.0, theta = 0.9, shape = 0.2, sigma = 3.0, phiu = 0.05)
  expect_true(all(diff(val) <= 0))
})

## ---- dmix3 / Smix3 ----------------------------------------------------------

test_that("dmix3 returns numeric vector of correct length", {
  x   <- 1:15
  val <- dmix3(x, v = 4L, u = 10L,
               alpha1 = 2.0, theta1 = 0.9,
               alpha2 = 1.5, theta2 = 0.7,
               shape = 0.1, sigma = 2.0,
               phi1 = 0.3, phi2 = 0.5, phiu = 0.2)
  expect_type(val, "double")
  expect_length(val, length(x))
})

test_that("dmix3 values are non-negative", {
  val <- dmix3(1:20, v = 5L, u = 12L,
               alpha1 = 2.0, theta1 = 0.9,
               alpha2 = 1.5, theta2 = 0.7,
               shape = 0.1, sigma = 2.0,
               phi1 = 0.3, phi2 = 0.5, phiu = 0.2)
  expect_true(all(val >= 0))
})

test_that("Smix3 returns numeric vector of correct length", {
  x   <- 1:15
  val <- Smix3(x, v = 4L, u = 10L,
               alpha1 = 2.0, theta1 = 0.9,
               alpha2 = 1.5, theta2 = 0.7,
               shape = 0.1, sigma = 2.0,
               phi1 = 0.3, phi2 = 0.5, phiu = 0.2)
  expect_type(val, "double")
  expect_length(val, length(x))
})

test_that("Smix3 values are in [0, 1]", {
  val <- Smix3(1:20, v = 5L, u = 12L,
               alpha1 = 2.0, theta1 = 0.9,
               alpha2 = 1.5, theta2 = 0.7,
               shape = 0.1, sigma = 2.0,
               phi1 = 0.3, phi2 = 0.5, phiu = 0.2)
  expect_true(all(val >= 0 & val <= 1))
})

test_that("Smix3 is non-increasing", {
  val <- Smix3(1:20, v = 5L, u = 12L,
               alpha1 = 2.0, theta1 = 0.9,
               alpha2 = 1.5, theta2 = 0.7,
               shape = 0.1, sigma = 2.0,
               phi1 = 0.3, phi2 = 0.5, phiu = 0.2)
  expect_true(all(diff(val) <= 0))
})

## ---- marg_pow ---------------------------------------------------------------

test_that("marg_pow errors when lower <= 1", {
  df <- data.frame(x = 1:5, count = c(10, 5, 3, 2, 1))
  expect_error(marg_pow(df, lower = 1.0, upper = 3.0), "lower bound")
  expect_error(marg_pow(df, lower = 0.5, upper = 3.0), "lower bound")
})

test_that("marg_pow errors when lower >= upper", {
  df <- data.frame(x = 1:5, count = c(10, 5, 3, 2, 1))
  expect_error(marg_pow(df, lower = 2.0, upper = 2.0), "lower bound")
  expect_error(marg_pow(df, lower = 3.0, upper = 2.0), "lower bound")
})

test_that("marg_pow returns list with log_marginal and posterior", {
  df  <- data.frame(x = 1:5, count = c(50, 20, 10, 5, 3))
  res <- marg_pow(df, lower = 1.01, upper = 5.0, by = 0.1)
  expect_type(res, "list")
  expect_named(res, c("log_marginal", "posterior"))
  expect_type(res$log_marginal, "double")
  expect_s3_class(res$posterior, "data.frame")
  expect_named(res$posterior, c("alpha", "density"))
})

test_that("marg_pow posterior densities are non-negative", {
  df  <- data.frame(x = 1:5, count = c(50, 20, 10, 5, 3))
  res <- marg_pow(df, lower = 1.01, upper = 5.0, by = 0.1)
  expect_true(all(res$posterior$density >= 0))
})

test_that("marg_pow posterior alpha values are within [lower, upper]", {
  df  <- data.frame(x = 1:5, count = c(50, 20, 10, 5, 3))
  res <- marg_pow(df, lower = 1.5, upper = 4.0, by = 0.1)
  expect_true(all(res$posterior$alpha >= 1.5 & res$posterior$alpha <= 4.0))
})
