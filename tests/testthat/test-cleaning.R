test_that("conditional_change replaces matched strings", {
  expect_equal(crandep:::conditional_change("linkingto", "linkingto", "linking to"), "linking to")
  expect_equal(crandep:::conditional_change("imports", "linkingto", "linking to"), "imports")
})

test_that("conditional_change works on vectors", {
  x <- c("linkingto", "imports", "depends")
  result <- crandep:::conditional_change(x, "linkingto", "linking to")
  expect_equal(result, c("linking to", "imports", "depends"))
})

test_that("check_dep_word handles 'all' (case-insensitive)", {
  expected <- c("Depends", "Imports", "LinkingTo", "Suggests", "Enhances")
  expect_equal(crandep:::check_dep_word("all"), expected)
  expect_equal(crandep:::check_dep_word("ALL"), expected)
  expect_equal(crandep:::check_dep_word("All"), expected)
})

test_that("check_dep_word handles 'strong' (case-insensitive)", {
  expected <- c("Depends", "Imports", "LinkingTo")
  expect_equal(crandep:::check_dep_word("strong"), expected)
  expect_equal(crandep:::check_dep_word("STRONG"), expected)
})

test_that("check_dep_word normalises case for individual dependency types", {
  expect_equal(crandep:::check_dep_word("depends"), "Depends")
  expect_equal(crandep:::check_dep_word("imports"), "Imports")
  expect_equal(crandep:::check_dep_word("suggests"), "Suggests")
  expect_equal(crandep:::check_dep_word("enhances"), "Enhances")
})

test_that("check_dep_word handles LinkingTo variants", {
  expect_equal(crandep:::check_dep_word("LinkingTo"), "LinkingTo")
  expect_equal(crandep:::check_dep_word("linkingto"), "LinkingTo")
  expect_equal(crandep:::check_dep_word("linking to"), "LinkingTo")
  expect_equal(crandep:::check_dep_word("linking_to"), "LinkingTo")
})

test_that("check_dep_word errors on 'Reverse' prefix", {
  expect_error(crandep:::check_dep_word("Reverse imports"), "Reverse")
  expect_error(crandep:::check_dep_word("reverse depends"), "Reverse")
})

test_that("check_dep_word errors on invalid dependency word", {
  expect_error(crandep:::check_dep_word("foo"))
  expect_error(crandep:::check_dep_word("Links"))
})

test_that("check_dep_word handles multiple valid types", {
  result <- crandep:::check_dep_word(c("Imports", "Depends"))
  expect_equal(result, c("Imports", "Depends"))
})

test_that("get_dep_vec returns NA for NA input", {
  result <- crandep:::get_dep_vec(NA)
  expect_true(is.na(result))
})

test_that("get_dep_vec extracts package names from simple string", {
  result <- crandep:::get_dep_vec("dplyr, ggplot2, stringr")
  expect_equal(sort(result), sort(c("dplyr", "ggplot2", "stringr")))
})

test_that("get_dep_vec strips version specifications", {
  result <- crandep:::get_dep_vec("dplyr (>= 1.0.0), ggplot2 (>= 3.0)")
  expect_equal(sort(result), sort(c("dplyr", "ggplot2")))
})

test_that("get_dep_vec returns NA when only R is listed", {
  result <- crandep:::get_dep_vec("R (>= 3.4)")
  expect_true(is.na(result))
})

test_that("get_dep_vec drops R from a mixed list", {
  result <- crandep:::get_dep_vec("R (>= 3.4), dplyr, ggplot2")
  expect_false("R" %in% result)
  expect_true("dplyr" %in% result)
  expect_true("ggplot2" %in% result)
})
