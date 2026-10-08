context("cacheLocation")

test_that("Monte-Carlo simulations are only cached inside the temporary directory", {
  skip_if_not_installed("R.cache")
  
  # tests/testthat.R redirects the R.cache root path into tempdir(), without it
  # simulations are saved in the cache directory of the user, which is not allowed on CRAN
  root <- normalizePath(R.cache::getCacheRootPath(), winslash = "/", mustWork = FALSE)
  temp <- normalizePath(tempdir(), winslash = "/", mustWork = FALSE)
  expect_true(startsWith(root, temp))
  
  # a call with the default options saves on the file system,
  # the simulation has to be saved below the redirected root
  testfilter <- lowpassFilter::lowpassFilter(type = "bessel", param = list(pole = 4, cutoff = 0.1),
                                             sr = 1e4)
  before <- list.files(file.path(root, "stepR"))
  getCritVal(n = 137L, filter = testfilter, r = 1e3L)
  after <- list.files(file.path(root, "stepR"))
  expect_true(length(setdiff(after, before)) > 0L)
})
