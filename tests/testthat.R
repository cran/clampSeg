# keep Monte-Carlo simulations inside the session's temporary directory,
# setting the option before R.cache is loaded avoids that R.cache creates
# its default cache root in the user's home directory
options(R.cache.rootPath = file.path(tempdir(), "Rcache"))

library(testthat)
library(clampSeg)

test_check("clampSeg")
