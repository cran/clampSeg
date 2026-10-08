\name{createLocalList}
\alias{createLocalList}

\title{Precomputes quantities for hilde and improveSmallScales}
\description{Allows to precompute quantities that can be passed to \code{\link{hilde}} and \code{\link{improveSmallScales}}, all arguments have to coincide with the corresponding ones in those function calls. Doing so saves run time if those functions are called multiple times (with the same arguments) as then the quantities have to be computed only once.}

\usage{
createLocalList(filter, method = c("2Param", "LR"),
               lengths = if (method == "LR") 1:20 else 1:65)
}

\arguments{
\item{filter}{an object of class \code{\link{lowpassFilter}} giving the used analogue lowpass filter}
\item{method}{the testing method for short events, \code{"2Param"} allows for heterogeneous noise, \code{"LR"} assumes homogeneous noise}
\item{lengths}{a vector of integers giving the lengths on which tests will be performed to detect short events, should be chosen such that events on larger scales are already detected by the previous \code{\link{jsmurf}} step}
}

\value{
An object of class \code{"localList"}, contains computed quantities used by \code{\link{hilde}} and \code{\link{improveSmallScales}}
}

\seealso{\code{\link{hilde}}, \code{\link{improveSmallScales}}, \code{\link{lowpassFilter}}}

\examples{
\dontshow{
# keep Monte-Carlo simulations inside the session's temporary directory,
# setting the option before R.cache is loaded avoids that R.cache creates
# its default cache root in the user's home directory
savePathRcache <- getOption("R.cache.rootPath")
options(R.cache.rootPath = file.path(tempdir(), "Rcache"))
}
# the used filter of the gramicidin A recordings given by gramA
filter <- lowpassFilter(type = "bessel", param = list(pole = 4L, cutoff = 1e3 / 1e4),
                        sr = 1e4)

\donttest{
localList <- createLocalList(filter = filter, method = "LR")

# idealisation by HILDE using localList
# this call requires a Monte-Carlo simulation
# and therefore might last a few minutes,
# progress of the Monte-Carlo simulation is reported
# r = 100 instead of the defaults r = 1e4 for q1 and r = 1e3 for q2
# is used to reduce the run time,
# this is okay for illustration purposes, but for precise results
# the default number of Monte-Carlo simulations is recommended
idealisation <- hilde(gramA, filter = filter, family = "jsmurfPS", method = "LR",
                      startTime = 9, messages = 10, localList = localList, r = 100)

identical(hilde(gramA, filter = filter, family = "jsmurfPS", method = "LR",
                startTime = 9, r = 100), idealisation)

# for HILDE allowing heterogeneous noise, here for only few lengths
localList <- createLocalList(filter = filter, method = "2Param", lengths = c(3, 10, 25))
# localList has to be computed with the same filter, method and lengths
# as used in hilde, see the examples of hilde for HILDE allowing heterogeneous noise
}
\dontshow{
unlink(file.path(tempdir(), "Rcache"), recursive = TRUE, force = TRUE)
options(R.cache.rootPath = savePathRcache)
}
}

\keyword{nonparametric}
