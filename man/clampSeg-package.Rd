\name{clampSeg-package}
\alias{clampSeg-package}
\alias{clampSeg}
\docType{package}

\title{Idealisation of patch clamp recordings}

\description{Allows for idealisation (fitting) of patch clamp (ion channel) recordings by implementing the \bold{JU}mp \bold{L}ocal d\bold{E}convolution \bold{S}egmentation filter \acronym{JULES} (\cite{Pein et al.}, 2017) in the function \code{\link{jules}}. This non-parametric (model-free) segmentation method combines statistical multiresolution techniques with local deconvolution for idealising patch clamp recordings. In particular, also flickering (events on small time scales) can be detected and idealised.}

\details{
The main function of this package is \code{\link{jules}} which implements the \bold{JU}mp \bold{L}ocal d\bold{E}convolution \bold{S}egmentation filter \acronym{JULES} (\cite{Pein et al.}, 2017). It reconstructs the signal underlying the data which is assumed to be a step (piecewise constant) function, e.g. constant conductance levels are assumed. The signal is pertubed by (Gaussian) white noise and convolved with a lowpass filter, resulting in a smooth signal pertubed by correlated noise with known correlation structure. The recorded data points are modelled as sampled (digitised) recordings of this process. For more details on this model see (\cite{Pein et al.}, 2017, section II). A small example of such a recording, 3 seconds of a gramicidin A recording, is given by \code{\link{gramA}}.\cr
The filter can be created by the function \code{\link{lowpassFilter}}, currently only Bessel filters are supported. The critical value \code{q} in (\cite{Pein et al.}, 2017, (7)), the main parameter of JULES, can either be given by the user or be obtained by the function \code{\link{getCritVal}}, automatically called if required, in an universal manner by Monte-Carlo simulations such that (7) is a level \code{alpha}-test. The critical value \code{q}, or alternatively the significance level \code{alpha}, balances the risk of over- and underfitting. By default a small significance level of \code{alpha = 0.05} is chosen to guarantee that additional artificial changes are only be detected with a small probability. The critical value \code{q} and the Monte-Carlo simulations depend on the number of data points and the filter.\cr
Since a Monte-Carlo simulation lasts potentially much longer (up to several hours or days if the number of observations is in the millions) than the main calculations, multiple possibilities for saving and loading the simulations are offered. Simulations can either be saved in the workspace in the variable \code{critValStepRTab} or persistently on the file system for which the package \code{\link[R.cache]{R.cache}} is used. Moreover, storing in and loading from variables and \link[=readRDS]{RDS} files is supported. The simulation, saving and loading can be controlled by the argument \code{option}. By default simulations will be saved in the workspace and on the file system. For more details and for how simulation can be removed see the documentation of \code{\link{getCritVal}}.\cr 
The detection and estimation step of \acronym{JULES} can be obtained separately by the functions \code{\link{stepDetection}} and \code{\link{deconvolveLocally}}, respectively.
}


\encoding{utf8}
\references{
Pein, F., Tecuapetla-Gómez, I., \enc{Schütte}{Schuette}, O. Steinem, C., and Munk, A. (2017) Fully-automatic multiresolution idealization for filtered ion channel recordings: flickering event detection. In preparation.

Hotz, T., \enc{Schütte}{Schuette}, O., Sieling, H., Polupanow, T., Diederichsen, U., Steinem, C., and Munk, A. (2013) Idealizing ion channel recordings by a jump segmentation multiresolution filter. \emph{IEEE Transactions on NanoBioscience} \bold{12}(4), 376--386.

Frick, K., Munk, A. and Sieling, H. (2014) Multiscale change-point inference. With discussion and rejoinder by the authors. \emph{Journal of the Royal Statistical Society, Series B} \bold{76}(3), 495--580.

Pein, F., Sieling, H. and Munk, A. (2016) Heterogeneous change point inference. \emph{Journal of the Royal Statistical Society, Series B}, early view.
}

\keyword{package,nonparametric}

\seealso{\code{\link{jules}}, \code{\link{critVal}}, \code{\link{lowpassFilter}}, \code{\link{gramA}}, \code{\link{deconvolveLocally}}, \code{\link{stepDetection}}}

\examples{
## idealisation of the gramicidin A recordings given by gramA
# the used filter
filter <- lowpassFilter(type = "bessel", param = list(pole = 4L, cutoff = 1e3 / 1e4),
                        sr = 1e4)

# the corresponding time points
time <- 9 + seq(along = gramA) / filter$sr

# plot of the data as in (Pein et al., 2017, figure 1 lower panel)
plot(time, gramA, pch = ".", col = "grey30", ylim = c(20, 50),
     ylab = "Conductance in pS", xlab = "Time in s")
\donttest{
# idealisation by JULES
# this call requires a Monte-Carlo simulation
# and therefore might last a few minutes,
# progress of the Monte-Carlo simulation is reported
idealisation <- jules(gramA, filter = filter, startTime = 9, messages = 100)

# this second call should be much faster
# as the previous Monte-Carlo simulation will be loaded
jules(gramA, filter = filter, startTime = 9)

# add idealisation to the plot
lines(idealisation, col = "#FF0000", lwd = 3)

# much larger significance level alpha for a larger detection power,
# but also with the risk of detecting additional artefacts
# in this example much more changes are detected,
# most of them are probably artefacts, but for instance the event at 11.36972
# might be an additional small event that was missed before
jules(gramA, filter = filter, alpha = 0.9, startTime = 9)

# getCritVal was called in jules, can be called explicitly
# for instance outside of a for loop to save computation time
q <- getCritVal(length(gramA), filter = filter)
identical(jules(gramA, q = q, filter = filter, startTime = 9), idealisation)

# both steps of JULES can be called separately
fit <- stepDetection(gramA, filter = filter, startTime = 9)
identical(deconvolveLocally(fit, data = gramA, filter = filter, startTime = 9),
          idealisation)

# more detailed output
each <- jules(gramA, filter = filter, startTime = 9, output = "each")
every <- jules(gramA, filter = filter, startTime = 9, output = "every")

identical(idealisation, each$idealization)
idealisationEvery <- every$idealization[[3]]
attr(idealisationEvery, "noDeconvolution") <- attr(every$idealization,
                                                   "noDeconvolution")
identical(idealisation, idealisationEvery)

identical(each$fit, fit)
identical(every$fit, fit)


## zoom into a single event, (Pein et al., 2017, figure 2 lower left panel)
plot(time, gramA, pch = 16, col = "grey30", ylim = c(20, 50),
     xlim = c(10.40835, 10.4103), ylab = "Conductance in pS", xlab = "Time in s")

# relevant part of the idealisation
cps <- idealisation$leftEnd[8:9]
levels <- idealisation$value[7:9]
t <- seq(cps[1] - 0.0009, cps[2] + 0.0023, 1e-6)

# idealisation
lines(t, ifelse(t < cps[1], rep(levels[1], length(t)),
                ifelse(t < cps[2], rep(levels[2], length(t)),
                       rep(levels[3], length(t)))),
      col = "#FF0000", lwd = 3)

# idealisation convolved with the filter
lines(t, levels[1] * (1 - filter$truncatedStepfun(t - cps[1])) + 
        levels[2] * (filter$truncatedStepfun(t - cps[1]) - 
                       filter$truncatedStepfun(t - cps[2])) + 
        levels[3] * filter$truncatedStepfun(t - cps[2]),
      col = "#770000", lwd = 3)

# fit prior to the deconvolution step
# does not fit the recorded data points appropriately
cps <- fit$leftEnd[8:9]
levels <- fit$value[7:9]
t <- seq(cps[1] - 0.0009, cps[2] + 0.0023, 1e-6)

# fit
lines(t, ifelse(t < cps[1], rep(levels[1], length(t)),
                ifelse(t < cps[2], rep(levels[2], length(t)),
                       rep(levels[3], length(t)))),
      col = "blue", lwd = 3)

# fit convolved with the filter
lines(t, levels[1] * (1 - filter$truncatedStepfun(t - cps[1])) + 
        levels[2] * (filter$truncatedStepfun(t - cps[1]) - 
                       filter$truncatedStepfun(t - cps[2])) + 
        levels[3] * filter$truncatedStepfun(t - cps[2]),
      col = "darkblue", lwd = 3)


## zoom into a single jump
plot(time, gramA, pch = 16, col = "grey30", ylim = c(20, 50),
     xlim = c(9.6476, 9.6496), ylab = "Conductance in pS", xlab = "Time in s")

# relevant part of the idealisation
cp <- idealisation$leftEnd[2]
levels <- idealisation$value[1:2]
t <- seq(cp - 0.0009, cp + 0.0023, 1e-6)

# idealisation
lines(t, ifelse(t < cp, rep(levels[1], length(t)), rep(levels[2], length(t))),
      col = "#FF0000", lwd = 3)

# idealisation convolved with the filter
lines(t, levels[1] * (1 - filter$stepfun(t - cp)) + levels[2] * filter$stepfun(t - cp),
      col = "#770000", lwd = 3)

# idealisation with a wrong filter
# does not fit the recorded data points appropriately
wrongFilter <- lowpassFilter(type = "bessel",
                             param = list(pole = 6L, cutoff = 0.2),
                             sr = 1e4)
                             
# Monte-Carlo simulation depend on the number of observations and on the filter
# hence a simulation is required again (if called for the first time)
# to save some time the number of iterations is reduced to r = 1e3
# hence the critical value is computed with less precision
# In general, r = 1e3 is enough for a first impression
# for a detailed analysis r = 1e4 is suggested
idealisationWrong <- jules(gramA, filter = wrongFilter, startTime = 9, 
                           r = 1e3, messages = 100)

# relevant part of the idealisation
cp <- idealisationWrong$leftEnd[2]
levels <- idealisationWrong$value[1:2]
t <- seq(cp - 0.0012, cp + 0.0023, 1e-6)

# idealisation
lines(t, ifelse(t < cp, rep(levels[1], length(t)), rep(levels[2], length(t))),
      col = "blue", lwd = 3)

# idealisation convolved with the filter
lines(t, levels[1] * (1 - filter$stepfun(t - cp)) + levels[2] * filter$stepfun(t - cp),
      col = "darkblue", lwd = 3)


# simulation for a larger number of observations can be used (nq = 3e4)
# does not require a new simulation as the simulation from above will be used
# (if the previous call was executed first)
jules(gramA[1:2.99e4], filter = wrongFilter, startTime = 9, 
      nq = 3e4, r = 1e3, messages = 100)

# simulation of type "vectorIncreased" for n1 observations can only be reused 
# for n2 observations if as.integer(log2(n1)) == as.integer(log2(n2))
# no simulation is required, since a simulation of type "matrixIncreased"
# will be loaded from the fileSystem
# this call also saves a simulation of type "vectorIncreased" in the workspace  
jules(gramA[1:1e4], filter = filter, startTime = 9,
      nq = 3e4, messages = 100, r = 1e3)
# here a new simulation is required
# (if no appropriate simulation is saved from a call outside of this file)
jules(gramA[1:1e3], filter = filter, startTime = 9,
      nq = 3e4, messages = 100, r = 1e3,
      options = list(load = list(workspace = c("vector", "vectorIncreased"))))

# the above calls saved and (attempted to) load Monte-Carlo simulations
# in the following call the simulations will neither be saved nor loaded
jules(gramA, filter = filter, startTime = 9, messages = 100, r = 1e3,
      options = list(load = list(), save = list()))

# only simulations of type "vector" and "vectorInceased" will only be in and
# loaded from the workspace, but no simulations of type "matrix" and
# "matrixIncreased" on the file system
jules(gramA, filter = filter, startTime = 9, messages = 100,
      options = list(load = list(workspace = c("vector", "vectorIncreased")),
                     save = list(workspace = c("vector", "vectorIncreased"))))

# explicit Monte-Carlo simulations, not recommended
stat <- stepR::monteCarloSimulation(n = length(gramA), , family = "mDependentPS",
                                    filter = filter, output = "maximum",
                                    r = 1e3, messages = 100)
jules(gramA, filter = filter, startTime = 9, stat = stat)

# with given standard deviation
sd <- stepR::sdrobnorm(gramA, lag = filter$len + 1)
identical(jules(gramA, filter = filter, startTime = 9, sd = sd), idealisation)

# with less regularisation of the correlation matrix
jules(gramA, filter = filter, startTime = 9, regularization = 0.5)

# with estimation of the level of long segments by the mean
# but requiring 30 observations for it
jules(gramA, filter = filter, startTime = 9,
      localEstimate = mean, thresholdLongSegment = 30)

# with one refinement step less, but with a larger grid
# progress of the deconvolution is reported
# potential warning for no deconvolution is suppressed
jules(gramA, filter = filter, startTime = 9,
      gridSize = c(1 / filter$sr, 1 / 10 / filter$sr),
      windowFactorRefinement = 2, report = TRUE,
      suppressWarningNoDeconvolution = TRUE)
}}