\name{lowpassFilter}
\alias{lowpassFilter}
\alias{print.lowpassFilter}

\title{Lowpass filtering}

\description{
Create lowpass filter.
}

\usage{
lowpassFilter(type = c("bessel"), param, sr = 1, len = NULL, shift = 0.5)
\method{print}{lowpassFilter}(x, \dots)
}

\arguments{
\item{type}{a string specifying the type of the filter, currently only Bessel filters are supported}
\item{param}{a \code{\link{list}} specifying the parameters of the filter depending on \code{type}. For \code{"bessel"} the entries \code{pole} and \code{cutoff} have to be specified and no other named entries are allowed. \code{pole} has to be a single integer giving the number of poles (order). \code{cutoff} has to be a single positive numeric not larger than \code{1} giving the normalized cutoff frequency, i.e. the cutoff frequency (in the temporal domain) of the filter divided by the sampling rate}
\item{sr}{a single numeric giving the sampling rate}
\item{len}{a single integer giving the filter length of the truncated and digitised filter, see \cite{Value} for more details. By default (\code{NULL}) chosen such that the autocorrelation function is below \code{1e-3} at \code{len / sr} and all lager lags \code{(len + i) / sr}, with \code{i} a positive integer}
\item{shift}{a single numeric between \code{0} and \code{1} giving a shift for the digitised filter, i.e. kernel and step are obtained at \code{(0:len + shift) / sr} from the corresponding functions}
\item{x}{the object}
\item{\dots}{for generic methods only}
}

\value{
An object of \code{\link{class}} \code{lowpassFilter}, i.e. a \code{\link{list}} that contains
\describe{
\item{\code{"type"}, \code{"param"}, \code{"sr"}, \code{"len"}}{the corresponding arguments}
\item{\code{"kernfun"}}{the kernel function of the filter, obtained as the Laplace transform of the corresponding transfer function}
\item{\code{"stepfun"}}{the step-response of the filter, i.e. the antiderivative of the filter kernel}
\item{\code{"acfun"}}{the autocorrelation function, i.e. the convolution of the filter kernel with itself}          
\item{\code{"truncatedKernfun"}}{the kernel function of the at \code{len / sr} truncated filter, i.e. \code{kernfun} truncated and rescaled such that the new kernel still integrates to \code{1}}
\item{\code{"truncatedStepfun"}}{the step-response of the at \code{len / sr} truncated filter, i.e. the antiderivative of the kernel of the truncated filter}
\item{\code{"truncatedAcfun"}}{the autocorrelation function of the at \code{len / sr} truncated filter, i.e. the convolution of the kernel of the truncated filter with itself}
\item{\code{"kern"}}{the digitised filter kernel normalised to one, i.e. \code{kernfun((0:len + shift) / sr) / sum(kernfun((0:len + shift) / sr))}}
\item{\code{"step"}}{the digitised step-response of the filter, i.e. \code{stepfun((0:len + shift) / sr)}}
\item{\code{"acf"}}{the discrete autocorrelation, i.e. \code{acfun(0:len / sr)}}
}
}

\seealso{\code{\link{filter}}}

\author{This function is a modified and extended version of the \code{\link{dfilter}} function in the \code{\link{stepR}} package written by Thomas Hotz. New code is written by Florian Pein and Inder Tecuapetla-Gómez.}

\examples{
# the filter used for the gramicidin A recordings given by gramA
filter <- lowpassFilter(type = "bessel", param = list(pole = 4L, cutoff = 1e3 / 1e4),
                        sr = 1e4)

# filter kernel, truncated version
plot(filter$kernfun, xlim = c(0, 20 / filter$sr))
t <- seq(0, 20 / filter$sr, 0.01 / filter$sr)
# truncated version looks very similar
lines(t, filter$truncatedKernfun(t), col = "red")

# filter$len (== 11) is chosen such that filter$acf < 1e-3 for it and all larger lags
plot(filter$acfun, xlim = c(0, 20 / filter$sr), ylim = c(-0.003, 0.003))
abline(h = 0.001, lty = "22")
abline(h = -0.001, lty = "22")

abline(v = (filter$len - 1L) / filter$sr, col = "grey")
abline(v = filter$len / filter$sr, col = "red")
\donttest{
## zoom into a single jump of the idealisation
## we suggest to do this for every new measurement setup once 
## to control whether the correct filter is assumed
# idealisation by JULES (might take some time if not called somewhere before,
# please see its documentation for more details)
idealisation <- jules(gramA, filter = filter, startTime = 9, messages = 100)

## zoom into a single jump
plot(9 + seq(along = gramA) / filter$sr, gramA, pch = 16, col = "grey30",
     ylim = c(20, 50), xlim = c(9.6476, 9.6496), ylab = "Conductance in pS",
     xlab = "Time in s")

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
# the needed Monte-Carlo simulation depends on the number of observations and the filter
# hence a new simulation is required (if called for the first time)
idealisationWrong <- jules(gramA, filter = wrongFilter, startTime = 9, messages = 100)

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
}

# filter with sr == 1
filter <- lowpassFilter(type = "bessel", param = list(pole = 4L, cutoff = 1e3 / 1e4))

# filter kernel and its truncated version
plot(filter$kernfun, xlim = c(0, 20 / filter$sr))
t <- seq(0, 20 / filter$sr, 0.01 / filter$sr)
# truncated version looks very similar
lines(t, filter$truncatedKernfun(t), col = "red")
# digitised filter
points((0:filter$len + 0.5) / filter$sr, filter$kern, col = "red", pch = 16)

# without a shift
filter <- lowpassFilter(type = "bessel", param = list(pole = 4L, cutoff = 1e3 / 1e4),
                        shift = 0)
# filter$kern starts with zero
points(0:filter$len / filter$sr, filter$kern, col = "blue", pch = 16)

# much shorter filter
filter <- lowpassFilter(type = "bessel", param = list(pole = 4L, cutoff = 1e3 / 1e4),
                        len = 4L)
points((0:filter$len + 0.5) / filter$sr, filter$kern, col = "darkgreen", pch = 16)
}

\keyword{ts}
