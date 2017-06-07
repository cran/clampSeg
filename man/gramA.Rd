\name{gramA}
\docType{data}
\alias{gramA}
\alias{gramicidin A}
\alias{gramicidin}
\alias{gramicidinA}

\title{Patch clamp recording of gramicidin A}

\description{
3 second part of a patch clamp recording of gramicidin A with solvent-free lipid bilayers using the Port-a-Patch measured in the Steinem lab (Institute of Organic and Biomolecular Chemistry, University of Goettingen). All rights reserved by them. The recorded data points are a conducatance trace in pico Siemens and were recorded at a sampling rate of 10 kHz using a 1 kHz 4-pole Bessel filter. More details of the recording can be found in (\cite{Pein et al.}, 2017, section V A) and a plot in the examples or in (\cite{Pein et al.}, 2017, figure 1 lower panel).
}

\usage{gramA}

\format{A \code{\link{numeric}} vector containing 30,000 values.}

\encoding{utf8}
\references{
Pein, F., Tecuapetla-Gómez, I., \enc{Schütte}{Schuette}, O. Steinem, C., and Munk, A. (2017) Fully-automatic multiresolution idealization for filtered ion channel recordings: flickering event detection. In preparation.
}

\examples{
# the recorded data points
gramA

# the used filter
filter <- lowpassFilter(type = "bessel", param = list(pole = 4L, cutoff = 1e3 / 1e4),
                        sr = 1e4)

# the corresponding time points
time <- 9 + seq(along = gramA) / filter$sr

# plot of the data as in (Pein et al., 2017, figure 1 lower panel)
plot(time, gramA, pch = ".", col = "grey30", ylim = c(20, 50),
     ylab = "Conductance in pS", xlab = "Time in s")
}

\keyword{datasets}
