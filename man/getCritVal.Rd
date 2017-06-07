\name{getCritVal}
\alias{getCritVal}

\title{Critical value}
\description{Computes the critical value \code{q} in (\cite{Pein et al.}, 2017, (7)) based on a Monte-Carlo simulation such that (7) is a level \code{alpha} test. Since a Monte-Carlo simulation lasts potentially much longer (up to several hours or days if the number of observations is in the millions) than the main calculations, the simulations are by default saved in the workspace and on the file system such that a second call that require the same Monte-Carlo simulation will be much faster. For more details, in particular to which arguments the Monte-Carlo simulations are specific, see Section \cite{Simulating, saving and loading of Monte-Carlo simulations} below. Progress of a Monte-Carlo simulation can be reported by the argument \code{messages} and the saving can be controlled by the argument \code{option}.}

\usage{
getCritVal(n, alpha = 0.05, filter, r = 1e4, nq = n, options = NULL, stat = NULL,
           messages = NULL)
}

\arguments{
\item{n}{a positive integer giving the number of observations}
\item{alpha}{a probability, i.e. a single numeric between 0 and 1, giving the significance level. Its choice is a trade-off between data fit and parsimony of the estimator. In other words, this argument balances the risks of missing changes and detecting additional artefacts. For more details on this choice see (Frick et al., 2014, section 4) and (Pein et al., 2016, section 3.4)}
\item{filter}{an object of class \code{\link{lowpassFilter}} giving the used analogue lowpass filter}
\item{r}{a positive integer giving the required number of Monte-Carlo simulations if they will be simulated or loaded from the workspace or the file system}
\item{nq}{a positive integer larger than or equal to \code{n} giving the (increased) number of observations for the Monte-Carlo simulation. See Section \cite{Simulating, saving and loading of Monte-Carlo simulations} for more details}
\item{options}{a \code{\link{list}} specifying how Monte-Carlo simulations will be simulated, saved and loaded. For more details see Section \cite{Simulating, saving and loading of Monte-Carlo simulations}}
\item{stat}{an object of class \code{"MCSimulationVector"} or \code{"MCSimulationMaximum"}, usually computed by \code{\link{monteCarloSimulation}}. Has to be simulated for at least the given number of observations \code{n} and for the given \code{filter}. If missing it will be loaded and if not found simulated accordingly to the given \code{options}. For more details see Section \cite{Simulating, saving and loading of Monte-Carlo simulations}}
\item{messages}{a positive integer or \code{NULL}, in each \code{messages} iteration a message will be given in order to show the progress of the simulation, if \code{NULL} no message will be given}
}

\value{
A single numeric giving the critical value q in (\cite{Pein et al.}, 2017, (7)).
}

\section{Simulating, saving and loading of Monte-Carlo simulations}{
Since a Monte-Carlo simulation lasts potentially much longer (up to several hours or days if the number of observations is in the millions) than the main calculations, this function offers multiple possibilities to save and load the simulations. The simulation, saving and loading can be controlled by the argument \code{option}. This argument has to be a \code{\link{list}} or \code{NULL}. For the \code{\link{list}} the following named entries are allowed: \code{"simulation"}, \code{"save"}, \code{"load"}, \code{"envir"} and \code{"dirs"}. All missing entries will be set to their default option.\cr
Each Monte-Carlo simulation is specific to the number of observations and the used filter. Monte-Carlo simulations can also be performed for a (slightly) larger number of observations \eqn{n_q}{nq} given in the argument \code{nq}, which avoids extensive resimulations for only a little bit varying number of observations at price of a (slightly) smaller detection power. We recommend to not use a \code{nq} more than two times larger than the number of observations \code{n}.\cr
Objects of the following types can be simulated, saved and loaded:
\itemize{
\item \code{"vector"}: an object of class \code{"MCSimulationMaximum"} for \code{n} observations, i.e. a numeric vector of length \code{r}
\item \code{"vectorIncreased"}: an object of class \code{"MCSimulationMaximum"} for \code{nq} observations, i.e. a numeric vector of length \code{r} 
\item \code{"matrix"}: an object of class \code{"MCSimulationVector"} for \code{n} observations, i.e. a matrix of dimensions \code{as.integer(log2(n)) + 1L} and \code{r}
\item \code{"matrixIncreased"}: an object of class \code{"MCSimulationVector"} for \code{nq} observations, i.e. a matrix of \code{as.integer(log2(n)) + 1L} and \code{r}
}
Objects of class \code{"MCSimulationVector"} and objects of class \code{"MCSimulationMaximum"} lead to the same result (if the number of observations is the same), but an object of class \code{"MCSimulationVector"} requires much more storage space and has slightly larger saving and loading times. However, simulations of type \code{"vectorIncreased"}, i.e. objects of class \code{"MCSimulationMaximum"} with \code{nq} observations, have to be resimulated if \code{as.integer(log2(n1)) != as.integer(log2(n2))} when the saved simulation was computed with \code{n == n1} and the simulation now is required for \code{n == n2} and \code{nq >= n1} and \code{nq >= n2}. All in all, if all data sets in the analysis have the same number of observations simulations of type \code{"vector"} are recommended. If they have a slightly different number of observations it is recommend to set \code{nq} to the largest number and to use simulations for an increased number of observations: If \code{as.integer(log2(n))} is the same for all data sets type \code{"vectorIncreased"} is recommend, if they differ type \code{"matrixIncreased"} avoids a resimulation at the price of a larger object to be stored and loaded.\cr
The simulations can either be saved in the workspace in the variable \code{critValStepRTab} or persistently on the file system for which the package \code{\link[R.cache]{R.cache}} is used. Loading from the workspace is faster, but either the user has to save the workspace manually or in a new session simulations have to be performed again. Moreover, storing in and loading from variables and \link[=readRDS]{RDS} files is supported.

\subsection{options$envir and options$dirs}{
For loading from / saving in the workspace the variable \code{critValStepRTab} in the \code{\link{environment}} \code{options$envir} will be looked for and if missing in case of saving also created there. Moreover, the variable(s) specified in \code{options$save$variable} (explained in the Subsection \cite{Saving: options$save}) will be assigned to this \code{\link{environment}}. By default the \link[=.GlobalEnv]{global environment} \code{.GlobalEnv} is used, i.e. \code{options$envir == .GlobalEnv}.\cr
For loading from / saving on the file system \code{\link[R.cache]{loadCache}(key = keyList, dirs = options$dirs)} and \code{\link[R.cache]{saveCache}(stat, key = attr(stat, "keyList"), dirs = options$dirs)} are called, respectively. In other words, \code{options$dirs} has to be a \code{\link{character}} \code{\link{vector}} constituting the path to the cache subdirectory relative to the cache root directory as returned by \code{\link[R.cache]{getCacheRootPath}}(). If \code{options$dirs == ""}, the path will be the cache root path. By default the subdirectory \code{"stepR"} is used, i.e. \code{options$dirs == "stepR"}. Missing directories will be created.
}

\subsection{Simulation: options$simulation}{
Whenever Monte-Carlo simulations have to be performed, i.e. when \code{stat == NULL} and the required Monte-Carlo simulation could not be loaded, the type specified in \code{options$simulation} will be simulated by \code{\link{monteCarloSimulation}}. In other words, \code{options$simulation} must be a single string of the following: \code{"vector"}, \code{"vectorIncreased"}, \code{"matrix"} or \code{"matrixIncreased"}. By default (\code{options$simulation == NULL}), an object of class \code{"MCSimulationVector"} for \code{nq} observations will be simulated, i.e. \code{options$simulation} \code{== "matrixIncreased"}. For this choice please recall the explanations regarding computation time and flexibility at the beginning of this section.
}

\subsection{Loading: options$load}{
Loading of the simulations can be controlled by the entry \code{options$load} which itself has to be a \code{\link{list}} with possible entries: \code{"RDSfile"}, \code{"workspace"}, \code{"package"} and \code{"fileSystem"}. Missing entries disable the loading from this option.
Whenever a Monte-Carlo simulation is required, i.e. when the variable \code{q} is not given, it will be searched for at the following places in the given order until found:
\enumerate{
\item in the variable \code{stat},
\item in \code{options$load$RDSfile} as an \link[=readRDS]{RDS} file, i.e. the simulation will be loaded by \preformatted{readRDS(options$load$RDSfile).} In other words, \code{options$load$RDSfile} has to be a \code{\link{connection}} or the name of the file where the \R object is read from,
\item in the workspace or on the file system in the following order: \code{"vector"}, \code{"matrix"}, \code{"vectorIncreased"} and finally of \code{"matrixIncreased"}. For each option it will first be looked in the workspace and then on the file system. All searches can be disabled by not specifying the corresponding string in \code{options$load$workspace} and \code{options$load$fileSystem}. In other words, \code{options$load$workspace} and \code{options$load$fileSystem} have to be vectors of strings containing none, some or all of \code{"vector"}, \code{"matrix"}, \code{"vectorIncreased"} and \code{"matrixIncreased"},
\item if all other options fail a Monte-Carlo simulation will be performed.
}
By default (if \code{options$load} is missing / \code{NULL}) no \link[=readRDS]{RDS} file is specified and all other options are enabled, i.e. \preformatted{
options$load <- list(workspace = c("vector", "vectorIncreased",
                                   "matrix", "matrixIncreased"),
                     fileSystem = c("vector", "vectorIncreased",
                                    "matrix", "matrixIncreased"),
                     RDSfile = NULL).}
}

\subsection{Saving: options$save}{
Saving of the simulations can be controlled by the entry \code{options$save} which itself has to be a \code{\link{list}} with possible entries: \code{"workspace"}, \code{"fileSystem"}, \code{"RDSfile"} and \code{"variable"}. Missing entries disable the saving in this option.\cr
All available simulations, no matter whether they are given by \code{stat}, loaded, simulated or in case of \code{"vector"} and \code{"vectorIncreased"} computed from \code{"matrix"} and \code{"matrixIncreased"}, respectively, will be saved in all options for which the corresponding type is specified. Here we say a simulation is of type \code{"vectorIncreased"} or \code{"matrixIncreased"} if the simulation is not performed for \code{n} observations. More specifically, a simulation will be saved:
\enumerate{
\item in the workspace or on the file system if the corresponding string is contained in \code{options$save$workspace} and \code{options$save$fileSystem}, respectively. In other words, \code{options$save$workspace} and \code{options$save$fileSystem} have to be vectors of strings containing none, some or all of \code{"vector"}, \code{"matrix"}, \code{"vectorIncreased"} and \code{"matrixIncreased"},
\item in a variable named by \code{options$save$variable} in the \code{\link{environment}} \code{options$envir}. Hence, \code{options$save$variable} has to be a vector of one or two containing variable names (character vectors). If \code{options$save$variable} is of length two a simulation of type \code{"vector"} or \code{"vectorIncreased"} (only one can occur at one function call) will be saved in \code{options$save$variable[1]} and \code{"matrix"} or \code{"matrixIncreased"} (only one can occur at one function call) will be saved in \code{options$save$variable[2]}. If \code{options$save$variable} is of length one both will be saved in \code{options$save$variable} which means if both occur at the same call only \code{"vector"} or \code{"vectorIncreased"} will be saved. Each saving can be disabled by not specifying \code{options$save$variable} or by passing \code{""} to the corresponding entry of \code{options$save$variable}.
}
By default (if \code{options$save} is missing) \code{"vector"} and \code{"vectorIncreased"} will be saved in the workspace and \code{"matrixIncreased"} on the file system, i.e. \preformatted{
options$save <- list(workspace = c("vector", "vectorIncreased"),
                     fileSystem = c("matrix", "matrixIncreased"),
                     RDSfile = NULL, variable = NULL).}
Simulations can be removed from the workspace by removing the variable \code{critValStepRTab}, i.e. by calling \code{remove(critValStepRTab, envir = envir)}, with \code{envir} the used environment, and from the file system by deleting the corresponding subfolder, i.e. by calling \preformatted{unlink(file.path(R.cache::getCacheRootPath(), dirs), recursive = TRUE),} with \code{dirs} the corresponding subdirectory.
}
}

\seealso{\code{\link{jules}}, \code{\link{lowpassFilter}}, \code{\link{stepDetection}}}

\encoding{utf8}
\references{
Pein, F., Tecuapetla-Gómez, I., \enc{Schütte}{Schuette}, O. Steinem, C., and Munk, A. (2017) Fully-automatic multiresolution idealization for filtered ion channel recordings: flickering event detection. In preparation.

Frick, K., Munk, A. and Sieling, H. (2014) Multiscale change-point inference. With discussion and rejoinder by the authors. \emph{Journal of the Royal Statistical Society, Series B} \bold{76}(3), 495--580.

Pein, F., Sieling, H. and Munk, A. (2016) Heterogeneous change point inference. \emph{Journal of the Royal Statistical Society, Series B}, early view.
} 

\examples{
# the for the recording of the gramA data set used filter
filter <- lowpassFilter(type = "bessel", param = list(pole = 4L, cutoff = 1e3 / 1e4),
                        sr = 1e4)
\donttest{
# critical value for jules or stepDetection
# this call requires a Monte-Carlo simulation at the first time
# and therefore might last a few minutes,
# progress of the Monte-Carlo simulation is reported
q <- getCritVal(length(gramA), filter = filter, messages = 100)

# this second call should be much faster
# as the previous Monte-Carlo simulation will be loaded
getCritVal(length(gramA), filter = filter)

# much larger significance level alpha for a larger detection power,
# but also with the risk of detecting additional artefacts
getCritVal(length(gramA), filter = filter, alpha = 0.9)

# medium significance level alpha for a tradeoff between detection power
# and the risk to detect additional artefacts
getCritVal(length(gramA), filter = filter, alpha = 0.5)}

# critical values depend on the number of observations and on the filter
# also a new Monte-Carlo simulation is required
getCritVal(100, filter = filter, messages = 500)
\donttest{
otherFilter <- lowpassFilter(type = "bessel",
                             param = list(pole = 6L, cutoff = 0.2),
                             sr = 1e4)
getCritVal(100, filter = otherFilter, messages = 500)

# simulation for a larger number of oberservations can be used (nq = 100)
# does not require a new simulation as the simulation from above will be used
# (if the previous call was executed first)
getCritVal(90, filter = filter, nq = 100)

# simulation of type "vectorIncreased" for n1 observations can only be reused 
# for n2 observations if as.integer(log2(n1)) == as.integer(log2(n2))
# no simulation is required, since a simulation of type "matrixIncreased"
# will be loaded from the fileSystem
# this call also saved a simulation of type "vectorIncreased" in the workspace
getCritVal(30, filter = filter, nq = 100)
# here a new simulation is required
# (if no appropriate simulation is saved from a call outside of this file)
getCritVal(10, filter = filter, nq = 100, messages = 500,
           options = list(load = list(workspace = c("vector", "vectorIncreased")))) 

# the above calls saved and (attempted to) load Monte-Carlo simulations
# in the following call the simulations will neither be saved nor loaded
# to save some time the number of iterations is reduced to r = 1e3
# hence the critical value is computed with less precision
# In general, r = 1e3 is enough for a first impression
# for a detailed analysis r = 1e4 is suggested
getCritVal(100, filter = filter, messages = 100, r = 1e3,
           options = list(load = list(), save = list()))

# simulations will only be saved in and loaded from the workspace,
# but not on the file system
getCritVal(100, filter = filter, messages = 100, r = 1e3,
           options = list(load = list(workspace = c("vector", "vectorIncreased")),
                          save = list(workspace = c("vector", "vectorIncreased"))))

# explicit Monte-Carlo simulations, not recommended
stat <- stepR::monteCarloSimulation(n = 100, , family = "mDependentPS",
                                    filter = filter, output = "maximum",
                                    r = 1e3, messages = 100)
getCritVal(100, filter = filter, stat = stat)
}}

\keyword{nonparametric}
