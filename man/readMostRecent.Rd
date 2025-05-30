% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readMostRecent.R
\name{readMostRecent}
\alias{readMostRecent}
\alias{readMostRecent.csv}
\alias{readMostRecent.excel}
\alias{readMostRecent.RData}
\alias{getMostRecent}
\alias{writeMostRecent.csv}
\title{Read most recent version of a file}
\usage{
readMostRecent(
  path,
  extension = ".csv",
  pattern = NULL,
  returnDate = FALSE,
  sep = ",",
  dec = ".",
  ...
)

readMostRecent.csv(path, returnDate = F, fileEncoding = NULL, ...)

readMostRecent.excel(path, returnDate = F, ...)

readMostRecent.RData(path, pattern = ".RData", returnDate = F, ...)

getMostRecent(path, pattern)

writeMostRecent.csv(x, file)
}
\arguments{
\item{path}{Vector of individual detection frequencies.}

\item{extension}{Integer specifying the number of realizations to generate. 
Only n = 1 is supported.}

\item{pattern}{Vector of indices of traps where the detections in \emph{x}
were recorded; from the \emph{detIndices} object returned by the \code{\link{getSparseY}} function.}

\item{returnDate}{Number of traps with at least one detection recorded in 
\emph{x}; from the \emph{detNums} object returned by the \code{\link{getSparseY}} function.}

\item{sep}{Vector of the number of trials (zero or more) for each trap (\emph{trapCoords}).}

\item{dec}{Baseline detection probability (scalar) used in the half-normal 
detection function. For trap-specific baseline detection probabilities use 
argument \emph{p0Traps} (vector) instead.}

\item{...}{additional optional parameters.}

\item{fileEncoding}{Vector of baseline detection probabilities for each trap 
used in the half-normal detection function. When \emph{p0Traps} is used, 
\emph{p0} should not be provided.}

\item{x}{Scale parameter of the half-normal detection function.}

\item{file}{Individual activity center x- and y-coordinates scaled to the 
habitat (see (\code{\link{scaleCoordsToHabitatGrid}}).}
}
\value{
The data loaded
}
\description{
The \code{readMostRecent} functions provide a suite of utility functions to 
identify and load the most recent file matching a set of characteristics.
}
\author{
Pierre Dupont
}
