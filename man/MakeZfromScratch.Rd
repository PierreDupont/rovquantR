% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MakeZfromScratch.R
\name{MakeZfromScratch}
\alias{MakeZfromScratch}
\title{Create initial values for individual states}
\usage{
MakeZfromScratch(data.alive, data.dead = NULL, samplingMonths = 1:12)
}
\arguments{
\item{data.alive}{A \code{sf dataframe}. object containing individual detections.}

\item{data.dead}{A \code{sf dataframe}.}

\item{samplingMonths}{A \code{numeric}.}
}
\value{
A \code{Matrix} containing the set of initial states z.
}
\description{
\code{MakeZfromScratch} creates a sequence of states for each individual
that is compatile with its detection history, i.e. "unborn" before the first detection,
"alive" between the first and last detection and "dead" after that.
 the dimensions of an object along
}
\examples{
\dontrun{MakeZfromScratch()}

}
