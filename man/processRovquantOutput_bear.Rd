% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processRovquantOutput_bear.R
\name{processRovquantOutput_bear}
\alias{processRovquantOutput_bear}
\title{RovQuant OPSCR bear output processing}
\usage{
processRovquantOutput_bear(
  data.dir = getwd(),
  working.dir = NULL,
  nburnin = 0,
  niter = 100,
  extraction.res = 5000,
  overwrite = FALSE
)
}
\arguments{
\item{data.dir}{A \code{path}}

\item{working.dir}{A \code{path}}

\item{nburnin}{An \code{integer} denoting the number of iterations to be removed from each MCMC as burnin.}

\item{niter}{An \code{integer} denoting the number of MCMC iterations to be used for density extraction.}

\item{extraction.res}{A \code{integer} denoting the raster resolution for density extraction.}
}
\value{
A \code{.RData} file with the clean NGS and dead recovery data objects
for the species and period specified.
A \code{html} report summarizing the data cleaning process
Additional \code{.png} images that can be reused somewhere else.
}
\description{
\code{processRovquantOutput_bear} calls a custom Rmarkdown template that combines 
and processes MCMC outputs from NIMBLE models and produces figures,
tables and rasters of interest (e.g. population density maps)
}
\author{
Pierre Dupont
}
