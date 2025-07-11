% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cleanRovBaseData.R
\name{cleanRovbaseData}
\alias{cleanRovbaseData}
\title{Data set clean-up.}
\usage{
cleanRovbaseData(
  data.dir = "./Data",
  working.dir = NULL,
  species = c("bear", "wolf", "wolverine"),
  years = NULL,
  two.sex = TRUE,
  sampling.months = NULL,
  rename.list = NULL,
  legal.dead = NULL,
  print.report = TRUE,
  Rmd.template = NULL,
  output.dir = NULL,
  overwrite = FALSE
)
}
\arguments{
\item{data.dir}{the \code{path} pointing to the directory containing the raw data from Rovbase.}

\item{working.dir}{the \code{path} pointing to the working directory. By default, the cleaned data will be stored in a subfolder of this working directory called 'data'.}

\item{species}{A \code{character} string with the name of the focal species
("bear", "wolf", or "wolverine").}

\item{years}{A \code{numeric} vector containing the years of interest. 
Only data for those years will be cleaned and returned.}

\item{two.sex}{A \code{logical} determining whether the analysis will be done by sex (two.sex = T) or both together (two.sex = F).}

\item{sampling.months}{(Optional) A \code{list} containing the sampling period months. If the sampling period overlaps two calendar years, the list should contain one element per year (e.g. samplingMonths <- list(c(11,12), c(1,2,3,4))) for a sampling period extending from November to April of the following year.}

\item{rename.list}{(Optional) A named \code{character} vector used to rename columns in the raw Rovbase files.}

\item{print.report}{A \code{logical} denoting whether to print out a \code{.html} report summarizing the cleaning process or not.}

\item{Rmd.template}{(Optional) The \code{path} to the \code{.rmd} template to be used for cleaning the data. By default, the \code{.rmd} template provided with the \code{rovquantR} package is used.}

\item{output.dir}{(Optional) the \code{path} pointing to the directory where the \code{.html} report will be printed.
By default, the \code{.html} report describing the content of the clean data will 
be placed in a subfolder of the working directory (\code{working.dir}) called 'reports'.}

\item{overwrite}{A \code{logical} (default = FALSE) to force overwriting of previously existing clean data.
If FALSE, the function checks for any pre-existing clean data files and ask whether to overwrite it or not.}
}
\value{
This function returns:
\enumerate{
\item A \code{.RData} file with the clean NGS and dead recovery data objects 
for the species and period specified. The clean data file is saved as an \code{.RData} 
file named using the species name and the date of extraction of the raw Rovbase data 
to facilitate replicability (e.g. 'CleanData_bear_2024-08-10.RData').
\item A \code{.html} report summarizing the data cleaning process. 
The \code{.RData} report is using the same naming convention as the clean \code{.RData} (e.g. 'CleanData_bear_2024-08-10.html').
\item Additional \code{.png} images and summary \code{.csv} tables that can be reused somewhere else.
}
}
\description{
\code{cleanRovbaseData} identifies and loads the most recent RovBase data available 
for the specified species in the specified \code{data.dir} location, and conducts a set of data cleaning steps that include:
\itemize{
 \item{removing un-identified samples}
 \item{checking sex-assignment}
 \item{removing samples flagged as unusable by RovData/Mike}
}
 
Additionally, it can produce a \code{html} report describing the content of the data in terms of number of samples, individuals, etc...
}
\author{
Pierre Dupont
}
