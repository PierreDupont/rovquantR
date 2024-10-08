% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/translateForeignCharacters.R
\name{translateForeignCharacters}
\alias{translateForeignCharacters}
\title{Scandinavian character translation}
\usage{
translateForeignCharacters(data, dir.translation = NULL)
}
\arguments{
\item{data}{a character vector where matches are sought, or an object which can be coerced by as.character to a character vector.}

\item{dir.translation}{(optional) the path to the file containing special characters and their "non-special" translation.}
}
\value{
the same object than \code{data} with all Scandinavian characters translated.
}
\description{
Replaces Scandinavian characters for easier processing.
}
\examples{
translateForeignCharacters("bjørn")

}
\author{
Richard Bischof and Pierre Dupont
}
