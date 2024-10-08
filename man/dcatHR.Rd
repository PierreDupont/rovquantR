% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dcatHR.R
\name{dcatHR}
\alias{dcatHR}
\alias{rcatHR}
\title{Density and random generation of a categorical distribution describing OPSCR
 state transitions with one alive and two dead states}
\usage{
dcatHR(x, z, gamma, mhH, mhW, log = 0)

rcatHR(n, z, gamma, mhH, mhW)
}
\arguments{
\item{x}{Scalar, individual arrival state \code{z_{i,t+1}}.}

\item{z}{Scalar, individual departure state \code{z_{i,t}}.}

\item{gamma}{scalar, probability to transition from state 1 to 2.}

\item{mhH}{scalar, probability to transition from state 2 to 3.}

\item{mhW}{scalar, probability to transition from state 2 to 4.}

\item{log}{Logical argument, specifying whether to return the log-probability of the distribution.
\code{dcatHR} gives the (log) probability density of \code{x}. 
\code{rcatHR} gives a randomly generated individual states conditional on the departure state \code{z}.}

\item{n}{Integer specifying the number of realizations to generate.  Only n = 1 is supported.}
}
\description{
The \code{dcatHR} distribution is a NIMBLE custom distribution which can be used 
to model and simulate individual state transitions. This function can be used 
to model transitions from one alive and two dead states, corresponding to 
two competing mortality causes.
}
\details{
The possible transitions are:
If \code{z_{i,t} = 1}, individual i can be recruited (transition to state 2) with probability \code{gamma}, so that \code{z_{i,t+1} ~ dcat(1-gamma, gamma, 0, 0, 0)} where gamma represents the probability of an unborn individual to be recruited.
If \code{z_{i,t} = 2}, individual i can die from one cause of mortality (e.g. culling) and transition to \code{z_{i,t+1} = 3} with log-hazard rate \code{mhH}, or die from the second mortality cause with log-hazard rate \code{mhW} and transition to \code{z_{i,t+1} = 4}. If the individual does not die it can survive and remain in state 2
Individuals in dead states \code{z_{i,t} = 3 or 4}, transition to \code{z_{i,t+1} = 4}, the absorbing state, with probability 1.
}
\examples{
# Use the distribution in R
z <- 2
gamma <- 0.2
mhH <- 0.4
mhW <- 0.1

## Sample random state
zPlusOne <- rcatHR( z  = z,
                    gamma = gamma,
                    mhH = mhH,
                    mhW = mhW)
                                 
## Calculate probability of transition  
dcatHR( x = zPlusOne,
       z = z,
       gamma = gamma,
       mhH = mhH,
       mhW = mhW,
       log = TRUE)

}
\author{
Pierre Dupont
}
