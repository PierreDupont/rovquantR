#' Density and random generation of a categorical distribution describing OPSCR
#'  state transitions with one alive and two dead states  
#'
#' The \code{dcatHR} distribution is a NIMBLE custom distribution which can be used 
#' to model and simulate individual state transitions. This function can be used 
#' to model transitions from one alive and two dead states, corresponding to 
#' two competing mortality causes.
#' 
#' The possible transitions are:
#' \itemize{
#' \item If z_{i,t} = 1, individual i can be recruited (transition to state 2) with probability \code{gamma}, so that z_{i,t+1} ~ dcat(1-gamma, gamma, 0, 0, 0) where gamma represents the probability of an unborn individual to be recruited.
#' \item If z_{i,t} = 2, individual i can die from one cause of mortality (e.g. culling) and transition to z_{i,t+1}=3 with log-hazard rate ´\code{mhH}, or die from the second mortality cause with log-hazard rate \code{mhW} and transition to z_{i,t+1}=4. If the individual does not die it can survive and remain in state 2
#' \item Individuals in dead states (z_{i,t} = 3 or 4) transition to z_{i,t+1} = 4, the absorbing state, with probability 1.
#' } 
#' 
#' @name dcatHR
#'
#' @param x Scalar, individual arrival state (corresponding to z_{i,t+1}).
#' @param n Integer specifying the number of realizations to generate.  Only n = 1 is supported.
#' @param z Scalar, individual departure state (i.e. z_{i,t}).
#' @param gamma scalar, probability to transition from state 1 to 2.
#' @param mhH scalar, probability to transition from state 2 to 3.
#' @param mhW scalar, probability to transition from state 2 to 4.
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#
#' \code{dcatHR} gives the (log) probability density of \code{x}. 
#' \code{rcatHR} gives a randomly generated individual states conditional on the departure state \code{z}.  
#' 
#' @author Pierre Dupont
#'
#' @import nimble
#'
#' @examples
#' # Use the distribution in R
#' 
#' z <- 2
#' gamma <- 0.2
#' mhH <- 0.4
#' mhW <- 0.1
#' 
#' 
#' ## Sample random state
#' zPlusOne <- rcatHR( z  = z,
#'                     gamma = gamma,
#'                     mhH = mhH,
#'                     mhW = mhW)
#'                                  
#' ## Calculate probability of transition  
#'dcatHR( x = zPlusOne,
#'        z = z,
#'        gamma = gamma,
#'        mhH = mhH,
#'        mhW = mhW,
#'        log = TRUE)
#' 
#' @export
NULL

#' @rdname dcatHR
#' @export
#### 1.Density function ####
dcatHR <- nimbleFunction(run = function( 
    x = double(0),
    z = double(0),
    gamma = double(0),
    mhH = double(0),
    mhW = double(0),
    log = integer(0, default = 0)
){
  # Return type declaration
  returnType(double(0))
  
  if(z == 1){
    logLikelihood <- dcat(x, prob = c(1 - gamma, gamma), log=1)
    if(log == 1){return(logLikelihood)}else{return(exp(logLikelihood))}
  }
  
  if(z == 2){
    mhH1 <- exp(mhH)
    mhW1 <- exp(mhW)
    
    h <- (1-exp(-(mhH1+mhW1)))* (mhH1/(mhH1+mhW1))
    w <- (1-exp(-(mhH1+mhW1)))* (mhW1/(mhH1+mhW1))
    phi <- 1-h-w
    
    logLikelihood <- dcat(x, prob = c(0, phi, h, w), log=1)
    if(log == 1){return(logLikelihood)}else{return(exp(logLikelihood))}
  }
  
  if(z == 3 | z == 4){
    logLikelihood <- dcat(x, prob = c(0, 0, 0, 1), log=1)
    if(log == 1){return(logLikelihood)}else{return(exp(logLikelihood))}
  }
})


#### 2.Sampling function ####
rcatHR <- nimbleFunction(run = function( 
    n = integer(0),
    z = double(0),
    gamma = double(0),
    mhH = double(0),
    mhW = double(0)
){
  # Return type declaration
  returnType(double(0))

  if(z == 1){
    state <- rcat(1, prob = c(1 - gamma, gamma))
    return(state)
  }
  
  if(z == 2){
    mhH1 <- exp(mhH)
    mhW1 <- exp(mhW)
    
    h <- (1-exp(-(mhH1+mhW1)))* (mhH1/(mhH1+mhW1))
    w <- (1-exp(-(mhH1+mhW1)))* (mhW1/(mhH1+mhW1))
    phi <- 1-h-w
    
    state <- rcat(1, prob = c(0, phi, h, w))
    return(state)
  }
  
  if(z == 3 | z == 4){
    state <- 4
    return(state)
  }
  
})

#### 3.Registration ####
registerDistributions(list(
  dcatHR = list(
    BUGSdist = "dcatHR(z, gamma, mhH, mhW)",
    types = c( "value = double(0)",
               "z = double(0)",
               "gamma = double(0)",
               "mhH = double(0)",
               "mhW = double(0)"
    ),
    pqAvail = FALSE,
    discrete = TRUE
  )))

