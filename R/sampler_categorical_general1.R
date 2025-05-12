#' MCMC Sampling Algorithms
#'
#' Details of the MCMC sampling algorithms provided with the NIMBLE MCMC engine; HMC samplers are in the \code{nimbleHMC} package and particle filter samplers are in the \code{nimbleSMC} package.
#'
#'
#' @param model (uncompiled) model on which the MCMC is to be run
#' @param mvSaved \code{modelValues} object to be used to store MCMC samples
#' @param target node(s) on which the sampler will be used
#' @param control named list that controls the precise behavior of the sampler, with elements specific to \code{samplertype}.  The default values for control list are specified in the setup code of each sampling algorithm.  Descriptions of each sampling algorithm, and the possible customizations for each sampler (using the \code{control} argument) appear below.
#'
#' @name sampler_categorical_general1
#'
#' @return After an MCMC algorithm has been configured and built, the value of the proposal standard deviation of a RW_block sampler can be modified using the setScale method of the sampler object.  This use the scalar argument to will modify the current value of the proposal standard deviation, as well as modifying the initial (pre-adaptation) value which the proposal standard deviation is reset to, at the onset of a new MCMC chain.
#'
#' Operating analogous to the setScale method, the RW_block sampler also has a setPropCov method.  This method accepts a single matrix-valued argument, which will modify both the current and initial (used at the onset of a new MCMC chain) values of the multivariate normal proposal covariance.
#'
#' Note that modifying elements of the control list may greatly affect the performance of this sampler. In particular, the sampler can take a long time to find a good proposal covariance when the elements being sampled are not on the same scale. We recommend providing an informed value for \code{propCov} in this case (possibly simply a diagonal matrix that approximates the relative scales), as well as possibly providing a value of \code{scale} that errs on the side of being too small. You may also consider decreasing \code{adaptFactorExponent} and/or \code{adaptInterval}, as doing so has greatly improved performance in some cases. 
#'
#' @section General categorical sampler:
#'
#' The sampler_categorical_general1 sampler is used to sample individual states in demographic models with transition probabilities depending on hazard rates.
#' This sampler is not assigned as a default sampler by \code{configureMCMC} and so can only be used if manually added to an MCMC configuration.
#' 
#' 
#' @seealso \code{\link{configureMCMC}} \code{\link{addSampler}} \code{\link{buildMCMC}} \code{\link{runMCMC}}
#'
#' @author Daniel Turek
#'
#' @import nimble 
#' 
#' @rdname sampler_categorical_general1
#' @export
sampler_categorical_general1 <- nimbleFunction(
  name = 'sampler_categorical1',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## node list generation
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes  <- model$getDependencies(target)
    calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    ## numeric value generation
    k <- control[["numCategories"]]
    if(is.null(k)) stop("Must provide control argument numCategories")
    probs <- numeric(k)
    logProbs <- numeric(k)
    ## checks
    if(length(targetAsScalar) > 1)  stop('cannot use categorical sampler on more than one target node')
   #[CM] if(model$getDistribution(target) != 'dCatZwh') stop('can only use categorical sampler on node with dcat distribution')
  },
  run = function() {
    currentValue <- model[[target]]
    logProbs[currentValue] <<- getLogProb(model, calcNodes)
    for(i in 1:k) {
      if(i != currentValue) {
        model[[target]] <<- i
        logProbPrior <- calculate(model, target)
        if(logProbPrior == -Inf) {
          logProbs[i] <<- -Inf
        } else {
          if(is.nan(logProbPrior)) {
            logProbs[i] <<- -Inf
          } else {
            logProbs[i] <<- logProbPrior + calculate(model, calcNodesNoSelf)
            if(is.nan(logProbs[i])) logProbs[i] <<- -Inf
          }
        }
      }
    }
    logProbs <<- logProbs - max(logProbs)
    probs <<- exp(logProbs)
    newValue <- rcat(1, probs) #rcat normalizes the probs internally
    if(newValue != currentValue) {
      model[[target]] <<- newValue
      calculate(model, calcNodes)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    } else nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(
    reset = function() { }
  )
)


