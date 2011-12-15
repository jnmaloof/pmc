#objects_ouch.R

## Ouch is really already object oriented, but needs a few additional methods.  
## Also needs to tweak the "update" method (since ouch simulate returns a list of replicates)

#' S3 generic to provide a log-likelihood
#' @param x a model fit that can return a loglikelihood
#' @keywords internal
loglik <- function(x, ...) UseMethod("loglik")

#' Method to extract the log likelihood 
#' @return the log likelihood  
#' @method loglik hansentree
#' @S3method loglik hansentree
#' @keywords internal
loglik.hansentree <- function(object) object@loglik


#' Method to extract the log likelihood 
#' @return the log likelihood  
#' @method loglik hansentree
#' @S3method loglik hansentree
#' @keywords internal
loglik.browntree <- function(object) object@loglik


#' S3 generic to provide the parameters 
#' @param x a model fit of interest 
#' @keywords internal
getParameters <- function(x, ...) UseMethod("getParameters")

#' Method to extract the parameters 
#' @return a list of parameters
#' @method getParameters hansentree
#' @S3method getParameters hansentree
#' @keywords internal
getParameters.hansentree <- function(object){
	c(sigma=object@sigma, unlist(object@theta), sqrt.alpha=object@sqrt.alpha)
}

#' Method to extract the parameters 
#' @return a list of parameters
#' @method getParameters browntree
#' @S3method getParameters browntree
#' @keywords internal
getParameters.browntree <- function(object){
	c(sigma=object@sigma, unlist(object@theta))
}


## ouch cannot simulate before it has fit to data.  
## this makes a tree that can be simulated directly from set parameters
## not strictly needed pmc
make_browntree <- function(tree, sigma, theta){
	class(tree) <- "browntree"
	tree@nchar <- as.integer(1)
	tree@sigma <- sigma
	tree@theta <- list(theta)
	tree@data <- list(numeric(tree@nnodes))
	tree
}
