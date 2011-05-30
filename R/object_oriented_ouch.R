#objects_ouch.R

## Ouch is really already object oriented, but needs a few additional methods.  
## Also needs to tweak the "update" method (since ouch simulate returns a list of replicates)
loglik.hansentree <- function(object) object@loglik
loglik.browntree <- function(object) object@loglik
getParameters.hansentree <- function(object){
	c(sigma=object@sigma, unlist(object@theta), sqrt.alpha=object@sqrt.alpha)
}
getParameters.browntree <- function(object){
	c(sigma=object@sigma, unlist(object@theta))
}






## ouch cannot simulate before it has fit to data.  this makes a tree that can be simulated directly from set parameters
## not strictly needed pmc
make_browntree <- function(tree, sigma, theta){
	class(tree) <- "browntree"
	tree@nchar <- as.integer(1)
	tree@sigma <- sigma
	tree@theta <- list(theta)
	tree@data <- list(numeric(tree@nnodes))
	tree
	}
