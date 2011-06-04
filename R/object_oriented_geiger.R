# First we need to teach geiger to treat its fits as proper objects 
# that can have methods associated with them


#  A geiger function they don't bother exposing.  Completely hides the second Brownian parameter!!!
phylogMean<-function(phyvcv, data){
	o<-rep(1, length(data))
	ci<-solve(phyvcv)
	m1<-solve(t(o) %*% ci %*% o)
	m2<-t(o) %*% ci %*% data
	m1 %*% m2
}
	


## should take "fit" output option as post wrapper to a fit continuous,
## along with tree, data, and model.  If that fit is null, it can
## make the call to fitContinuous
fitContinuous_object <- function(tree, data, model="BM", fit=NULL, ...){
# Args
#   fit: output from fitContinuous.  Will rerun if not provided
  if(is.null(fit))
  	fit <- fitContinuous(tree, data, model=model, ...)
	class(fit) <- "fitContinuous"
	fit$model <- model
	fit$tree <- tree
	fit$data <- data
	trans_tree <- transformTree(fit) 
	fit$root <- phylogMean( (fit[[1]]$beta)*vcv.phylo(trans_tree), data)
  fit$options <- list(opts=...)
	fit
}

update.fitContinuous <- function(fit, data){
  
	fitContinuous_object(tree=fit$tree, data=data, model=fit$model,
                       fit=NULL, bounds=fit$options$opts)
## This is not the best way to pass bounds
}

# rTraitCont might not be ideal here, should warn if alpha is large!
simulate.fitContinuous <- function(fit){
	tree <- transformTree(fit) 
	if(fit$model != "OU")
		data <- rTraitCont(tree, model="BM", sigma=sqrt(fit[[1]]$beta), root=fit$root)
	if(fit$model == "OU")
		data <- rTraitCont(fit$tree, model="OU", sigma=sqrt(fit[[1]]$beta), alpha=fit[[1]]$alpha, theta=fit$root, root=fit$root)

	## Former method used geiger's sim.char instead of ape's new rTraitCont.  Does the same thing
#	out<-sim.char(tree, as.matrix(fit[[1]]$beta), 1)
#	data <- out[,,1] + fit$root # probably not ok for the ou transform.  need to also account for theta 
#	names(data) <- names(fit$data)
	data
}

# then this'll work
loglik.fitContinuous <- function(fit){
	fit[[1]]$lnl
}

# then this'll work
getParameters.fitContinuous <- function(fit){
	unlist(fit[[1]])
}


# Note that this is not multivariate
transformTree <- function(fit){
	if(fit$model == "BM") out <- fit$tree
	else if(fit$model == "OU") out <- ouTree(fit$tree, fit[[1]]$alpha)
	else if(fit$model == "lambda") out <- lambdaTree(fit$tree, fit[[1]]$lambda)
	else if(fit$model == "kappa") out <- kappaTree(fit$tree, fit[[1]]$kappa)
	else if(fit$model == "delta") out <- deltaTree(fit$tree, fit[[1]]$delta)
	else print(paste("Transform model ", fit$model, " not recognized"))
	out
}

