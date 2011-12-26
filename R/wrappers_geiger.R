# file: object_oriented_geiger.R 
# author: Carl Boettiger <cboettig@gmail.com>
# license: BSD 
# Description: A set of functions to make Geiger's fits into objects that can have 
#   methods associated with them using R's S3 Classes



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
fitContinuous_object <- function(tree, data, model="BM", bounds=NULL, meserr=NULL, fit=NULL){
  tmp <- treedata(tree,data)
  data <- tmp$data
  tree <- tmp$phy
  fit <- fitContinuous(phy=tree, data=data, model=model, bounds=bounds, meserr=meserr)
  for(i in 1:length(fit)){
    ## This makes lots of copies of the data if you have multiple traits in a data frame
    ## For effiency, please don't run large trait matrix simultaneously!!!!!
    fit[[i]]$tree <- tree 
    fit[[i]]$data <- data
    fit[[i]]$model <- model
    fit[[i]]$bounds <- bounds
    fit[[i]]$meserr <- meserr
  	trans_tree <- transformTree(fit[[i]])
    if(is(data, "data.frame")){
	    fit[[i]]$root <- phylogMean( (fit[[i]]$beta)*vcv.phylo(trans_tree), data[[i]] )
    } else {
 	    fit[[i]]$root <- phylogMean( (fit[[i]]$beta)*vcv.phylo(trans_tree), data)
    }
#    class(fit[[i]]) <- "fitContinuous"
  }
	class(fit) <- "fitContinuous"
  fit
}

#' update a fitContinuous model 
#' @return the MLE estimated model for the given data
#' @method update fitContinuous
#' @S3method update fitContinuous
update.fitContinuous <- function(fit, data){
	fitContinuous_object(tree=fit[[1]]$tree, data=data, model=fit[[1]]$model,
                       bounds=fit[[1]]$bounds, meserr=fit[[1]]$meserr)
}


#' simulate method
#' @return simulated dataset
#' @method simulate fitContinuous
#' @S3method simulate fitContinuous
#' @details rTraitCont might not be ideal here, should warn if alpha is large!
#'  Currently only designed to simulate from the first fit
simulate.fitContinuous <- function(fit){
  data <- data.frame(NULL)
  for(i in 1:(dim(fit)[2])){
    if(fit[[i]]$model != "OU"){
      tree <- transformTree(fit[[1]]) 
      data[[i]] <- rTraitCont(tree, model="BM", sigma=sqrt(fit[[i]]$beta),
                              root=fit[[i]]$root)
    } else if(fit[[i]]$model == "OU"){
      data[[i]] <- rTraitCont(fit$tree, model="OU", sigma=sqrt(fit[[i]]$beta),
                              alpha=fit[[i]]$alpha, theta=fit[[i]]$root, 
                              root=fit[[i]]$root)
    }
  }
	data
}

#' Method to extract the log likelihood 
#' @return the log likelihood  
#' @method loglik fitContinuous 
#' @S3method loglik fitContinuous  
loglik.fitContinuous <- function(fit){
	fit[[1]]$lnl
}


#' Method to extract the parameters 
#' @return a list of parameters
#' @method getParameters fitContinuous  
#' @S3method getParameters fitContinuous  
getParameters.fitContinuous <- function(fit){
	unlist(fit[[1]])
}


# Note that this is not multivariate, and is actually passed a fit[[i]] object
# This preserves the "vector" approach of geiger (which isn't actually 
# a vector/multivariate analysis so is stupid notation!) but is memory-inefficient 
transformTree <- function(fit){
	if(fit$model == "BM") out <- fit$tree
	else if(fit$model == "OU") out <- ouTree(fit$tree, fit$alpha)
	else if(fit$model == "lambda") out <- lambdaTree(fit$tree, fit$lambda)
	else if(fit$model == "kappa") out <- kappaTree(fit$tree, fit$kappa)
	else if(fit$model == "delta") out <- deltaTree(fit$tree, fit$delta)
	else if(fit$model == "EB") out <- exponentialchangeTree(fit$tree, a=fit$a)
	else print(paste("Transform model ", fit$model, " not recognized"))
	out
}

