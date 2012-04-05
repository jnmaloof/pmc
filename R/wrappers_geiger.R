# file: object_oriented_geiger.R 
# author: Carl Boettiger <cboettig@gmail.com>
# license: BSD 
# Description: A set of functions to make Geiger's fits into objects that can have 
#   methods associated with them using R's S3 Classes



#  A geiger function they don't bother exposing. Copied here to access the 
# fitted value for the phylogenentic mean under BM.  
phylogMean<-function(phyvcv, data){
	o<-rep(1, length(data))
	ci<-solve(phyvcv)
	m1<-solve(t(o) %*% ci %*% o)
	m2<-t(o) %*% ci %*% data
	m1 %*% m2
}
	

#' wrapper for fitContinuous
#' 
#' Runs the fitContinuous function and returns and object with necessary
#' information to run methods such as 'update' and 'simulate'.  Call is the
#' same as the original fitContinuous function.
#' @param tree a phylo object phylogenetic tree
#' @param data ape data phylo format
#' @param model the type of model to fit
#' @param bounds optional list of bounds for the parameter fits
#' @param meserr measurment error, if any
#' @return an object of class fitContinuous 
## should take "fit" output option as post wrapper to a fit continuous,
## along with tree, data, and model.  If that fit is null, it can
## make the call to fitContinuous
fitContinuous_object <- function(tree, data, model="BM", bounds=NULL, 
                                 meserr=NULL){
  tmp <- treedata(tree,data)
  data <- tmp$data
  tree <- tmp$phy
  fit <- fitContinuous(phy=tree, data=data, model=model, 
                       bounds=bounds, meserr=meserr)
  for(i in 1:length(fit)){
    ## This makes lots of copies of the data if you have multiple 
    ## traits in a data frame.  For effiency, please don't run large
    ## numbers of traits in matrix form, but sapply over them instead.
    fit[[i]]$tree <- tree 
    fit[[i]]$data <- data
    fit[[i]]$model <- model
    fit[[i]]$bounds <- bounds
    fit[[i]]$meserr <- meserr
  	trans_tree <- transformTree(fit[[i]])

    if(model != "white"){
      if(is(data, "data.frame")){
        fit[[i]]$root <- phylogMean((fit[[i]]$beta) * 
                                     vcv.phylo(trans_tree), data[[i]])
      } else {
        fit[[i]]$root <- phylogMean((fit[[i]]$beta) * 
                                     vcv.phylo(trans_tree), data)
      }
    } else { ## "white" model
      fit[[i]]$root <- fit[[i]]$mean
      fit[[i]]$beta <- fit[[i]]$nv
    }
  }
	class(fit) <- "fitContinuous"
  fit
}

#' update a fitContinuous model
#' @param object the output of fitContinuous
#' @param ... data on which update will be based
#' @return the MLE estimated model for the given data
#' @method update fitContinuous
#' @S3method update fitContinuous
update.fitContinuous <- function(object, ...){
	fitContinuous_object(tree=object[[1]]$tree, ...,
                       model=object[[1]]$model,
                       bounds=object[[1]]$bounds, 
                       meserr=object[[1]]$meserr)
}


#' simulate method
#' @param object a fitContinuous object
#' @param nsim number of sims (currently always 1) 
#' @param seed an optional seed for the simulations (not implemented)
#' @param ... additional arguments, not implemented for fitContinuous simulations
#' @return simulated dataset
# These next lines would generate the the S3 method line for the NAMESPACE,
# but this causes a conflict with the S4 method so it has been removed.  
# instead we export this function
#' @method simulate fitContinuous
#' @S3method simulate fitContinuous
#' @details intended as an internal function, though an available S3 method
#' @import ape
# properly speaking we don't want to export this, but it solves the s4 conflict
#' @export
simulate.fitContinuous <- function(object, nsim=1, seed=NULL, ...){
# rTraitCont might not be ideal here, should warn if alpha is large!
#  Currently only designed to simulate from the first object
  data <- object[[1]]$data 
  for(i in 1:length(object)){
    if(object[[i]]$model != "OU"){
      tree <- transformTree(object[[i]]) 
    # transform the tree and apply the BM simulation method
      sim <- rTraitCont(tree, model="BM", sigma=sqrt(object[[i]]$beta),
                              root.value=object[[i]]$root)
      # Make sure the simulation species names are in the order given by the data
      # (rTraitCont gives them in the order given in the tree, which need not match!)
      # Thanks to Jason Weir for catching this.  
      data[, i] <- sim[match(rownames(data), names(sim))]
    } else if(object[[i]]$model == "white"){ 
      data[,i] <- rnorm(length(data[,i]), mean=mean(object[[i]]$data), sd = sqrt(object[[i]]$beta))
    } else { 
      # use untransfromed tree & with the OU simulation method
      data[,i] <- rTraitCont(object[[i]]$tree, model="OU", 
                             sigma=sqrt(object[[i]]$beta),
                             alpha=object[[i]]$alpha, theta=object[[i]]$root, 
                             root.value=object[[i]]$root)
    }
  }
  # Keep labels for traits and species on data 
  names(data) <- names(object[[1]]$data)
  rownames(data) <- rownames(object[[1]]$data)
	data
}

#' Method to extract the log likelihood
#' @param fit a fitContinuous object
#' @return the log likelihood, currently only for a single/first trait
#' @method loglik fitContinuous 
#' @S3method loglik fitContinuous  
#' @keywords internal
loglik.fitContinuous <- function(fit){
	fit[[1]]$lnl
}


#' Method to extract the parameters 
#' @param fit a fitContinuous object
#' @return a list of parameters, currently only for a single/first trait
#' @method getParameters fitContinuous  
#' @S3method getParameters fitContinuous  
#' @keywords internal
getParameters.fitContinuous <- function(fit){
	who <- match(c("lnl", "beta", "alpha", "root", "lambda", "k", "kappa", "delta", "a"),  names(fit[[1]]))
  pars <- fit[[1]][na.omit(who)]
  out <- as.numeric(pars)
  names(out) <- names(pars) 
  out
}


# Note that this is not multivariate, and is actually passed a fit[[i]] object
# This preserves the "vector" approach of geiger (which isn't actually 
# a vector/multivariate analysis so is stupid notation!) but is memory-inefficient 
# @keywords internal
transformTree <- function(fit){
	if(fit$model == "BM" | fit$model == "white") out <- fit$tree
	else if(fit$model == "OU") out <- ouTree(fit$tree, fit$alpha)
	else if(fit$model == "lambda") out <- lambdaTree(fit$tree, fit$lambda)
	else if(fit$model == "kappa") out <- kappaTree(fit$tree, fit$kappa)
	else if(fit$model == "delta") out <- deltaTree(fit$tree, fit$delta)
	else if(fit$model == "EB") out <- exponentialchangeTree(fit$tree, a=fit$a)
	else print(paste("Transform model ", fit$model, " not recognized"))
	out
}



#' Method to grab the phylogeny 
#' @return the phylogeny of the object 
#' @method get_phy fitContinuous 
#' @S3method get_phy fitContinuous
#' @keywords internal
get_phy.fitContinuous <- function(x, ...) x[[1]]$tree

#' Method to grab the data 
#' @return the trait data of the object 
#' @method get_data fitContinuous 
#' @S3method get_data fitContinuous
#' @keywords internal
get_data.fitContinuous <- function(x, ...) x[[1]]$data


## Add fitDiscrete functionality !! ###

