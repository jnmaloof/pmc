#pmc.R


#' Performs a phylogenetic monte carlo between modelA and modelB
#'
#' Simulates data under each model and returns the distribution of 
#' likelihood ratio, L(B)/L(A), under for both simulated datasets.
#' @param tree A phylogenetic tree.  Can be phylo (ape) or ouch tree
#' @param data The data matrix
#' @param modelA a model from the list, or a custom model, see details
#' @param modelB any other model from the list, or custom model, see details
#' @param optionsA a named list of arguments for fitting model A
#' @param optionsB a named list of arguments for fitting model B
#' @param nboot number of bootstrap replicates to use
#' @param ... additional arguments
#' @return 
#' list with the nboot likelihood ratios obtained from fitting both models
#' to data simulated by model A, and the nboot likelihood ratios obtained
#' by fitting both models to simulations from model B, and the likelihood 
#' ratio between the original MLE estimated models from the data.  
#' return has object class pmc. 
#' @details Possible models are all models from fitContinuous, ouch's hansen model
#' @examples
#' require(geiger) # just for the sample data
#' data(geopsiza)
#' attach(geospiza)
#' out <- pmc(geospiza.tree, geospiza.data, "BM", "lambda")
#' ## FIXME add examples that require conversion of data formats
#' @import snowfall reshape
#' @export
pmc <- function(tree, data, 
                modelA = c("BM", "OU", "lambda", "kappa", "delta", "EB", "white", "trend", "hansen"), 
                modelB = c("BM", "OU", "lambda", "kappa", "delta", "EB", "white", "trend", "hansen"), 
                optionsA=list(), optionsB=list(), nboot=20, ...){

  A <- pmc_fit(tree, data, modelA, optionsA)
  B <- pmc_fit(tree, data, modelB, optionsB)
  lr_orig <- -2*(loglik(A) - loglik(B)) 

  reps <- sfLapply(1:nboot, function(i){
    ## Do the A sims
    simA <- simulate(A)
    if (is(A, "ouchtree"))
      simA <- simA$rep.1
    AfitA <- update(A, simA, ...)
    BfitA <- update(B, simA, ...)
    lrA <- -2*(loglik(AfitA) - loglik(BfitA)) 

    ## Do the B sims
    simB <- simulate(B)
    if (is(B, "ouchtree"))
      simB <- simB$rep.1
    AfitB <- update(A, simB, ...)
    BfitB <- update(B, simB, ...)
    lrB <- -2*(loglik(A) - loglik(B)) 
    list(lrA=lrA, parsAA=getParameters(AfitA), parsBA=getParameters(BfitA),
         lrB=lrB, parsAB=getParameters(AfitB), parsBB=getParameters(BfitB))
  })
  ## some reformating 
  dat <- melt(reps)
  names(dat) <- c("value", "variable", "rep")
  null = subset(dat, variable=="lrA")$value
  test = subset(dat, variable=="lrB")$value
  par_dists = subset(dat, !(variable %in% c("lrA", "lrB")))

  output <- list(lr=lr_orig, null=null, test=test, par_dists=par_dists, 
                 A=A, B=B, call=match.call()) 
  class(output) <- "pmc"
  output
}

#' Fit any model used in PMC 
#' @param tree a phylogenetic tree. can be ouch or ape format
#' @param data trait data in ape or ouch format
#' @param options whatever additional options would be provided 
#' to the model fit, see details
#' @return a pmc_model object, anything that has methods "simulate", 
#' "update", getParameters, and getLikelihood
#' @details options should include all parameters required by the fit method
#' Currently methods avialable are fitContinuous (see geiger package) and 
#' hansen (see the ouch package).  
#' @examples 
#' ## a geiger example
#' require(geiger) # just to load the data
#' data(geospiza)
#' attach(geospiza)
#' lambdaFit<-pmc_fit(geospiza.tree, geospiza.data, model="lambda") 
#' Or a single trait at a time 
#' lambdaFit<-pmc_fit(geospiza.tree, geospiza.data[1], model="lambda") 
#' ## an ouch example
#' require(ouch) # just for the data, 
#' data(bimac)
#' tree <- with(bimac,ouchtree(node,ancestor,time/max(time),species))
#' ou.2 <- pmc_fit(data=log(bimac['size']),tree, model="hansen", 
#'                 list(regimes=bimac['OU.2'],sqrt.alpha=1,sigma=1))
#' @import geiger ouch
#' @export
pmc_fit <- function(tree, data, model, options=list()){
  # Figure out if we need ape/geiger based formats or ouch formats
  fitContinuous_types <- c("BM", "OU", "lambda", "kappa", 
                           "delta", "EB", "white", "trend")

  if(model %in% fitContinuous_types){
    type <- "fitContinuous"  
  } else if(model %in% c("brown", "hansen")){
    type <- "hansen"
  } else {
    stop(paste("Model", model, "not recognized"))
  }
  ## Run a fitContinuous fit ##
  if(type == "fitContinuous"){
    # first, check data formats
    if(is(tree, "ouchtree")){
      # assumes data is same order as nodelabels
      names(data) <- tree@nodelabels
      data <- data[data!=NA]
      tree <- convert(tree)
    }
    args <- c(list(tree=tree, data=data, model=model), options) 
    object <- do.call(fitContinuous_object, args) 
  } else if(type == "hansen"){
    ## Fit an ouch object (hansen) 
    if(is(tree, "phylo")){
      tmp <- format_data(tree, traits)
      tree <- tmp$tree
      data <- tmp$data
    }
    args <- c(list(data=data, tree=tree), options) 
   if(model == "hansen")
       object <- do.call(hansen, args)
   else if(model == "brown")
       object <- do.call(brown, args)
  } else {
    stop("Error: format not recognized.")
  }
  object
}
