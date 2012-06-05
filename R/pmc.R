# pmc.R


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
#' ratio between the original MLE estimated models from the data. Also returns
#' each model object (A and B), and the estimates of each parameter for each
#' replicate in the element "par_dists", (see details), and call (the pmc 
#' function call used, for reference).  The returned object has S3 class pmc. 
#'
#' @details Possible models are all models from fitContinuous & ouch
#' 
#' The return value includes all parameters estimated under the structure
#' "par_dists", a data frame with columns "value", (the numerical value of
#' the parameter estimated) "parameter" (a factor indicating the name of the
#' parameter in the model, i.e. lambda), "comparison", (one of AA, AB, BA, or
#' BB, where the first letter indicates the class of model estimated, and 
#' the second indicates the model used to simulate the data on which the 
#' estimate is based. Hence AA is the bootstrap of model A, and BB the bootstrap
#' of model B.  Note that the cross-comparisons can also be informative).  
#' Finally "rep" indicates the replicate number for the simulation.  See examples
#' for plotting and calculating statistics from this data frame. 
#' @examples
#' require(geiger) # just for the sample data
#' data(geospiza)
#' attach(geospiza)
#' out <- pmc(geospiza.tree, geospiza.data[1], "BM", "lambda", nboot=5)
#' plot_pars(out$par_dists) # show the parameters  
#' ## Ex. mixing methods from packages -- data formats handled automatically 
#' ## Load Libraries ##
#' require(TreeSim) # to simulate a sample phyologeny
#' require(snowfall) # optional for parallelization
#' 
#' ## Simulate Data ##
#' simtree <- sim.bd.taxa(n=60, numbsim=1, lambda=1, mu=0, frac=1,
#'                        complete=FALSE, stochsampling=FALSE)[[1]][[1]] 
#' # simulates the EB model (by tranforming the tree)
#' dat <- rTraitCont(exponentialchangeTree(simtree, a=-0.5), sigma=5)
#' 
#' ## run PMC ## 
#' ## Specify the models and the options needed. 
#' regimes <- format_data(simtree, dat)$noregimes
#' eb_v_ou <- pmc(simtree, dat, modelA="EB", modelB="hansen", 
#'                optionsB=list(sqrt.alpha=1, sigma=1, regimes=regimes), 
#'                nboot=5) # small nboot just for example
#' plot(eb_v_ou) 
#' 
#' @import snowfall 
#' @import ouch
#' @export
# @enhance ouch
# @suggest TreeSim
pmc <- function(tree, data, 
                modelA = c("BM", "OU", "lambda", "kappa", "delta", "EB",
                           "white", "trend", "hansen", "brown"), 
                modelB = c("BM", "OU", "lambda", "kappa", "delta", "EB",
                           "white", "trend", "hansen", "brown"), 
                optionsA=list(), optionsB=list(), nboot=20, ...){

  A <- pmc_fit(tree, data, modelA, optionsA)
  B <- pmc_fit(tree, data, modelB, optionsB)
  lr_orig <- -2*(loglik(A) - loglik(B)) 

  reps <- sfLapply(1:nboot, function(i){

    null_sims <- simulate_and_update(A,B)
    AfitA <- null_sims[[1]]
    BfitA <- null_sims[[2]]

    ## Do the B sims
    test_sims <- simulate_and_update(B,A)
    BfitB <- test_sims[[1]] 
    AfitB <- test_sims[[2]]

    lrA <- -2*(loglik(AfitA) - loglik(BfitA)) 
    lrB <- -2*(loglik(AfitB) - loglik(BfitB))
    list(AA=as.list(c(lr=loglik(AfitA), getParameters(AfitA))), 
         BA=as.list(c(lr=loglik(BfitA), getParameters(BfitA))),
         AB=as.list(c(lr=loglik(AfitB), getParameters(AfitB))), 
         BB=as.list(c(lr=loglik(BfitB), getParameters(BfitB))))
  })
  ## some reformating, for convenience 
  dat <- melt(reps)
  names(dat) <- c("value", "parameter", "comparison", "rep")
  class(dat$value) <- "numeric"
  null = -2*(subset(dat, parameter=="lr" & comparison=="AA")$value -
             subset(dat, parameter=="lr" & comparison=="BA")$value)
  test = -2*(subset(dat, parameter=="lr" & comparison=="AB")$value - 
             subset(dat, parameter=="lr" & comparison=="BB")$value)
  output <- list(lr=lr_orig, null=null, test=test, par_dists=dat, 
                 A=A, B=B, call=match.call()) 
  class(output) <- "pmc"
  output
}


#' Internal helper function
#' Simulates under model A, updates both A and B based on that data
#' @param A model with simulate & update methods
#' @param B another model with simulate & update methods
#' @return a list of with fit of A on data simulated under A,
#' the fit on B on data simulated under A, and the simulated data itself
simulate_and_update <- function(A,B){
    error <- 1
    while(error){ ## Some sims don't converge. 

      simA <- simulate(A)
      simA <- matchformats(A, B, simA)
      AfitA <- try( update(A, data=simA[[1]]) )
      BfitA <- try( update(B, data=simA[[2]]) )

      if(!is(AfitA, "try-error") && !is(BfitA, "try-error") )
        error <- 0
      else
        warning("model fitting failed on data, simulating again")
    }
    list(AfitA, BfitA, simA)
}

#' Internal helper function to match data formats
#' @param A model with simulate & update methods
#' @param B another model with simulate & update methods
#' @param sim simulations produced by model A
#' @return a list of two data sets.  The first is in a format
#' appropriate to update model A, the second is in the format
#' approriate to update model B
matchformats <- function(A, B, sim){
  sim_for_updating_A <- sim
  sim_for_updating_B <- sim
  # reformat output from ouch's simulate method
  if (is(A, "ouchtree")){
    sim_for_updating_A <- sim$rep.1     
    sim_for_updating_B <- sim$rep.1
  }
  # The cross-comparison needs to check the data-formats match
  if(is(A, "ouchtree") & !is(B, "ouchtree")){
    ## this gets data into geiger format
    apeformatted <- format_data(get_phy(A), sim_for_updating_A)
    data <-apeformatted$data

    ## We just don't want data frame outputs to pass to geiger
    if(is(data, "data.frame")){
      tmp <- data[[1]] # we want numeric data
      names(tmp) <- rownames(data)
      data <- tmp
#    tmp <- na.exclude(data)
#    tmp2 <- as.numeric(tmp)
#    names(tmp2) <- names(tmp)
#    data <- tmp2
    }
    sim_for_updating_B <- data
  }
  # The cross-comparison needs to check the data-formats match
  if(!is(A, "ouchtree") & is(B, "ouchtree")){
    sim_for_updating_B <- format_data(get_phy(A),
                                      sim_for_updating_A)$data
  }
list(sim_for_updating_A, sim_for_updating_B)
}

#' plot the distributions
#' @param x a pmc object
#' @param ... Additional arguments: 
#'  A= a name for the first model in the pmc pairwise comparison
#'  B= a name for the second model in the pairwise comparison
#' @import ggplot2
#' @import reshape2
#' @method plot pmc
#' @S3method plot pmc
plot.pmc <- function(x, ...){
  df <- data.frame(x$null, x$test)
  colnames(df) <- c("null", "test")
  dat <- melt(df)
  # FIXME strange that ggplot wants 'value' and 'variable' to be unquoted!
  # Makes these appear to be undefined values & creates a warning in check  
  ggplot(dat) + geom_density(aes(value, fill=variable), alpha=.7) +
       geom_vline(x=x$lr, lwd=1, lty=2)
}

#' plot the parameter distributions
#' @param object a pmc object fit
#' @return a ggplot2 plot object
#' @export
plot_pars <- function(object){
## add lines for the estimated parameter values of A & B
  ggplot(object) + geom_boxplot(aes(comparison, value)) + facet_wrap(~parameter, scales="free_y") 
}



#' Fit any model used in PMC 
#'
#' The fitting function used by pmc to generalize fitting to any model
#' @param tree a phylogenetic tree. can be ouch or ape format
#' @param data trait data in ape or ouch format
#' @param model the name of the model to fit, see details for a list of 
#' currently supported types
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
#' ## Or a single trait at a time: 
#' lambdaFit<-pmc_fit(geospiza.tree, geospiza.data[1], model="lambda") 
#' ## an ouch example
#' require(ouch) # just for the data, 
#' data(bimac)
#' tree <- with(bimac,ouchtree(node,ancestor,time/max(time),species))
#' ou.3 <- pmc_fit(data=log(bimac['size']),tree, model="hansen", 
#'                 list(regimes=bimac['OU.3'],sqrt.alpha=1,sigma=1))
#' @import geiger 
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
      if(is(data, "data.frame"))
        data <- data[[1]] # we want numeric data
      names(data) <- tree@nodelabels
      tmp <- na.exclude(data)
      tmp2 <- as.numeric(tmp)
      names(tmp2) <- names(tmp)
      data <- tmp2
      tree <- convert(tree)
    }
    args <- c(list(tree=tree, data=data, model=model), options) 
    object <- do.call(fitContinuous_object, args) 
  } else if(type == "hansen"){
    ## Fit an ouch object (hansen) 
    if(is(tree, "phylo")){
      tmp <- format_data(tree, data)
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
