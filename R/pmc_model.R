# Example get_sim_input wrapper 
pureBirth2sim.bd.taxa <- function(m){ 
  list(n = m$fit_input$phy$Nnode+1, numbsim=1,
         lambda=m$output$b, mu=0, frac = 1, 
                complete = FALSE, stochsampling = FALSE)
}
# example get_sim_output wrapper
TreeSim2laser <- function(x){
  branching.times(x[[1]])
}



pmc_model <- function(fit_method, simulate_method, fit_input, 
                      data_name = NULL, get_sim_input, get_sim_output,
                      get_loglik)
# create an object that we can use for pmc methods
# Args
#   fit_method: a function (or character string with the function name)
#               used to generate the fit, i.e. pureBirth
#   simulate_method: a function (or function name) to create simulated data
#   fit_input: a named list of the arguments of fit_method function
#   data_name: a character or index of where the data (such as is produced
#              by the simulation method) is given in the fit_input list.  If
#              not specified, we will try and use the first slot. 
#   get_sim_input: a function which creates a named list of inputs to the 
#                   simulation function from a pmc_model object; see Details
#   get_sim_output: a function which reformat the output of the simulation
#                   to match that taken in fit_input
#   get_loglik: a function to get the loglikelihood back from the output of 
#               the fit_method function
#
# Examples:
#      pb <- function(fit_input){
#            pmc_model(pureBirth, sim.bd.taxa, 
#                     fit_input = fit_input, 
#                     data_name = "x",
#                     get sim input = function(m) list(n = length(m$fit_input$x), 
#                         numbsim=1, lambda=m$fit_results$r1, mu=0, frac = 1, 
#                         complete = FALSE, stochsampling = FALSE),
#                     get_sim_output = function(x) branching.times(x[[1]][[1]]),
#                     get_loglik = function(fit_results) fit_results$LH ) 
#      }
#       require(pmc); data(geospiza) fit_input = list(x=branching.times(geospiza$geospiza.tree))
#       m <- pb(fit_input)
#       x <- simulate(m)
#       m2<-update(m,x)
#       print(loglik(m2))
{
  fit_results <- do.call(fit_method, fit_input)

  if(is.null(data_name))
    data_name <- 1

  out <- list(fit_results=fit_results, fit_input=fit_input, fit_method=fit_method, 
              simulate_method=simulate_method, get_sim_input=get_sim_input,
              get_sim_output=get_sim_output, get_loglik=get_loglik)
  class(out) <- "pmc_model"
  out
}



simulate.pmc_model <- function(m){
# simulate a pmc_model object
# Args: m: a pmc_model object, see example
  x <- do.call(m$simulate_method, m$get_sim_input(m))
  m$get_sim_output(x)
}

update.pmc_model <- function(m, x){
# update a pmc model fit
# Args
#   m: an object of class pmc_model (contains given data and the methods)
#   x: the data the model will use to update, can come from simulate.pmc_model
  m$fit_input[m$data_name] <- x
  m$fit_returns <- do.call(m$fit_method, m$fit_input) 
  m
}

getParameters.pmc_model <- function(m) m$fit_results

loglik.pmc_model <- function(m) m$get_loglik(m$fit_results)

pmc <- function(null, test, nboot=100, cpu=1, getParNames=FALSE){
  ## parallelization 
	if(cpu>1 & !sfIsRunning()){ 	
		sfInit(parallel=TRUE, cpu=cpu) 
		sfLibrary(pmc)
		sfExportAll()
	} else if(cpu<2 & !sfIsRunning()){
    sfInit()
  }


  ## Perform the bootstrapping
  null_sim <- sfSapply(1:nboot, function(i){
    data <- simulate(null)
    null <- update(null, data)
    test <- update(test, data)
    lr <- -2*(loglik(null) - loglik(test)) 
	 list(lr, getParameters(null), getParameters(test))
	})
	test_sim <- sfSapply(1:nboot, function(i){
		data <- simulate(test)
		null <- update(null, data)
		test <- update(test, data)
		lr <- -2*(loglik(null) - loglik(test))
		list(lr, getParameters(null), getParameters(test))
	})


  ## Formatting output
	null_dist <- unlist(null_sim[1,]) 
	test_dist <- unlist(test_sim[1,]) 
	null_bootstrap_pars <- matrix( unlist(null_sim[2,]), ncol=nboot)
	test_bootstrap_pars <- matrix( unlist(test_sim[3,]), ncol=nboot)
	null_sim_test_pars  <- matrix( unlist(null_sim[3,]), ncol=nboot)
  test_sim_null_pars  <- matrix( unlist(test_sim[2,]), ncol=nboot)
  if(GetParNames){ 
    rownames(null_bootstrap_pars) <- names(getParameters(null))
    rownames(test_bootstrap_pars) <- names(getParameters(test))
    rownames(null_sim_test_pars) <- names(getParameters(test))
    rownames(test_sim_null_pars) <- names(getParameters(null))
  }
  output <- list(null_dist=null_dist, test_dist=test_dist, nboot=nboot, 
                 null=null, test=test, null_par_dist = null_bootstrap_pars, 
                 test_par_dist = test_bootstrap_pars, null_sim_test_pars = 
                 null_sim_test_pars, test_sim_null_pars=test_sim_null_pars)
	class(output) <- "pow"
	output

}
