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
                      get_loglik, get_parameters=NULL)
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
#   get_parameters: an optional function that returns a numeric listing the 
#                   model parameters given the results of fit.  If left as NULL,
#                   will just include all fit results
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
#                     get_parameters = function(fit_results) fit_results$r1
#      }
#       require(pmc); data(geospiza) fit_input = list(x=branching.times(geospiza$geospiza.tree))
#       m <- pb(fit_input)
#       x <- simulate(m)
#       m2<-update(m,x)
#       print(loglik(m2))
{
  fit_results <- do.call(fit_method, fit_input)
  if(is.null(get_parameters))
    get_parameters <- function(fit_results) fit_results
  if(is.null(data_name))
    data_name <- 1

  out <- list(fit_results=fit_results, fit_input=fit_input, 
              fit_method=fit_method, simulate_method=simulate_method,
              get_sim_input=get_sim_input, get_sim_output=get_sim_output,
              get_loglik=get_loglik, get_parameters=get_parameters)
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

getParameters.pmc_model <- function(m) m$get_parameters(m$fit_results)

loglik.pmc_model <- function(m) m$get_loglik(m$fit_results)

