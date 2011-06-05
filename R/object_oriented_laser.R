#object_oriented_laser.R
pureBirth_fit <- function(fit_input){
      pmc_model(pureBirth, sim.bd.taxa, 
               fit_input = fit_input, 
               data_name = "x",
               get_sim_input = function(m) 
                                  list(n = length(m$fit_input$x), numbsim=1,
                                  lambda=m$fit_results$r1, mu=0, frac = 1, 
                                  complete = FALSE, stochsampling = FALSE),
               get_sim_output = function(x) branching.times(x[[1]][[1]]),
               get_loglik = function(fit_results) fit_results$LH ) 
 
}


# need a helper function to convert between a,r and b,d:
swp <- function(r, a){
    b = r/(1-a)
      d = r*a/(1-a)
        c(b,d)
}

bd_fit <- function(fit_input){
      pmc_model(bd, sim.bd.taxa, 
               fit_input = fit_input, 
               data_name = "x",
               get_sim_input = function(m) 
                                  list(n = length(m$fit_input$x), numbsim=1,
                                  lambda = swp(m$fit_results$r, m$fit_results$a)[1], 
                                  mu = swp(m$fit_results$r, m$fit_results$a)[2], 
                                  frac = 1, complete = FALSE, stochsampling = FALSE),
               get_sim_output = function(x) branching.times(x[[1]][[1]]),
               get_loglik = function(fit_results) fit_results$LH ) 
}


## LASER doesn't actually provide simulation methods for most of its fit methods
