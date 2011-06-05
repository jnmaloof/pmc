# pmc_model makes it easier to write wrappers for phylogenetic methods that don't have update and simulate methods.  


# use the pmc_model to define a function that fits under "pureBirth" and simulates under "sim.bd.taxa"
pureBirth_fit <- function(fit_input){
      pmc_model(pureBirth, sim.bd.taxa, 
               fit_input = fit_input, 
               data_name = "x",
               function(m) list(n = length(m$fit_input$x), numbsim=1,
                                  lambda=m$fit_results$r1, mu=0, frac = 1, 
                                  complete = FALSE, stochsampling = FALSE),
               function(x) branching.times(x[[1]][[1]]),
               function(fit_results) fit_results$LH ) 
 
}



# ouch fits and simulations
ouch <- function(fit_input){
      pmc_model(hansen, simulate, 
               fit_input = fit_input, #list(data=, tree=, regimes=, sqrt.al=, sigma=)
               data_name = "data", # is first, can ignore
               function(m) list(m$fit_results), # already has a method to simulate these
               function(x) x$rep.1,
               function(fit_results) fit_results@loglik ) 
 
}

require(ouch)
data(bimac)
tree <- with(bimac,ouchtree(node,ancestor,time/max(time),species))
m <- ouch( list(data=log(bimac['size']),tree=tree,regimes=bimac['OU.1'],sqrt.alpha=1,sigma=1) )
x <- simulate(m)
m2 <- update(m,x)

#fit_input <- formals(hansen)
#sapply(fit_input, function(x) is.name(x)) # which ones are empty?

require(pmc)
data(geospiza)
fit_input = list(x=branching.times(geospiza$geospiza.tree))
m <- pureBirth_fit(fit_input)
x <- simulate(m)
m2<-update(m,x)
print(loglik(m2))
