# need a helper function to convert between a,r and b,d:
swp <- function(r, a){
    b = r/(1-a)
      d = r*a/(1-a)
        c(b,d)
}


get_lambdas <- function(res){
## Pull the lambda values out from the treepar output
 res[[2]][[i]]
}
get_mus <- function(res){
}


treepar_fit <- function(fit_input){
    pmc_model(bd.shifts.optim, sim.rateshifts.taxa,
              fit_input=fit_input, data="x",
              get_sim_input = function(m)
                list(n = length(m$fit_input$x), numbsim=1,
                     lambda = get_lambda(m$fit_results),
                     mu <- get_mus(m$fit_results),
                     frac=1, times, complete = FALSE),
              get_sim_output = function(x) sort(getx(x[[1]][[1]]),decr=T),
              get_loglik = function(fit_results) fit_results
              )
                
}


# Treepar output is something of a nightmare:
#res[[2]][[i]]: Maximum likelihood parameter estimates for i-1 shifts (i
#          in 1:m): First entry is the (-log likelihood) value. The next
#          i entries are the turnover (extinction/speciation) estimates,
#          for the successive intervals going back in time. The next i
#          entries are the diversification rate estimates
#          (speciation-extinction). The next i-1 entries are the
#          sampling estimates (if ME=TRUE). The last i-1 entries are the
#          shift times. (Note: if ME=TRUE and all==FALSE, the second
#          entry is the turnover, the third the diversification rate,
#          followed by the sampling estimates).
#

