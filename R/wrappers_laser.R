## LASER doesn't actually provide simulation methods for most of its fit methods
# methods not written 

fitLaser <- function(object, method=c("pureBirth", "bd"), ...){
  branching.times(object)
}

update.fitLaser <- function(object, ...){ } 
simulate.fitLaser <- function(object, nsim=1, seed=NULL, ...){ }
getParameters.fitLaser<- function(fit){ }
loglik.fitLaser <- function(fit){}
get_data.fitLaser <- function(x, ...){}
get_phy.fitLaser <- function(x, ...){ }

# need a helper function to convert between a,r and b,d:
swp <- function(r, a){
    b = r/(1-a)
      d = r*a/(1-a)
        c(b,d)
}

