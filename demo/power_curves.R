#power_curves.R
require(pmc)
require(TreeSim)
require(ouch)

nboot <- 200
cpu <- 4

alpha  <- c(.01, .1, 1, 1.5, 2, 5, 10, 20)
n      <- c(5, 10, 20, 50, 100)
lambda <- c(.25, .5, .75, 1)

# size simulations
size <- lapply(1:length(n), function(i){
	simtree <- sim.bd.taxa(n=n[i], numbsim=1, lambda=1, mu=0, frac=1, complete=FALSE, stochsampling=FALSE)[[1]][[1]] 
	treepower(ape2ouch(simtree), nboot=nboot, cpu=cpu, alpha=alpha)
})
## Shape simulations
N <- 50 ## number of taxa fixed to 50
shape <- lapply(1:length(lambda), function(i){
	simtree <- sim.bd.taxa(n=N, numbsim=1, lambda=1, mu=0, frac=1, complete=FALSE, stochsampling=FALSE)[[1]][[1]]
	simtree <- lambdaTree(simtree, lambda[i])
	treepower(ape2ouch(simtree), nboot=nboot, cpu=cpu, alpha=alpha)
})
## Do the Anoles tree for comparison
data(bimac) # ouch package Anolis sizes (from N. Lesser Antilles)
tree <- with(bimac,ouchtree(node,ancestor,time/max(time),species))
anoles <- treepower(tree, nboot=nboot, cpu=cpu, alpha=alpha )


############## Figure 6 #################
# n is a vec of number of taxa used in each sim: 
plot_size <- function(){
    k <- length(n)
  plot(1,1, type='n', xlim=c(min(alpha), max(alpha)), ylim = c(0,1),
       main="Power by tree size", log="x", xlab="alpha", ylab="power")
  for(i in 1:k ){ ## skip the last 2, which haven't converged
    points(alpha, size[[i]]$power, pch=16, col=i)
    lines(alpha, size[[i]]$power, col=i)
  }
  points(alpha, anoles$power, pch=16, col="purple")
  lines(alpha, anoles$power, col="purple", lwd=4)
  legend("topleft", c(paste(n[1:3], "taxa"), "23 (anoles)"), col=c(1:k,"purple"), pch=16  ) 

}
plot_size()

plot_shape <- function(){
  plot(1,1, type='n', xlim=c(min(alpha), max(alpha)), ylim = c(0,1),
  main="Power by tree topology", log="x", xlab="alpha", ylab="power")
    k <- length(lambda)
  for(i in 1:length(n)){
    points(alpha, shape[[i]]$power, pch=16, col=i)
    lines(alpha, shape[[i]]$power,  col=i)
  }
  points(alpha, anoles$power, pch=16, col="purple")
  lines(alpha, anoles$power, col="purple", lwd=4)
  legend("topleft", c(paste(lambda, "lambda"), "anoles"), col=c(1:k, "purple"), pch=16  ) 
}
plot_shape()

