#power_curves.R

###############
require(socialR)
script <- "power_curves.R"
gitaddr <- gitcommit(script)
tags="phylogenetics"
tweet_errors(script, tags=tags)
###############



require(pmc)
require(TreePar)
require(ouch)

nboot <- 500
cpu <- 16

alpha  <- c(.01, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 2, 5, 10, 20, 50)
n      <- c(10,  20, 40, 60, 80, 100)
lambda <- c(.01, .1, .4, .6, .8, 1)


size <- lapply(1:length(n), function(i){
	simtree <- sim.bd.taxa(n=n[i], numbsim=1, lambda=1, mu=0, frac=1, complete=FALSE, stochsampling=FALSE)[[1]][[1]] 
	treepower(ape2ouch(simtree), nboot=nboot, cpu=cpu, alpha=alpha)
})

## number of taxa
N <- 50
shape <- lapply(1:length(lambda), function(i){
	simtree <- sim.bd.taxa(n=N, numbsim=1, lambda=1, mu=0, frac=1, complete=FALSE, stochsampling=FALSE)[[1]][[1]]
	simtree <- lambdaTree(simtree, lambda[i])
	treepower(ape2ouch(simtree), nboot=nboot, cpu=cpu, alpha=alpha)
})


save(file="power_curves.Rdat", list=ls() )

## Do the Anoles tree for comparison
data(bimac) # ouch package Anolis sizes (from N. Lesser Antilles)
tree <- with(bimac,ouchtree(node,ancestor,time/max(time),species))
anoles <- treepower(tree, nboot=nboot, cpu=cpu, alpha=alpha )

save(file="power_curves.Rdat", list=ls() )


############## Figure 6 #################
# n is a vec of number of taxa used in each sim: 
plot_size <- function(){
    k <- length(n)-2 ## skip the last 2, which haven't converged
  plot(1,1, type='n', xlim=c(min(alpha), max(alpha)), ylim = c(0,1), main="Power by tree size", log="x", xlab="alpha", ylab="power")
  for(i in 1:k ){
    points(alpha, size[[i]]$power, pch=16, col=i)
    lines(alpha, size[[i]]$power, col=i)
  }
  points(alpha, anoles$power, pch=16, col="purple")
  lines(alpha, anoles$power, col="purple", lwd=4)
  legend("topleft", c(paste(n[1:k], "taxa"), "23 (anoles)"), col=c(1:k,"purple"), pch=16) 

}


plot_shape <- function(){
  plot(1,1, type='n', xlim=c(min(alpha), max(alpha)), ylim = c(0,1), main="Power by tree topology", log="x", xlab="alpha", ylab="power")
    k <- length(lambda)
  for(i in 1:length(n)){
    points(alpha, shape[[i]]$power, pch=16, col=i)
    lines(alpha, shape[[i]]$power,  col=i)
  }
  points(alpha, anoles$power, pch=16, col="purple")
  lines(alpha, anoles$power, col="purple", lwd=4)
  legend("topleft", c(paste(lambda, "lambda"), "anoles"), col=c(1:k, "purple"), pch=16  ) 
}

png("powercurve_size.png")
plot_size()
dev.off()

png("powercurve_shape.png")
plot_shape()
dev.off()

upload("powercurve_size.png", script=script, gitaddr=gitaddr, tags=tags)
upload("powercurve_shape.png", script=script, gitaddr=gitaddr, tags=tags)

