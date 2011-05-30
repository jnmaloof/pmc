#object_oriented_laser.R

#	Note that ape also has a fitting function called birthdeath
#	Laser 
#

fit_speciation <- function(tree, model=c("bd", "pureBirth", "birthdeath"), ...){
	model=match.arg(model)
	if(model=="bd"){ 
		fit <- bd(branching.times(tree), ...)
		b = fit$r/(1-fit$a)
		d = fit$r*fit$a/(1-fit$a)
		out <- list(tree = tree, model=model, loglik=fit$LH, b=b, d=d, r=fit$r, a=fit$a, aic=fit$aic, k=2)
	} else if (model =="pureBirth"){
		fit <- pureBirth(branching.times(tree), ...)
		out <- list(tree = tree, model=model, loglik=fit$LH, b=fit$r1, d=0, aic=fit$aic, k=1)
	} else if (model == "birthdeath"){ # ape fn
		fit <- birthdeath(tree)
		fit$a <- fit$para[1]
		fit$r <- fit$para[2]
		b = fit$r/(1-fit$a)
		d = fit$r*fit$a/(1-fit$a)
		out <- list(tree = tree, model=model, loglik=-fit$dev/2, b=b, d=d, r=fit$r, a=fit$a, k=1)
	}
	class(out) <- "speciation"
	out
}

simulate.speciation <- function(fit){
#	birthdeathSim(b=fit$b, d=fit$d, CladeSize=fit$tree$Nnode+1)i
	sim.bd.taxa(n=fit$tree$Nnode+1, numbsim=1, lambda=fit$b, mu=fit$d, frac=1, complete=FALSE, stochsampling=FALSE)[[1]]
}
update.speciation <- function(fit, tree){
	fit <- fit_speciation(tree, model=fit$model)
}
loglik.speciation <- function(fit) fit$loglik

getParameters.speciation <- function(fit){
	if(fit$model=="bd") out <- c(llik=fit$loglik, b=fit$b, d=fit$d, r=fit$r, a=fit$a)
	if(fit$model=="pureBirth") out <- c(llik=fit$loglik, b=fit$b)
	out
}

