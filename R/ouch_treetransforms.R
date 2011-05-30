
# Need a fast version of this function which can be called 
# directly on an ouch tree.
scaleTree <- function(tree, model, par){
	if(model == "BM") out <- tree
	else if(model == "OU") out <- ouTree(tree, par)
	else if(model == "lambda") out <- lambdaTree(tree, par)
	else if(model == "kappa") out <- kappaTree(tree, par)
	else if(model == "delta") out <- deltaTree(tree, par)
	else print(paste("Transform model ", model, " not recognized"))
	out
}


fitcts <- function(tree, data, model){
	require(phyloniche)
	# Takes ape tree and data formats.  
	# data must have names (not rownames) matching tips.  
	# should convert the tree to an ouchtree only once
	# and have a function to transform the ouchtree.  
	scaletree_loglik <- 
	function(par){
		tree <- scaleTree(tree, model, par)
		input <- ape2ouch_all(tree, data)
		bm <- brown(input$data, input$tree)
		-bm@loglik
	}
	o <- optimize(f = scaletree_loglik, interval=c(0,1))
	bm <- brown(input$data, input$tree)
	output <- list(tree=tree, data=data, model=model, scalepar=o$minimum, loglik=-o$objective, sigma=bm@sigma, browntree=bm, root=bm@theta[[1]])
	class(output) <- "fitcts"
	output
}

simulate.fitcts <- function(fit){
		tree <- transformTree(fit) 
	if(fit$model != "OU")
		data <- rTraitCont(tree, model="BM",
							sigma=fit[[1]]$beta, 
							root=fit$root)
	if(fit$model == "OU")
		data <- rTraitCont(fit$tree, model="OU", 
							sigma=fit[[1]]$beta,
							alpha=fit[[1]]$alpha,
							theta=fit$root, root=fit$root)
	data
}
update.fitcts <- function(fit, data){
	fitcts(fit$tree, data, fit$model)
}

loglik.fitcts <- function(fit) fit$loglik



