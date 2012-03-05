## powertest.R
# power curve
#' function to create the power curve for a tree 
#'
#' compares BM against OU of increasing value of alpha
#' @param tree an ouchtree class phylogeny
#' @param nboot number of bootstraps to use in power simulation
#' @param cpu number of cpu's available
#' @param threshold to use for the power calculation
#' @param alpha a sequence of alpha values to try
#' @param method currently ignored as only comparisons to OU through hansen are allowed
#' @return a power object, see examples in the vignette
#' @details function can take a while to run all comparisons.  
#' @export
treepower <- function(tree, nboot = 100, cpu = 2, threshold = .95, alpha = seq(0.001, 30, length=100), method="hansen"){
	## Gotta get templates for the models, do so by fitting some dummy data
	data <- c(rep(NA, tree@nnodes-tree@nterm), rnorm(tree@nterm))
	names(data) <- tree@nodes
	null <- brown(data,tree)
	regimes <- as.factor(rep("ns", tree@nnodes))
	names(regimes) <- tree@nodes
	test <- hansen(data,tree,regimes,1,1)

	## are we in parallel?
	if(cpu>1){ 	
		sfInit(parallel=TRUE, cpus=cpu) 
		sfLibrary(pmc)
		sfExportAll()
	} else sfInit()

	null_dist <- sfSapply(1:nboot, function(i){
		data <- simulate(null)$rep.1
		null <- update(null, data)
		test <- update(test, data)
		-2*(null@loglik - test@loglik) 
	})

	## Actually do the bootstraps for each alpha
	test_dist <- sfLapply(1:length(alpha), function(i){
		test@sqrt.alpha <- sqrt(alpha[i])
		sapply(1:nboot, function(i){
				data <- simulate(test)$rep.1
				null <- update(null, data)
				test <- update(test, data)
				-2*(null@loglik - test@loglik) 
		})
	})

	## Power calculation
	power <- sapply(1:length(alpha), function(i){
		threshold_tail <- sort(null_dist)[ round(threshold*nboot) ]
		sum(test_dist[[i]] > threshold_tail)/nboot
	})

	## format the output
    list(null_dist=null_dist, test_dist=test_dist, power=power, alpha=alpha, nboot=nboot, threshold=threshold, tree=tree)
}


