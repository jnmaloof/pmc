#blomberg_power.R

blomberg_power <- function(phy, nboot=200, cpu=2, threshold =.95){

	initBM <- function(phy, sigma2 = 10, root = 0){
		bm <- list(Trait1=list(beta=sigma2), tree=phy, root=root, model="BM")
		class(bm) <- "fitContinuous"
		bm
	}

	K_dist_bm <- sfSapply(1:nboot, function(i){
		bm <- initBM(phy)
		x <- simulate(bm)  ## simulate and estimate blomberg's stat nboot times
		out <- phylosignal(x, phy, reps=100) ## blomberg reshuffles reps times
		out$K
	})


	K_dist_null <- sfSapply(1:nboot, function(i){
		bm <- initBM(phy)
		x <- simulate(bm) ## Will replace this data with independent normals of same statistics
		y <- rnorm(length(x), mean=mean(x), sd=sd(x))
		out <- phylosignal(y, phy, reps=100)
		out$K
	})

	
	threshold_tail <- sort(K_dist_null)[ round(threshold*nboot) ]
	power <- sum(K_dist_bm > threshold_tail)/nboot

	overlap <- K_dist_bm %*% K_dist_null
	out <- list(overlap=overlap, power=power, null_dist =K_dist_null, test_dist=K_dist_bm, phy=phy)	

	plt <- function(){
		plot(density(K_dist_bm), col="darkred", lwd=3, ylim=c(0,3), main=paste("overlap = ", overlap))
		lines(density(K_dist_null), col="darkblue", lwd=3)
	}
	out
}


