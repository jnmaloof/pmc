require(pmc)
## SIMULATED EXAMPLE 
simtree <- sim.bd.taxa(n=600, numbsim=1, lambda=1, mu=0, frac=1, complete=FALSE, stochsampling=FALSE)[[1]] 
## Set simulated tree to acutally look like lambda~.6 data.  
dat <- rTraitCont(lambdaTree(simtree, .6), sigma=2)
 dat <- rTraitCont(simtree, sigma=2)
 dat <- rnorm(length(simtree$tip.label) )
data<-list(phy=simtree, data=dat)

# Okay, fit the models
bm <-  fitContinuous_object(data$phy, data$data)
lambda <- fitContinuous_object(data$phy, data$data, model="lambda")
lambda[[1]]$lambda <- 0.6
bm_v_lambda <- montecarlotest(bm, lambda, nboot = 1000, cpu=16)
save(list=ls(), file="simtree_lambda_dist.Rdat")

o <- confidenceIntervals.pow(bm_v_lambda)
# display the confidence intervals for lambda (the test model)
o[["test"]][,"lambda"]
o[["test"]][,"beta"] # beta = sigma^2 is name used by geiger



## Figure 1b
cairo_pdf("simtree_lambda.pdf", width=4, height=4)
hist(bm_v_lambda$test_par_dist[3,], col=rgb(0,0,1,.5), border="white", breaks=15, main="", xlab="Estimated lambda")
abline(v=lambda[[1]][3], lwd=3, lty=2, col="darkred") #True value
text(lambda[[1]][3], 300, "True lambda", pos=2)
dev.off()


