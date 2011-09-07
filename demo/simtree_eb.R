require(pmc)
## SIMULATED EXAMPLE 
simtree <- sim.bd.taxa(n=600, numbsim=1, lambda=1, mu=0, frac=1, complete=FALSE, stochsampling=FALSE)[[1]][[1]] 

exptree <- exponentialchangeTree(simtree, a=-2)

dat <- rTraitCont(exptree, sigma=2)
dat <- rTraitCont(simtree, sigma=2)
dat <- rnorm(length(simtree$tip.label) )
data<-list(phy=simtree, data=dat)

# Okay, fit the models
bm <-  fitContinuous_object(data$phy, data$data)
eb <- fitContinuous_object(data$phy, data$data, model="EB")
eb[[1]]$a <- -2
bm_v_eb <- montecarlotest(bm, eb, nboot = 1000, cpu=16)
save(list=ls(), file="simtree_lambda_dist.Rdat")

o <- confidenceIntervals.pow(bm_v_eb)
# display the confidence intervals for lambda (the test model)
#o[["test"]][,"lambda"]
#o[["test"]][,"beta"] # beta = sigma^2 is name used by geiger



## Figure 1b
png("simtree_eb.png")
hist(bm_v_eb$test_par_dist[3,], col=rgb(0,0,1,.5), border="white", breaks=15, main="", xlab="Estimated a")
abline(v=eb[[1]]$a, lwd=3, lty=2, col="darkred") #True value
dev.off()

require(socialR)
upload("simtree_eb.png", script="simtree_eb.R")
