require(pmc)
## GEOSPIZA EXAMPLE
data(geospiza)
bm_v_lambda <- pmc(geospiza$geospiza.tree, geospiza$geospiza.data, "BM", "lambda")


# get the data as a single named numeric
traits <- geospiza$geospiza.data[[1]]
names(traits) <- rownames(geospiza$geospiza.data)
## get rid of tips we don't have data for
data <- treedata(geospiza$geospiza.tree,traits)
names(data$data) <- rownames(data$data)

# theoretical maximum lambda
C <- vcv.phylo(data$phy)
maxLambda <- max(C)/max(C[upper.tri(C)])
# replicates of maxLambda > 1 won't bootstrap, result in singular matrices
bounds <- list(lambda=c(0,1))

# Okay, fit the models
bm <-  fitContinuous_object(data$phy, data$data, bounds=bounds)
lambda <- fitContinuous_object(data$phy, data$data, model="lambda", 
                               bounds=bounds)

lambda[[1]]$lambda <- 0.6
bm_v_lambda <- montecarlotest(bm, lambda, nboot = 50, cpu=2)

# display the confidence intervals for lambda (the test model)
o <- confidenceIntervals.pow(bm_v_lambda)
o[["test"]][,"lambda"]
o[["test"]][,"beta"] # beta = sigma^2 is name used by geiger

save(list=ls(), file="geospiza_lambda.Rdat")
## FIGURE 1a
hist(bm_v_lambda$test_par_dist[3,], col=rgb(0,0,1,.5),
     border="white", breaks=15, main="", xlab="Estimated lambda")
abline(v=lambda[[1]][3], lwd=3, lty=2, col="darkred") #True value
text(lambda[[1]][3], 300, "True lambda", pos=2)
## Figure 1c
hist(bm_v_lambda$test_par_dist[2,], col=rgb(0,0,1,.5),
     border="white", breaks=15, main="", xlab="Estimated sigma")
abline(v=lambda[[1]][2], lwd=3, lty=2, col="darkred") #True value
text(lambda[[1]][2], 300, "True sigma", pos=4)


