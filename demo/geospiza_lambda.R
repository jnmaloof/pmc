require(pmc)
## GEOSPIZA EXAMPLE
data(geospiza)
# get the data as a single named numeric
traits <- geospiza$geospiza.data[[1]]
names(traits) <- rownames(geospiza$geospiza.data)
## get rid of tips we don't have data for
data <- treedata(geospiza$geospiza.tree,traits)
names(data$data) <- rownames(data$data)


# theoretical maximum lambda
C<-vcv.phylo(data$phy)
maxLambda<-max(C)/max(C[upper.tri(C)])
bounds = list(lambda=c(0,maxLambda))
print(bounds)

# Okay, fit the models
bm <-  fitContinuous_object(data$phy, data$data, bounds)
lambda <- fitContinuous_object(data$phy, data$data, model="lambda", bounds)
lambda[[1]]$lambda <- 0.6
bm_v_lambda <- montecarlotest(bm, lambda, nboot = 1000, cpu=16)

o <- confidenceIntervals.pow(bm_v_lambda)
# display the confidence intervals for lambda (the test model)
o[["test"]][,"lambda"]
o[["test"]][,"beta"] # beta = sigma^2 is name used by geiger


require(socialR)
social_plot({
hist(bm_v_lambda$test_par_dist[3,], col=rgb(0,0,1,.5), border="white", breaks=15, main="", xlab="Estimated lambda")
abline(v=lambda[[1]][3], lwd=3, lty=2, col="darkred") #True value
text(lambda[[1]][3], 300, "True lambda", pos=2)}, tags="phylogenetics", description="Fig1a with different lambda bounds")




save(list=ls(), file="geospiza_lambda.Rdat")
## FIGURE 1a
cairo_pdf("geospiza_lambda.pdf", width=4, height=4)
hist(bm_v_lambda$test_par_dist[3,], col=rgb(0,0,1,.5), border="white", breaks=15, main="", xlab="Estimated lambda")
abline(v=lambda[[1]][3], lwd=3, lty=2, col="darkred") #True value
text(lambda[[1]][3], 300, "True lambda", pos=2)
dev.off()

## Figure 1c
cairo_pdf("geospiza_sigma.pdf", width=4, height=4)
hist(bm_v_lambda$test_par_dist[2,], col=rgb(0,0,1,.5), border="white", breaks=15, main="", xlab="Estimated sigma")
abline(v=lambda[[1]][2], lwd=3, lty=2, col="darkred") #True value
text(lambda[[1]][2], 300, "True sigma", pos=4)
dev.off()







