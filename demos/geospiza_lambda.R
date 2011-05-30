require(pmc)
## GEOSPIZA EXAMPLE
data(geospiza)
# get the data as a single named numeric
traits <- geospiza$geospiza.data[[1]]
names(traits) <- rownames(geospiza$geospiza.data)
## get rid of tips we don't have data for
data <- treedata(geospiza$geospiza.tree,traits)
names(data$data) <- rownames(data$data)

# Okay, fit the models
bm <-  fitContinuous_object(data$phy, data$data)
lambda <- fitContinuous_object(data$phy, data$data, model="lambda")
lambda[[1]]$lambda <- 0.6
bm_v_lambda <- montecarlotest(bm, lambda, nboot = 1000, cpu=16)
save(list=ls(), file="geospiza_lambda.Rdat")


## FIGURE 1a
#cairo_pdf("geospiza_lambda.pdf", width=4, height=4)


require(socialR)

social_plot({
par(mfrow=c(1,2))
hist(bm_v_lambda$test_par_dist[3,], col=rgb(0,0,1,.5), border="white", breaks=15, main="", xlab="Estimated lambda")
abline(v=lambda[[1]][3], lwd=3, lty=2, col="darkred") #True value
text(lambda[[1]][3], 300, "True lambda", pos=2)
#dev.off()

## Figure 1c
#cairo_pdf("geospiza_sigma.pdf", width=4, height=4)
hist(bm_v_lambda$test_par_dist[2,], col=rgb(0,0,1,.5), border="white", breaks=15, main="", xlab="Estimated sigma")
abline(v=lambda[[1]][2], lwd=3, lty=2, col="darkred") #True value
text(lambda[[1]][2], 300, "True sigma", pos=4)
#dev.off()
}, tag="phylogenetics")


