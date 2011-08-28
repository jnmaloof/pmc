require(pmc)
data(geospiza)


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


# parallel by hand, annoying but more flexible
require(snowfall)
sfInit(parallel=TRUE, cpu=16)
sfExportAll()
sfLibrary(pmc)
sfLibrary(geiger)


nboot <- 45
Asim <- sfLapply(1:nboot, function(i) compare_models(bm, lambda))
Bsim <- sfLapply(1:nboot, function(i) compare_models(lambda, bm))
bm_v_lambda <- collect(Asim, Bsim, bm, lambda)

png("test.png")
plot(bm_v_lambda)
dev.off()

require(socialR)
upload("test.png", script="example.R", gitaddr=gitcommit(script="example.R"))

## parallel automatically ## 
mc <- montecarlotest(bm, lambda, nboot = nboot, cpu=16)
png("mc.png")
plot(mc)
dev.off()

upload("mc.png", script="example.R", gitaddr=gitcommit(script="example.R"))



