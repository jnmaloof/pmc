###############
require(socialR)
script <- "early_burst_example.R"
gitaddr <- gitcommit(script)
tags="phylogenetics"
tweet_errors(script, tags=tags)
###############

#early_burst_example.R
require(pmc)
require(TreeSim)
## SIMULATED EXAMPLE 
simtree <- sim.bd.taxa(n=60, numbsim=1, lambda=1, mu=0, frac=1, complete=FALSE, stochsampling=FALSE)[[1]][[1]] 
dat <- rTraitCont(exponentialchangeTree(simtree, a=-1.5), sigma=5)
#dat <- rTraitCont(simtree, sigma=5, alpha=10, model="OU")
data<-list(phy=simtree, data=dat)

# Okay, fit the models
#bm <-  fitContinuous_object(data$phy, data$data)
eb <- fitContinuous_object(data$phy, data$data, model="EB")
ou <- fitContinuous_object(data$phy, data$data, model="OU")



#bm_v_eb <- montecarlotest(bm, eb, nboot = 200, cpu=16)
ou_v_eb <- montecarlotest(ou, eb, nboot = 200, cpu=16)
png("eb.png")
plot(ou_v_eb)
dev.off()

save(list=ls(), file="early_burst.Rdat")

require(socialR)
upload("eb.png", script="early_burst_example.R", tags="phylogenetics")
