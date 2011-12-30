#early_burst_example.R
require(pmc)
require(TreeSim)
## SIMULATED EXAMPLE 
simtree <- sim.bd.taxa(n=60, numbsim=1, lambda=1, mu=0, frac=1, complete=FALSE, stochsampling=FALSE)[[1]][[1]] 
# simulates the EB model (by tranforming the tree)
dat <- rTraitCont(exponentialchangeTree(simtree, a=-0.5), sigma=5)
## run PMC
#ou_v_eb <- pmc(simtree, dat, modelA="EB", modelB="OU", nboot=2)
#plot(ou_v_eb)

ouch <- format_data(simtree, dat)
ape <- format_data(ouch$tree, ouch$dat)
A <- pmc_fit(simtree, dat, model="BM")

regimes <- format_data(simtree,dat)$noregimes
B <- pmc_fit(simtree, dat, model="hansen", options=list(sqrt.alpha=1, sigma=1, regimes=regimes))
#eb_v_ou <- pmc(simtree, dat, modelA="EB", modelB="hansen", optionsB=list(sqrt.alpha=1, sigma=1, regimes=regimes), nboot=10)
#plot(eb_v_ou)

