#early_burst_example.R

## Ex. mixing methods from packages -- data formats handled automatically 
## Load Libraries ##
require(pmc)
require(TreeSim) # to simulate a sample phyologeny
require(snowfall) # optional for parallelization

## Simulate Data ##
simtree <- sim.bd.taxa(n=60, numbsim=1, lambda=1, mu=0, frac=1,
                       complete=FALSE, stochsampling=FALSE)[[1]][[1]] 
# simulates the EB model (by tranforming the tree)
dat <- rTraitCont(exponentialchangeTree(simtree, a=-0.5), sigma=5)

## run PMC ## 
sfInit(par=T, cpu=2)
sfExportAll() # for parallelization
## Specify the models and the options needed.  
eb_v_ou <- pmc(simtree, dat, modelA="EB", modelB="hansen", 
               optionsB=list(sqrt.alpha=1, sigma=1, regimes=regimes), 
               nboot=2) # small nboot just for example
plot(eb_v_ou) 


## examples of more basic functions 
ouch <- format_data(simtree, dat)
ape <- format_data(ouch$tree, ouch$dat)
A <- pmc_fit(ape$tree, ape$data, model="BM")
regimes <- format_data(simtree,dat)$noregimes
B <- pmc_fit(simtree, dat, model="hansen", options=list(sqrt.alpha=1, sigma=1, regimes=regimes))



