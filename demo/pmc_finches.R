# Author: Carl Boettiger <cboettig@gmail.com>
# License: BSD

require(pmc)
# load the geospiza (Darwin's finches) data set from the "geiger" package
library(geiger)
data(geospiza)

# Some commands to run this in parallel 
require(snowfall) # load the parallelization package
sfInit(parallel=T, cpu=4) # specify how many cpus we have
sfLibrary(geiger) # export any libraries we've loaded 

# running pmc takes just one line:
bm_v_lambda <- pmc(geospiza$geospiza.tree, geospiza$geospiza.data, "BM", "lambda", nboot=200)

##### Data Plots & Summaries  ##############

## All the parameter information is available in par_dists
## For instance, we can histogram the bootstrapped estimates of lambda, 
## obtained by simulating under the lambda model (model B)
lambdas <- subset(bm_v_lambda$par_dist, 
           comparison=="BB" & parameter=="lambda")
p1 <- ggplot(lambdas) + geom_histogram(aes(value)) +
      geom_vline(xintercept=bm_v_lambda$B[[1]]$lambda)
print(p1)

betas <- subset(bm_v_lambda$par_dist, comparison=="BB" & parameter=="beta")
p2 <- ggplot(betas) + geom_histogram(aes(sqrt(value))) # beta == sigma^2
print(p2)

## We can also look at the confidence intervals of that parameter
cast(lambdas, comparison ~ parameter, function(x)
     quantile(x, c(.05, .95)), value=c("lower", "upper"))

## or of all the parameters at once 
cast(bm_v_lambda$par_dist, comparison ~ parameter, function(x)
     quantile(x, c(.05, .95)), value=c("lower", "upper"))
## Use variations of the subset command above to include 
## only certain parameters or comparisons

## We can also look at these on a plot 
p3 <- plot_pars(bm_v_lambda$par_dists) 
print(p3)
## (again, subset first if desired)


## Simulated larger tree example
require(pmc)
require(TreeSim)
simtree <- sim.bd.taxa(n=281, numbsim=1, lambda=1, mu=0, frac=1, complete=FALSE, stochsampling=FALSE)[[1]][[1]] 
simdat <- rTraitCont(lambdaTree(simtree, .6), sigma=2)
bm_v_lambda <- pmc(simtree, simdat, "BM", "lambda", nboot=200)
lambdas <- subset(bm_v_lambda$par_dist, 
           comparison=="BB" & parameter=="lambda")
p4 <- ggplot(lambdas) + geom_histogram(aes(value)) +
      geom_vline(xintercept=bm_v_lambda$B[[1]]$lambda)
print(p4)


