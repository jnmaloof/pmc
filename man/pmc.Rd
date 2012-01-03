\name{pmc}
\alias{pmc}
\title{Performs a phylogenetic monte carlo between modelA and modelB}
\usage{
  pmc(tree, data,
    modelA = c("BM", "OU", "lambda", "kappa", "delta", "EB", "white", "trend", "hansen", "brown"),
    modelB = c("BM", "OU", "lambda", "kappa", "delta", "EB", "white", "trend", "hansen", "brown"),
    optionsA = list(), optionsB = list(), nboot = 20, ...)
}
\arguments{
  \item{tree}{A phylogenetic tree.  Can be phylo (ape) or
  ouch tree}

  \item{data}{The data matrix}

  \item{modelA}{a model from the list, or a custom model,
  see details}

  \item{modelB}{any other model from the list, or custom
  model, see details}

  \item{optionsA}{a named list of arguments for fitting
  model A}

  \item{optionsB}{a named list of arguments for fitting
  model B}

  \item{nboot}{number of bootstrap replicates to use}

  \item{...}{additional arguments}
}
\value{
  list with the nboot likelihood ratios obtained from
  fitting both models to data simulated by model A, and the
  nboot likelihood ratios obtained by fitting both models
  to simulations from model B, and the likelihood ratio
  between the original MLE estimated models from the data.
  return has object class pmc.
}
\description{
  Simulates data under each model and returns the
  distribution of likelihood ratio, L(B)/L(A), under for
  both simulated datasets.
}
\details{
  Possible models are all models from fitContinuous & ouch
  Currently mixing models isn't supported.
}
\examples{
require(geiger) # just for the sample data
data(geospiza)
attach(geospiza)
out <- pmc(geospiza.tree, geospiza.data[1], "BM", "lambda", nboot=5)
plot_pars(out) # show the parameters
## Ex. mixing methods from packages -- data formats handled automatically
## Load Libraries ##
require(TreeSim) # to simulate a sample phyologeny
require(snowfall) # optional for parallelization

## Simulate Data ##
simtree <- sim.bd.taxa(n=60, numbsim=1, lambda=1, mu=0, frac=1,
                       complete=FALSE, stochsampling=FALSE)[[1]][[1]]
# simulates the EB model (by tranforming the tree)
dat <- rTraitCont(exponentialchangeTree(simtree, a=-0.5), sigma=5)

## run PMC ##
sfInit(parallel=TRUE, cpu=2)
sfExportAll() # for parallelization
sfLibrary(geiger)
## Specify the models and the options needed.
regimes <- format_data(simtree, dat)$noregimes
eb_v_ou <- pmc(simtree, dat, modelA="EB", modelB="hansen",
               optionsB=list(sqrt.alpha=1, sigma=1, regimes=regimes),
               nboot=5) # small nboot just for example
plot(eb_v_ou)
}

