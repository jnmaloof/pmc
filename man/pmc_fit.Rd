\name{pmc_fit}
\alias{pmc_fit}
\title{Fit any model used in PMC}
\usage{
  pmc_fit(tree, data, model, options = list())
}
\arguments{
  \item{tree}{a phylogenetic tree. can be ouch or ape
  format}

  \item{data}{trait data in ape or ouch format}

  \item{options}{whatever additional options would be
  provided to the model fit, see details}
}
\value{
  a pmc_model object, anything that has methods "simulate",
  "update", getParameters, and getLikelihood
}
\description{
  Fit any model used in PMC
}
\details{
  options should include all parameters required by the fit
  method Currently methods avialable are fitContinuous (see
  geiger package) and hansen (see the ouch package).
}
\examples{
## a geiger example
require(geiger) # just to load the data
data(geospiza)
attach(geospiza)
lambdaFit<-pmc_fit(geospiza.tree, geospiza.data, model="lambda")
Or a single trait at a time
lambdaFit<-pmc_fit(geospiza.tree, geospiza.data[1], model="lambda")
## an ouch example
require(ouch) # just for the data,
data(bimac)
tree <- with(bimac,ouchtree(node,ancestor,time/max(time),species))
ou.2 <- pmc_fit(data=log(bimac['size']),tree, model="hansen",
                list(regimes=bimac['OU.2'],sqrt.alpha=1,sigma=1))
}

