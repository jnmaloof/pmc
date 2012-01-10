\name{fitContinuous_object}
\alias{fitContinuous_object}
\title{wrapper for fitContinuous}
\usage{
  fitContinuous_object(tree, data, model = "BM",
    bounds = NULL, meserr = NULL)
}
\arguments{
  \item{tree}{a phylo object phylogenetic tree}

  \item{data}{ape data phylo format}

  \item{model}{the type of model to fit}

  \item{bounds}{optional list of bounds for the parameter
  fits}

  \item{meserr}{measurment error, if any}
}
\value{
  an object of class fitContinuous
}
\description{
  Runs the fitContinuous function and returns and object
  with necessary information to run methods such as
  'update' and 'simulate'.  Call is the same as the
  original fitContinuous function.
}

