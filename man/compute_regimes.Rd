\name{compute_regimes}
\alias{compute_regimes}
\title{internal method for computing the regimes used by format_data()}
\usage{
  compute_regimes(tree, traits, species_names, regimes)
}
\arguments{
  \item{tree}{ouch-formated phylogenetic tree}

  \item{traits}{ouch-formated data set}

  \item{species_names}{for all the nodes}

  \item{regimes}{list in ouch format}
}
\value{
  a list of the regimes as edge-labels for ape's phylo
  format
}
\description{
  internal method for computing the regimes used by
  format_data()
}
\keyword{internal}

