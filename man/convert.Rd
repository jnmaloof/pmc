\name{convert}
\alias{convert}
\title{toggles between ouch and ape format trees}
\usage{
  convert(ot, regimes = NULL)
}
\arguments{
  \item{ot}{a phylogenetic tree in ouch or ape format}

  \item{regimes}{if given in ape format, are appended to
  phylo format as phy$regimes.  If the ouchtree is a fitted
  hansen object, regimes will automatically be imported
  from it unless other regime choice is given here.}
}
\value{
  a phylogenetic tree in the opposite format
}
\description{
  toggles between ouch and ape format trees
}

