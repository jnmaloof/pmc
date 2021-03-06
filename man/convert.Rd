\name{convert}
\alias{convert}
\title{toggles between ouch and ape format trees}
\usage{
  convert(ot, regimes = NULL, safe = TRUE)
}
\arguments{
  \item{ot}{a phylogenetic tree in ouch or ape format}

  \item{regimes}{if given in ape format, are appended to
  phylo format as phy$regimes.  If the ouchtree is a fitted
  hansen object, regimes will automatically be imported
  from it unless other regime choice is given here.}

  \item{safe}{mode for going from ape to ouch.  After
  converting, this writes to a temporary nexus file and
  reads the tree back in, because phylo format does not
  have a unique specification for a unique tree, resulting
  in all kinds of silly problems when developers haven't
  been careful. Defaults to true, and will clean up after
  itself.}
}
\value{
  a phylogenetic tree in the opposite format
}
\description{
  toggles between ouch and ape format trees
}

