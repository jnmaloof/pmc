\name{treepower}
\alias{treepower}
\title{function to create the power curve for a tree}
\usage{
  treepower(tree, nboot = 100, cpu = 2, threshold = 0.95,
    alpha = seq(0.001, 30, length = 100),
    method = "hansen")
}
\arguments{
  \item{tree}{an ouchtree class phylogeny}

  \item{nboot}{number of bootstraps to use in power
  simulation}

  \item{cpu}{number of cpu's available}

  \item{threshold}{to use for the power calculation}

  \item{alpha}{a sequence of alpha values to try}

  \item{method}{currently ignored as only comparisons to OU
  through hansen are allowed}
}
\value{
  a power object, see examples in the vignette
}
\description{
  compares BM against OU of increasing value of alpha
}
\details{
  function can take a while to run all comparisons.
}

