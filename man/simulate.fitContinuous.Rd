\name{simulate.fitContinuous}
\alias{simulate.fitContinuous}
\title{simulate method}
\usage{
  \method{simulate}{fitContinuous} (object, nsim = 1,
    seed = NULL, ...)
}
\arguments{
  \item{object}{a fitContinuous object}

  \item{nsim}{number of sims (currently always 1)}

  \item{seed}{an optional seed for the simulations (not
  implemented)}

  \item{...}{additional arguments, not implemented for
  fitContinuous simulations}
}
\value{
  simulated dataset
}
\description{
  simulate method
}
\details{
  intended as an internal function, though an available S3
  method
}

