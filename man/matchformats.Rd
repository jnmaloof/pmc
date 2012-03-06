\name{matchformats}
\alias{matchformats}
\title{Internal helper function to match data formats}
\usage{
  matchformats(A, B, sim)
}
\arguments{
  \item{A}{model with simulate & update methods}

  \item{B}{another model with simulate & update methods}

  \item{sim}{simulations produced by model A}
}
\value{
  a list of two data sets.  The first is in a format
  appropriate to update model A, the second is in the
  format approriate to update model B
}
\description{
  Internal helper function to match data formats
}

