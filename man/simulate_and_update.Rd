\name{simulate_and_update}
\alias{simulate_and_update}
\title{Internal helper function
Simulates under model A, updates both A and B based on that data}
\usage{
  simulate_and_update(A, B)
}
\arguments{
  \item{A}{model with simulate & update methods}

  \item{B}{another model with simulate & update methods}
}
\value{
  a list of with fit of A on data simulated under A, the
  fit on B on data simulated under A, and the simulated
  data itself
}
\description{
  Internal helper function Simulates under model A, updates
  both A and B based on that data
}

