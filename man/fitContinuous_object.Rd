\name{fitContinuous_object}
\alias{fitContinuous_object}
\title{
an object oriented version of the Fit Continuous function in Geiger.  
}
\description{
Works as fitContinuous, but returns an object of type "fitContinuous" that has methods associated with it to "simulate", "update", get the "loglik" or the names of its parameters "GetParameterNames".  These objects can be passed directly to the montecarlotest function.
}
\usage{
fitContinuous_object(tree, data, model = "BM", ...)
}
\arguments{
  \item{tree}{
  A phylogenetic tree of class "phylo"
}
  \item{data}{
    data, see ?fitContinuous
}
  \item{model}{
    a model type, see ?fitContinuous
}
  \item{\dots}{
    other options taken by fitContinuous
}
}

\value{ A fitContinuous S3 object with associated methods
}
\author{
  Carl Boettiger
}


\examples{
  data(geospiza)
# get the data as a single named numeric
  traits <- geospiza$geospiza.data[[1]]
  names(traits) <- rownames(geospiza$geospiza.data)
## get rid of tips we don't have data for
  data <- treedata(geospiza$geospiza.tree,traits)
  names(data$data) <- rownames(data$data)

# Okay, fit the models
  bm <-  fitContinuous_object(data$phy, data$data)
  lambda <- fitContinuous_object(data$phy, data$data, model="lambda")
  lambda[[1]]$lambda <- 0.6
  \dontrun{
    bm_v_lambda <- montecarlotest(bm, lambda, nboot = 1000, cpu=16)
  }
}

