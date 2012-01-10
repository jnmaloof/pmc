# generics.R
# Author: Carl Boettiger <cboettig@gmail.com>
# License: BSD

#' S3 generic to provide a log-likelihood
#' @param x a model fit that can return a loglikelihood
#' @keywords internal
loglik <- function(x, ...) UseMethod("loglik")

#' S3 generic to extract the phylogeny  
#' @param x a model fit 
#' @keywords internal
get_phy <- function(x, ...) UseMethod("get_phy")

#' S3 generic to grab trait data
#' @param x a model fit 
#' @keywords internal
get_data <- function(x, ...) UseMethod("get_data")

#' S3 generic to provide the parameters 
#' @param x a model fit of interest 
#' @keywords internal
getParameters <- function(x, ...) UseMethod("getParameters")



