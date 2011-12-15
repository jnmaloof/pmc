#confidenceIntervals.R

ci <- function(dist, thresh=.95){
  n <- length(dist)
  lower <- (1-thresh)/2
  upper <- thresh+lower
  sort(dist)[c(round(lower*n), round(n*upper))]
}


ci_all_pars <-  function(par_dists, thresh=.95){ 
  sapply(1:dim(par_dists)[1], function(i) ci(par_dists[i,], thresh))
}


confidenceIntervals.pow <- function(pow, thresh=.95){
  test_ci <- ci_all_pars(pow$test_par_dist) 
  null_ci <- ci_all_pars(pow$null_par_dist)
  colnames(test_ci) <- names(getParameters(pow$test))
  rownames(test_ci) <- c(paste("lower", thresh, ":"), paste("upper", thresh, ":"))
  colnames(null_ci) <- names(getParameters(pow$null))
  rownames(null_ci) <-  c(paste("lower", thresh, ":"), paste("upper", thresh, ":"))
  
  list(test=test_ci, null=null_ci)
}


