\name{montecarlotest}
\alias{montecarlotest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Run the Phylogenetic Monte Carlo Test
}
\description{
}
\usage{
montecarlotest(null, test, nboot = 100, cpu = 2, threshold = 0.95, GetParNames = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{null}{
%%     ~~Describe \code{null} here~~
}
  \item{test}{
%%     ~~Describe \code{test} here~~
}
  \item{nboot}{
%%     ~~Describe \code{nboot} here~~
}
  \item{cpu}{
%%     ~~Describe \code{cpu} here~~
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
}
  \item{GetParNames}{
%%     ~~Describe \code{GetParNames} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (null, test, nboot = 100, cpu = 2, threshold = 0.95, 
    GetParNames = TRUE) 
{
    if (cpu > 1 & !sfIsRunning()) {
        sfInit(parallel = TRUE, cpu = cpu)
        sfLibrary(pmc)
        sfExportAll()
    }
    else if (cpu < 2 & !sfIsRunning()) {
        sfInit()
    }
    else {
    }
    null_sim <- sfSapply(1:nboot, function(i) {
        data <- simulate(null)
        if (is(data, "list")) 
            data <- data$rep.1
        null <- update(null, data)
        test <- update(test, data)
        lr <- -2 * (loglik(null) - loglik(test))
        list(lr, getParameters(null), getParameters(test))
    })
    test_sim <- sfSapply(1:nboot, function(i) {
        data <- simulate(test)
        if (is(data, "list")) 
            data <- data$rep.1
        null <- update(null, data)
        test <- update(test, data)
        lr <- -2 * (loglik(null) - loglik(test))
        list(lr, getParameters(null), getParameters(test))
    })
    null_dist <- unlist(null_sim[1, ])
    test_dist <- unlist(test_sim[1, ])
    null_bootstrap_pars <- matrix(unlist(null_sim[2, ]), ncol = nboot)
    if (GetParNames) 
        rownames(null_bootstrap_pars) <- names(getParameters(null))
    test_bootstrap_pars <- matrix(unlist(test_sim[3, ]), ncol = nboot)
    if (GetParNames) 
        rownames(test_bootstrap_pars) <- names(getParameters(test))
    null_sim_test_pars <- matrix(unlist(null_sim[3, ]), ncol = nboot)
    if (GetParNames) 
        rownames(null_sim_test_pars) <- names(getParameters(test))
    test_sim_null_pars <- matrix(unlist(test_sim[2, ]), ncol = nboot)
    if (GetParNames) 
        rownames(test_sim_null_pars) <- names(getParameters(null))
    threshold_tail <- sort(null_dist)[round(threshold * nboot)]
    power <- sum(test_dist > threshold_tail)/nboot
    lr <- -2 * (loglik(null) - loglik(test))
    p <- 1 - sum(null_dist < lr)/nboot
    print(paste("power is ", power, ", p = ", p))
    reverse_p <- sum(test_dist < lr)/nboot
    reverse_threshold_tail <- sort(-test_dist)[round(threshold * 
        nboot)]
    reverse_power <- sum(-null_dist > reverse_threshold_tail)/nboot
    print(paste("reverse test power ", reverse_power, ", reverse test p = ", 
        reverse_p))
    output <- list(power = power, null_dist = null_dist, test_dist = test_dist, 
        nboot = nboot, threshold = threshold, null = null, test = test, 
        p = p, lr = lr, reverse_p = reverse_p, reverse_power = reverse_power, 
        null_par_dist = null_bootstrap_pars, test_par_dist = test_bootstrap_pars, 
        null_sim_test_pars = null_sim_test_pars, test_sim_null_pars = test_sim_null_pars)
    class(output) <- "pow"
    output
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
