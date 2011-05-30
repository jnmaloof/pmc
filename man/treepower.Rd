\name{treepower}
\alias{treepower}

\title{
  Estimate power in a phylogenetic tree
 }
\description{
  Applies the montecarlotest function where the test model is an OU
  process, evaluated at several strenghts alpha, and the null midel is a 
  BM process. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
treepower(tree, nboot = 100, cpu = 2, threshold = 0.95, alpha = seq(1e-05, 100, length = 1000), method = "hansen")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
%%     ~~Describe \code{tree} here~~
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
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
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
function (tree, nboot = 100, cpu = 2, threshold = 0.95, alpha = seq(1e-05, 
    100, length = 1000), method = "hansen") 
{
    data <- c(rep(NA, tree@nnodes - tree@nterm), rnorm(tree@nterm))
    names(data) <- tree@nodes
    null <- brown(data, tree)
    regimes <- as.factor(rep("ns", tree@nnodes))
    names(regimes) <- tree@nodes
    test <- hansen(data, tree, regimes, 1, 1)
    if (cpu > 1) {
        sfInit(parallel = TRUE, cpu = cpu)
        sfLibrary(pmc)
        sfExportAll()
    }
    else sfInit()
    null_dist <- sfSapply(1:nboot, function(i) {
        data <- simulate(null)$rep.1
        null <- update(null, data)
        test <- update(test, data)
        -2 * (null@loglik - test@loglik)
    })
    test_dist <- sfLapply(1:length(alpha), function(i) {
        test@sqrt.alpha <- sqrt(alpha[i])
        sfSapply(1:nboot, function(i) {
            data <- simulate(test)$rep.1
            null <- update(null, data)
            test <- update(test, data)
            -2 * (null@loglik - test@loglik)
        })
    })
    power <- sapply(1:length(alpha), function(i) {
        threshold_tail <- sort(null_dist)[round(threshold * nboot)]
        sum(test_dist[[i]] > threshold_tail)/nboot
    })
    list(null_dist = null_dist, test_dist = test_dist, power = power, 
        alpha = alpha, nboot = nboot, threshold = threshold, 
        tree = tree)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
