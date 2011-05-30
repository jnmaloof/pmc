\name{treepalette}
\alias{treepalette}
\title{
 A simple function for autmatically coloring regimes on a phylogeny
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
treepalette(apetree, colormap = c("rainbow", "heat.colors", "terrain.colors", "topo.colors", "cm.colors"), rev = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{apetree}{
%%     ~~Describe \code{apetree} here~~
}
  \item{colormap}{
%%     ~~Describe \code{colormap} here~~
}
  \item{rev}{
%%     ~~Describe \code{rev} here~~
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
function (apetree, colormap = c("rainbow", "heat.colors", "terrain.colors", 
    "topo.colors", "cm.colors"), rev = FALSE) 
{
    colormap <- match.arg(colormap)
    if (colormap == "rainbow") 
        levels(apetree$regimes) <- rainbow(length(levels(apetree$regimes)))
    if (colormap == "heat.colors") 
        levels(apetree$regimes) <- heat.colors(length(levels(apetree$regimes)))
    if (colormap == "terrain.colors") 
        levels(apetree$regimes) <- terrain.colors(length(levels(apetree$regimes)))
    if (colormap == "topo.colors") 
        levels(apetree$regimes) <- topo.colors(length(levels(apetree$regimes)))
    if (colormap == "cm.colors") 
        levels(apetree$regimes) <- cm.colors(length(levels(apetree$regimes)))
    if (rev == TRUE) 
        levels(apetree$regimes) <- rev(levels(apetree$regimes))
    as.character(apetree$regimes)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
