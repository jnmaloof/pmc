\name{plot.pow}
\alias{plot.pow}
\title{
	Plots the output of the Monte Carlo simulation.  Lots of options to aid with visualizing and interpreting the results
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.pow(pow, main = "", legend = FALSE, type = "density", test_dist = TRUE, shade_power = FALSE, shade_p = FALSE, show_aic = FALSE, show_data = TRUE, shade = TRUE, shade_aic = FALSE, print_text = TRUE, show_text = c("p", "power", "reverse_p", "aic"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pow}{
%%     ~~Describe \code{pow} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
}
  \item{legend}{
%%     ~~Describe \code{legend} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{test_dist}{
%%     ~~Describe \code{test_dist} here~~
}
  \item{shade_power}{
%%     ~~Describe \code{shade_power} here~~
}
  \item{shade_p}{
%%     ~~Describe \code{shade_p} here~~
}
  \item{show_aic}{
%%     ~~Describe \code{show_aic} here~~
}
  \item{show_data}{
%%     ~~Describe \code{show_data} here~~
}
  \item{shade}{
%%     ~~Describe \code{shade} here~~
}
  \item{shade_aic}{
%%     ~~Describe \code{shade_aic} here~~
}
  \item{print_text}{
%%     ~~Describe \code{print_text} here~~
}
  \item{show_text}{
%%     ~~Describe \code{show_text} here~~
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
function(pow, main="", legend=FALSE, type="density", test_dist=TRUE, shade_power=FALSE, shade_p=FALSE, show_aic=FALSE, show_data=TRUE, shade=TRUE, shade_aic=FALSE, print_text=TRUE, show_text = c("p", "power", "reverse_p", "aic")){
	xlim <- c( min(pow$null_dist, pow$test_dist), max(pow$null_dist, pow$test_dist) )
	aic_line <- pow$aic_diff
	if(type != "hist"){
		nd <- density(pow$null_dist)
		td <- density(pow$test_dist)
		plot(nd, xlim=xlim, main=main, lwd=1, xlab=" Likelihood Ratio statistic")

		if(shade_p){
			shade_p <- which(nd$x > pow$lr)
			polygon(c(pow$lr,nd$x[shade_p]), c(0,nd$y[shade_p]), col=rgb(0,0,1,.5), border=rgb(0,0,1,.5))
		} else if(shade){
			polygon(nd$x, nd$y, col=rgb(0,0,1,.5), border=rgb(0,0,1,.5))
		}
		if(test_dist){
			lines(td, lwd=1)
			threshold_tail <- sort(pow$null_dist)[ round(pow$threshold*pow$nboot) ]
			if(shade_power){
				shade_power <- which(td$x > threshold_tail)
				polygon(c(threshold_tail, td$x[shade_power]), c(0,td$y[shade_power]), col=rgb(1,0,0,.5), border=rgb(1,0,0,.5))
			} else if(shade){
				polygon(td$x, td$y, col=rgb(1,0,0,.5), border=rgb(1,0,0,.5))
			}
		}
		if(shade_aic){
			shade_p <- which(nd$x > aic_line)
			polygon(c(aic_line,nd$x[shade_p]), c(0,nd$y[shade_p]), col=rgb(1,0,1,.5), border=rgb(1,0,1,.5))
		}
	} else {
	hist(pow$null_dist, xlim=xlim, lwd=3, breaks=20, col=rgb(1,0,0,.5), border="white", main=main, xlab=" Likelihood Ratio statistic")
		if(test_dist){
			hist(pow$test_dist, add=T, lwd=0, breaks=20, col=rgb(0,0,1,.5), border="white", main=main, xlab=" Likelihood Ratio statistic")
		}
	}
	if(show_data) abline(v=pow$lr, lwd=3, col="darkred", lty=2 )
	if(show_aic) abline(v=aic_line, lwd=3, col="darkblue", lty=1) 
	p <- pow$p
	aic_wrong <- sum(pow$null_dist > aic_line)/pow$nboot
	if(aic_line < pow$lr){ print("AIC prefers test model")} else { print("AIC prefers null model") }	
	print(paste("p = ",  p, ", power = ", pow$power))
	if(print_text){
		if(! is.na(match("p", show_text)) ) text(.85*par()$xaxp[2], .95*par()$yaxp[2], paste("p = ", p))
		if(! is.na(match("power", show_text)) ) text(.85*par()$xaxp[2], .9*par()$yaxp[2], paste("Power = ", pow$power))
		if(! is.na(match("aic", show_text)) ) text(.85*par()$xaxp[2], .85*par()$yaxp[2], paste("AIC wrong in ", aic_wrong*100, "\% of sims"))
		if(! is.na(match("reverse_p", show_text)) ) text(.85*par()$xaxp[2], .8*par()$yaxp[2], paste("p in reverse test ", pow$reverse_p))
	}
	if(legend) 
		legend("topleft", c("sim under Null", "sim under Test", "observed"), lty=c(1,2, 3), col=c("black", "black", "darkred"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
