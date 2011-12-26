#' flexible regime coloring for trees 
#' @param tree A phylo class or ouch-class tree
#' @param regimes in ouch format, (only needed if not already given in tree) 
#' @param colormap a standard colormap color, see the list
#' @param custom a custom colormap
#' @param rev logical, reverse the order of the colormap colors 
#' @return a set of colors than can be passed as edge.color in plot.phylo
#' @export
treepalette <- function(tree, regimes=NULL, colormap = c("rainbow", "heat.colors",
                        "terrain.colors", "topo.colors", "cm.colors", "gray"),
                        custom=NULL, rev=FALSE){

  if(!is.null(regimes) & is(tree,"ouchtree"))
    apetree <- convert(tree, regimes=regimes)
  else
    apetree <- tree

	colormap <- match.arg(colormap)
	if(colormap=="rainbow") 
    levels(apetree$regimes) <- rainbow(length(levels(apetree$regimes)))
	if(colormap=="heat.colors") 
    levels(apetree$regimes) <- heat.colors(length(levels(apetree$regimes)))
	if(colormap=="terrain.colors") 
    levels(apetree$regimes) <- terrain.colors(length(levels(apetree$regimes)))
	if(colormap=="topo.colors") 
    levels(apetree$regimes) <- topo.colors(length(levels(apetree$regimes)))
	if(colormap=="cm.colors") 
    levels(apetree$regimes) <- cm.colors(length(levels(apetree$regimes)))
	if(colormap=="gray"){
    n <-length(levels(apetree$regimes))+2
    levels(apetree$regimes) <- gray((1:n)/n)
  }
  if(!is.null(custom))
    levels(apetree$regimes) <- custom
	if(rev==TRUE) 
    levels(apetree$regimes) <- rev(levels(apetree$regimes))
	as.character(apetree$regimes)
} 


