\name{treepalette}
\alias{treepalette}
\title{flexible regime coloring for trees}
\usage{
  treepalette(tree, regimes = NULL,
    colormap = c("rainbow", "heat.colors", "terrain.colors", "topo.colors", "cm.colors", "gray"),
    custom = NULL, rev = FALSE)
}
\arguments{
  \item{tree}{A phylo class or ouch-class tree}

  \item{regimes}{in ouch format, (only needed if not
  already given in tree)}

  \item{colormap}{a standard colormap color, see the list}

  \item{custom}{a custom colormap}

  \item{rev}{logical, reverse the order of the colormap
  colors}
}
\value{
  a set of colors than can be passed as edge.color in
  plot.phylo
}
\description{
  flexible regime coloring for trees
}

