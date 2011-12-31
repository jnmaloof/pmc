## Illustrate plotting the trees 

require(ape)
require(pmc)
data(anoles)

## Convert is a custom pmc function that toggles between tree formats
ou4_tree <- convert(tree, regimes=anoles[["OU.4"]])
ou3_tree <- convert(tree, regimes=anoles[["OU.LP"]])
ou15_tree <- convert(tree, regimes=anoles[["OU.15"]])
## plot labels as a seperate panel without a visble tree
blank <- ou15_tree
blank$edge.length <- rep(0, length(blank$edge.length))
## Display the trees.  custom treepalette function for flexible coloring by regimes
par(mfrow=c(1,4), mar=c(5, 0.1, 4, 0.1))
plot(ou3_tree, edge.color = treepalette(ou3_tree), edge.width=5, cex=1.2, show.tip.label=FALSE)
mtext("(a) OU.3", cex=2)
plot(ou4_tree, edge.color = treepalette(ou4_tree), edge.width=5, cex=1.2, show.tip.label=FALSE)
mtext("(b) OU.4", cex=2)
plot(ou15_tree, edge.color = treepalette(ou15_tree), edge.width=5, cex=1.2, show.tip.label=FALSE)
mtext("(c) OU.15", cex=2)
plot(blank, edge.color = "white", edge.width=5, cex=1.8, show.tip.label=TRUE, label.offset = -.9)





