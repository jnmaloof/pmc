require(pmc)
data(anoles)

# pmc is parallel-ready using the snowfall package, which allows parallelization on multicore 
# processors (using SOCK method) and also across entire clusters using MPI or other methods.
# Here is how we'd enable it on a quad-core machine
require(snowfall)             # load package
sfInit(parallel=T, cpu=4)  # create the appropriate parrallelization mode
sfExportAll()                 # make all data available to all parallel nodes
#sfLibrary(pmc)                # also need to export the library

### Run all the pairwise Monte Carlo Likelihood tests
ou3v4 <- pmc(tree, log(anoles["size"]), modelA="hansen", modelB="hansen", 
             optionsA=list(regimes=anoles["OU.LP"], sqrt.alpha=1, sigma=1), 
             optionsB=list(regimes=anoles["OU.4"], sqrt.alpha=1, sigma=1), nboot=20)

ou3v15 <- pmc(tree, log(anoles["size"]), "hansen", "hansen", 
             list(regimes=anoles["OU.LP"], sqrt.alpha=1, sigma=1), 
             list(regimes=anoles["OU.15"], sqrt.alpha=1, sigma=1), nboot=20)
                   
ou1v3 <- pmc(tree, log(anoles["size"]), "hansen", "hansen", 
             list(regimes=anoles["OU.1"], sqrt.alpha=1, sigma=1), 
             list(regimes=anoles["OU.4"], sqrt.alpha=1, sigma=1), nboot=20)
 
ou0v1 <- pmc(tree, log(anoles["size"]), "brown", "hansen", 
             list(), 
             list(regimes=anoles["OU.1"], sqrt.alpha=1, sigma=1), nboot=20)


######################################
##             Plots               ###
######################################
a <- plot(ou3v4, A="OU.3", B="OU.4") 
b <- plot(ou3v15, A="OU.3", B="OU.15") 
c <- plot(ou1v3, A="OU.1", B="OU.3") 
d <- plot(ou0v1, A="BM", B="OU.1") 

# display on a grid of plots 
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y) 
print(a, vp = vplayout(1, 1)) 
print(b, vp = vplayout(1, 2)) 
print(c, vp = vplayout(2, 1)) 
print(d, vp = vplayout(2, 2))

## Illustrate plotting the trees 

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





