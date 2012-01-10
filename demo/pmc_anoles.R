# Author: Carl Boettiger <cboettig@gmail.com>
# License: BSD


require(pmc)
data(anoles)



#############################################
## Optional use of parallel computing #######
#############################################
require(snowfall)             # load package
sfInit(parallel=T, cpu=4)  #  parellel for 4 core, shared.  See snowfall manual for more.  
sfExportAll()                 # make all data available to all parallel nodes


#########################################################
### Run all the pairwise Monte Carlo Likelihood tests ###
#########################################################

ou3v4 <- pmc(tree, log(anoles["size"]), modelA="hansen", modelB="hansen", 
             optionsA=list(regimes=anoles["OU.LP"], sqrt.alpha=1, sigma=1), 
             optionsB=list(regimes=anoles["OU.4"], sqrt.alpha=1, sigma=1),
             nboot=20)

ou3v15 <- pmc(tree, log(anoles["size"]), "hansen", "hansen", 
             list(regimes=anoles["OU.LP"], sqrt.alpha=1, sigma=1), 
             list(regimes=anoles["OU.15"], sqrt.alpha=1, sigma=1),
             nboot=20)
                   
ou1v3 <- pmc(tree, log(anoles["size"]), "hansen", "hansen", 
             list(regimes=anoles["OU.1"], sqrt.alpha=1, sigma=1), 
             list(regimes=anoles["OU.LP"], sqrt.alpha=1, sigma=1),
             nboot=20)
 
ou0v1 <- pmc(tree, log(anoles["size"]), "brown", "hansen", 
             list(), 
             list(regimes=anoles["OU.1"], sqrt.alpha=1, sigma=1),
             nboot=20)



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


