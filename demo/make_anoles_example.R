
############################################################################
## Modify bimaculus example from OUCH to include two more regime options ###
############################################################################
require(ouch)
data(bimac)			# ouch package Anolis sizes (from N. Lesser Antilles)
tree <- with(bimac,ouchtree(node,ancestor,time/max(time),species))
## Full names please 
species <- as.character(bimac[["species"]])
species[23:45] <- c("pogus, St. Maarten", "schwartzi, St. Eustatius", "schwartzi, St. Christopher", "schwartzi, Nevis", "wattsi, Barbuda", "wattsi, Antigua", "bimaculatus, St. Eustatius", "bimaculatus, Nevis", "bimaculatus, St. Christopher", "leachi, Barbuda", "leachi, Antigua", "nubilus, Redonda", "sabanus, Saba", "gingivinus, St. Barthelemy", "gingivinus, Anguilla", "gingivinus, St. Maarten", "oculatus, Dominica", "ferreus, Marie Galante", "lividus, Monserrat", "marmoratus, Guadeloupe", "marmoratus, Desirade", "terraealtae, Illes de Saintes-1", "terraealtae, Illes de Saintes-2")
tree@nodelabels <- species

## 15 different regimes
regimes <- bimac['OU.4']
regimes <- as.character(regimes[,])
regimes[35:45] <- as.character(35:45)
regimes <- as.factor(regimes)
names(regimes) <- rownames(bimac['OU.4'])
OU.15 <- regimes 

# Set up and fit a  model where each of the four major clades have their own optimum 
regimes <- bimac['OU.LP']
regimes <- as.character(regimes[,])
regimes[34:38] <- "fourth"
regimes[13:16] <- "fourth"
regimes <- as.factor(regimes)
names(regimes) <- rownames(bimac['OU.LP'])
OU.4 <- regimes 

anoles <- data.frame(bimac[["size"]], bimac[["node"]], bimac[["ancestor"]], 
            bimac[["time"]], bimac[["species"]],
            bimac[["OU.1"]], bimac[["OU.LP"]], OU.4,OU.15)
colnames(anoles) <- c("size", "node", "ancestor", "time", "species", 
                      "OU.1", "OU.LP", "OU.4", "OU.15")
rownames(anoles) <- rownames(bimac)
write.csv(anoles, file="anoles.csv")

save(list=c("tree", "anoles"), file="anoles.rda")
