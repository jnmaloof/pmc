require(pmc)

nboot = 2000
cpus = 16

####### Anoles example ########	
data(bimac)			# ouch package Anolis sizes (from N. Lesser Antilles)
tree <- with(bimac,ouchtree(node,ancestor,time/max(time),species))

## The three models from Butler and King
bm <- brown(log(bimac['size']),tree)
ou1 <- hansen(log(bimac['size']),tree,bimac['OU.1'],1,1)
ouLP <- hansen(log(bimac['size']),tree,bimac['OU.LP'],1,1, reltol=1e-5)


## 15 different regimes
regimes <- bimac['OU.4']
regimes <- as.character(regimes[,])
regimes[35:45] <- as.character(35:45)
regimes <- as.factor(regimes)
names(regimes) <- rownames(bimac['OU.4'])
half <- hansen(log(bimac['size']),tree,regimes,1,1)

# Set up and fit a  model where each of the four major clades have their own optimum 
regimes <- bimac['OU.LP']
regimes <- as.character(regimes[,])
regimes[34:38] <- "fourth"
regimes[13:16] <- "fourth"
regimes <- as.factor(regimes)
names(regimes) <- rownames(bimac['OU.LP'])
ouLP4 <- hansen(log(bimac['size']),tree,regimes,1,1, method="L")


# Plot the trees of the model
## actual node labels 
full_labels <- bm@nodelabels
full_labels[23:45] <- c("pogus, St. Maarten", "schwartzi, St. Eustatius", "schwartzi, St. Christopher", "schwartzi, Nevis", "wattsi, Barbuda", "wattsi, Antigua", "bimaculatus, St. Eustatius", "bimaculatus, Nevis", "bimaculatus, St. Christopher", "leachi, Barbuda", "leachi, Antigua", "nubilus, Redonda", "sabanus, Saba", "gingivinus, St. Barthelemy", "gingivinus, Anguilla", "gingivinus, St. Maarten", "oculatus, Dominica", "ferreus, Marie Galante", "lividus, Monserrat", "marmoratus, Guadeloupe", "marmoratus, Desirade", "terraealtae, Illes de Saintes-1", "terraealtae, Illes de Saintes-2")
ouLP4@nodelabels <- full_labels
ouLP@nodelabels <- full_labels
half@nodelabels <- full_labels


## Convert is a custom pmc function that toggles between tree formats
ouLP4_tree <- convert(ouLP4)
ouLP_tree <- convert(ouLP)
half_tree <- convert(half)
## plot labels as a seperate panel without a visble tree
blank <- half_tree
blank$edge.length <- rep(0, length(blank$edge.length))


## Display the trees.  custom treepalette function for flexible coloring by regimes
## Refuses to render properly as a cairo_pdf (raster pdf ok)
png("tree_paintings.png", width=800, height=800)
par(mfrow=c(1,4), mar=c(5, 0.1, 4, 0.1))
plot(ouLP_tree, edge.color = treepalette(ouLP_tree), edge.width=5, cex=1.2, show.tip.label=FALSE)
text(par()$xaxp[2]/2, 1, "(a) OU.3", cex=2)
plot(ouLP4_tree, edge.color = treepalette(ouLP4_tree), edge.width=5, cex=1.2, show.tip.label=FALSE)
text(par()$xaxp[2]/2, 1, "(b) OU.4", cex=2)
plot(half_tree, edge.color = treepalette(half_tree), edge.width=5, cex=1.2, show.tip.label=FALSE)
text(par()$xaxp[2]/2, 1, "(c) OU.15", cex=2)
plot(blank, edge.color = "white", edge.width=5, cex=1.8, show.tip.label=TRUE, label.offset = -.9)
dev.off()


require(snowfall)
sfInit(parallel=TRUE, cpu=cpus)
sfExportAll()
sfLibrary(pmc)
sfLibrary(geiger)

### Run all the pairwise Monte Carlo Likelihood tests
### The actually time-consuming step
t1<- system.time(ouLP_v_ouLP4 <- montecarlotest(ouLP, ouLP4, nboot=nboot, cpu=cpus))
t2<- system.time(bm_v_ouLP  <- montecarlotest(bm, ouLP, nboot=nboot, cpu=cpus))
t3<- system.time(bm_v_ou1  <- montecarlotest(bm, ou1, nboot=nboot, cpu=cpus))
t4<- system.time(ou1_v_ouLP <-  montecarlotest(ou1, ouLP, nboot=nboot, cpu=cpus))
t5<- system.time(ouLP_v_half <-  montecarlotest(ouLP, half, nboot=nboot, cpu=cpus))

save(list=ls(), file="anoles_model_choice.Rdat")


######################### FIGURE 2 ####################################################################
#################### Plots with shading to demonstrate p and power. not sure if they are necessary ####
cairo_pdf("p.pdf", height=4, width=4)
plot(ouLP_v_ouLP4, show_text=FALSE, test=FALSE, shade_p=TRUE, shade=FALSE)
dev.off()
cairo_pdf("power.pdf", height=4, width=4)
plot(ouLP_v_ouLP4, show_text=FALSE, shade_power=TRUE, shade=FALSE)
dev.off()



############## FIGURE 3 #######################################################
## Display the trees.  custom treepalette function for flexible coloring by regimes
png("tree_paintings.png", width=1000, height=400)
#cairo_pdf("tree_paintings.pdf", width=8, height=8)
par(mfrow=c(1,4), mar=c(5, 0.1, 4, 0.1))
plot(ouLP_tree, edge.color = treepalette(ouLP_tree), edge.width=5, cex=1.2, show.tip.label=FALSE)
text(par()$xaxp[2]/2, 1, "(a) OU.3", cex=2)
plot(ouLP4_tree, edge.color = treepalette(ouLP4_tree), edge.width=5, cex=1.2, show.tip.label=FALSE)
text(par()$xaxp[2]/2, 1, "(b) OU.4", cex=2)
plot(half_tree, edge.color = treepalette(half_tree), edge.width=5, cex=1.2, show.tip.label=FALSE)
text(par()$xaxp[2]/2, 1, "(c) OU.15", cex=2)
plot(blank, edge.color = "white", edge.width=5, cex=1.8, show.tip.label=TRUE, label.offset = -.9)
dev.off()

############## FIGURE 4 ##########################
cairo_pdf("ou3_v_ou4.pdf", height=4, width=4)
plot(ouLP_v_ouLP4, show_text=FALSE)
legend("topright", c("OU.3 sims", "OU.4 sims", "obs"), pch=c(15,15,46), lty=c(0,0,2), col=c(rgb(0,0,1,.5), rgb(1,0,0,.5), "darkred"))
dev.off()

cairo_pdf("bm_v_ou3.pdf", height=4, width=4)
plot(bm_v_ouLP, show_text=FALSE)
legend("topright", c("BM sims", "OU.3 sims", "obs"), pch=c(15,15,46), lty=c(0,0,2), col=c(rgb(0,0,1,.5), rgb(1,0,0,.5), "darkred"))
dev.off()
cairo_pdf("bm_v_ou1.pdf", height=4, width=4)
plot(bm_v_ou1, show_text=FALSE)
legend("topright", c("BM sims", "OU.1 sims", "obs"), pch=c(15,15,46), lty=c(0,0,2), col=c(rgb(0,0,1,.5), rgb(1,0,0,.5), "darkred"))
dev.off()
cairo_pdf("ou1_v_ou3.pdf", height=4, width=4)
plot(ou1_v_ouLP, show_text=FALSE)
legend("topright", c("OU.1 sims", "OU.3 sims", "obs"), pch=c(15,15,46), lty=c(0,0,2), col=c(rgb(0,0,1,.5), rgb(1,0,0,.5), "darkred"))
dev.off()
cairo_pdf("ou3_v_ou15.pdf", height=4, width=4)
plot(ouLP_v_half, show_text=FALSE)
legend("topright", c("OU.3 sims", "OU.15 sims", "obs"), pch=c(15,15,46), lty=c(0,0,2), col=c(rgb(0,0,1,.5), rgb(1,0,0,.5), "darkred"))
dev.off()


######################## FIGURE 5 ##########################################################################
### Plots of the error rates ########
plot_error <- function(pow, null, test, info){
  plot(pow, shade_aic=T, show_aic=T, shade=F, show_data=F, show_text=F, info=info, legend=TRUE, main=paste(null, "vs", test))
# legend("topright", c(paste("AIC", "Type I err"), paste("AIC", "Type II err")), pch=15, col=c(rgb(1,.5,0,.5), rgb(1,1,0,.5)))
}

info <- "threshold"
#cairo_pdf("aic_errors.pdf", height=6, width=6)
par(mfrow=c(2,2))
plot_error(bm_v_ouLP, "BM", "OU.3", info)
plot_error(bm_v_ou1, "BM", "OU.1", info)
plot_error(ouLP_v_ouLP4, "OU.3", "OU.4", info)
plot_error(ouLP_v_half, "OU.3", "OU.15", info)
dev.off()





