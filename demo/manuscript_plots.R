require(pmc)


##################### FIGURE 1 ########################################
pdf("boettiger_figure1a.pdf")
  data(geospiza_lambda) # option to use saved data from the manuscript
  hist(bm_v_lambda$test_par_dist[3,], col=rgb(0,0,1,.5),
  border="white", breaks=15, main="", xlab="Estimated lambda")
  abline(v=lambda[[1]][3], lwd=3, lty=2, col="darkred") #True value
  text(lambda[[1]][3], 300, "True lambda", pos=2)
  dev.off()

pdf("boettiger_figure1b.pdf")
  data(simtree_lambda_dist)
  hist(bm_v_lambda$test_par_dist[3,], col=rgb(0,0,1,.5),
  border="white", breaks=15, main="", xlab="Estimated lambda")
  abline(v=lambda[[1]][3], lwd=3, lty=2, col="darkred") #True value
  text(lambda[[1]][3], 300, "True lambda", pos=2)
dev.off()

pdf("boettiger_figure1c.pdf")
  data(simtree_lambda_dist)
  hist(bm_v_lambda$test_par_dist[2,], col=rgb(0,0,1,.5),
  border="white", breaks=15, main="", xlab="Estimated sigma")
  abline(v=lambda[[1]][2], lwd=3, lty=2, col="darkred") #True value
  text(lambda[[1]][2], 300, "True sigma", pos=4)
dev.off()



##################### FIGURE 3 ########################################
pdf("boettiger_figure3.pdf")
  data(anoles_model_choice)
  layout(cbind(1,2,3,4), widths=c(1,1,1,1))
  par(mar=c(5,.1,.1,.1))
  plot(ouLP_tree, edge.color = treepalette(ouLP_tree), edge.width=5, show.tip.label=FALSE)
  mtext("(a) OU.3", 1, cex=1.2)
  plot(ouLP4_tree, edge.color = treepalette(ouLP4_tree), edge.width=5, show.tip.label=FALSE)
  mtext("(b) OU.4", 1, cex=1.2)
  plot(half_tree, edge.color = treepalette(half_tree), edge.width=5, show.tip.label=FALSE)
  mtext("(c) OU.15", 1, cex=1.2)
  plot(half_tree, edge.color ="white", edge.width=.1, cex=1., show.tip.label=TRUE, label.offset=-1.2)
dev.off()


##################### FIGURE 4 ########################################
cairo_pdf("boettiger_figure4.pdf")
  data(anoles_model_choice)
  par(mfrow=c(2,2))
  plot(bm_v_ouLP, show_text=FALSE, main="(a) BM vs. OU.3")
  legend("topright", c("BM sims", "OU.3 sims", "obs"), pch=c(15,15,46),
  lty=c(0,0,2), col=c(rgb(0,0,1,.5), rgb(1,0,0,.5), "darkred"), bty="n")

  plot(ouLP_v_half, show_text=FALSE, main="(b) OU.3 vs. OU.15")
  legend("topright", c("OU.3 sims", "OU.15 sims", "obs"), pch=c(15,15,46),
  lty=c(0,0,2), col=c(rgb(0,0,1,.5), rgb(1,0,0,.5), "darkred"), bty="n")

  plot(ouLP_v_ouLP4, show_text=FALSE, main="(c) OU.3 vs. OU.4")
  legend("topright", c("OU.3 sims", "OU.4 sims", "obs"), pch=c(15,15,46),
  lty=c(0,0,2), col=c(rgb(0,0,1,.5), rgb(1,0,0,.5), "darkred"), bty="n")

  plot(bm_v_ou1, show_text=FALSE, main="(d) BM vs. OU.1")
  legend("topright", c("BM sims", "OU.1 sims", "obs"), pch=c(15,15,46),
  lty=c(0,0,2), col=c(rgb(0,0,1,.5), rgb(1,0,0,.5), "darkred"), bty="n")
dev.off()



##################### FIGURE 5 ########################################
cairo_pdf("boettiger_figure5.pdf")
  data(anoles_model_choice)
  plot_error <- function(pow, main, info){
      plot(pow, shade_aic=T, show_aic=T, shade=F, show_data=F, show_text=F, info=info, legend=TRUE, main=main)
  }
  info <- "aic"
  par(mfrow=c(2,2))
  plot_error(bm_v_ouLP, "(a) BM vs. OU.3", info)
  plot_error(bm_v_ou1, "(b) BM vs. OU.1", info)
  plot_error(ouLP_v_ouLP4, "(c) OU.3 vs. OU.4", info)
  plot_error(ouLP_v_half, "(d) OU.3 vs. OU.15", info)
dev.off()



##################### FIGURE 6 ########################################
cairo_pdf("boettiger_figure6a.pdf")

data("power_curves")
# n is a vec of number of taxa used in each sim: 
k <- length(n)-2
plot(1,1, type='n', xlim=c(min(alpha), max(alpha)), ylim = c(0,1),
   main="Power by tree size", log="x", xlab="alpha", ylab="power")
for(i in 1:k ){ ## skip the last 2, which haven't converged
  points(alpha, size[[i]]$power, pch=16, col=i)
  lines(alpha, size[[i]]$power, col=i)
}
points(alpha, anoles$power, pch=16, col="purple")
lines(alpha, anoles$power, col="purple", lwd=4)
legend("topleft", c(paste(n[1:k], "taxa"), "23 (anoles)"), col=c(1:k,"purple"), pch=16  ) 



dev.off()


data("power_curves")
plot_shape <- function(){
  plot(1,1, type='n', xlim=c(min(alpha), max(alpha)), ylim = c(0,1),
  main="Power by tree topology", log="x", xlab="alpha", ylab="power")
    k <- length(lambda)
  for(i in 1:length(n)){
    points(alpha, shape[[i]]$power, pch=16, col=i)
    lines(alpha, shape[[i]]$power,  col=i)
  }
  points(alpha, anoles$power, pch=16, col="purple")
  lines(alpha, anoles$power, col="purple", lwd=4)
  legend("topleft", c(paste(lambda, "lambda"), "anoles"), col=c(1:k, "purple"), pch=16  ) 
}


cairo_pdf("boettiger_figure6b.pdf")
  plot_shape()
dev.off()
