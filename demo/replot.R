
###############
require(socialR)
script <- "replot.R"
gitaddr <- gitcommit(script)
tags="phylogenetics"
tweet_errors(script, tags=tags)
###############

require(pmc)
load("simtree_lambda_dist.Rdat")

png("simtree_eb.png")
hist(bm_v_eb$test_par_dist[3,], col=rgb(0,0,1,.5), border="white", breaks=15, main="", xlab="Estimated a")
abline(v=eb[[1]]$a, lwd=3, lty=2, col="darkred") #True value
dev.off()

require(socialR)
upload("simtree_eb.png", script="simtree_eb.R")

