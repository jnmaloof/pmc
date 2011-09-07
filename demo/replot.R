
###############
require(socialR)
script <- "replot.R"
gitaddr <- gitcommit(script)
tags="phylogenetics"
tweet_errors(script, tags=tags)
###############

require(pmc)
load("early_burst.Rdat")
png("bm_v_eb.png")
plot(bm_v_eb)
dev.off()

png("ou_v_eb.png")
plot(ou_v_eb)
dev.off()

upload("bm_v_eb.png", script=script, gitaddr=gitaddr, tags=tags)
upload("ou_v_eb.png", script=script, gitaddr=gitaddr, tags=tags)


