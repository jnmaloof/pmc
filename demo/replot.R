
###############
require(socialR)
script <- "replot.R"
gitaddr <- gitcommit(script)
tags="phylogenetics"
tweet_errors(script, tags=tags)
###############

require(pmc)
load("early_burst.Rdat")

upload("eb.png", script="early_burst_example.R")
