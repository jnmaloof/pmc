rm(list=ls())
###############
require(socialR)
script <- "early_burst_example.R"
gitaddr <- gitcommit(script)
tags="phylogenetics"
tweet_errors(script, tags=tags)
###############

#early_burst_example.R
require(pmc)
data(geospiza)
# get the data as a single named numeric
traits <- geospiza$geospiza.data[[1]]
names(traits) <- rownames(geospiza$geospiza.data)
## get rid of tips we don't have data for
data <- treedata(geospiza$geospiza.tree,traits)
names(data$data) <- rownames(data$data)

# Okay, fit the models
#bm <-  fitContinuous_object(data$phy, data$data)
eb <- fitContinuous_object(data$phy, data$data, model="EB")
ou <- fitContinuous_object(data$phy, data$data, model="OU")



#bm_v_eb <- montecarlotest(bm, eb, nboot = 200, cpu=16)
ou_v_eb <- montecarlotest(ou, eb, nboot = 200, cpu=16)
png("eb.png")
plot(ou_v_eb)
dev.off()

save(list=ls(), file="early_burst.Rdat")

upload("eb.png", script="early_burst_geospiza.R", tags="phylogenetics", comment="geospiza trait1, ou_v_eb")
