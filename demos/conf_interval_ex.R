require(pmc)
data(simtree_lambda_dist)
lambda <- fitContinuous(data$phy, data$data, model="lambda", hessian=TRUE)
out <- lambda[[1]]
ci <- sqrt(diag(solve(out$hessian)))

ci

o <- confidenceIntervals.pow(bm_v_lambda)
o[["test"]][,"lambda"]
o[["test"]][,"beta"] 

