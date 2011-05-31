require(pmc)

data(geospiza_lambda)
lambda <- fitContinuous(data$phy, data$data, model="lambda", hessian=TRUE)
out <- lambda[[1]]
ci <- sqrt(diag(solve(out$hessian)))
o <- confidenceIntervals.pow(bm_v_lambda)
print(ci)
print(paste("boot lambda:", o[["test"]][,"lambda"]))
print(paste("beta", o[["test"]][,"beta"]))




data(simtree_lambda_dist)
lambda <- fitContinuous(data$phy, data$data, model="lambda", hessian=TRUE)
out <- lambda[[1]]
ci <- sqrt(diag(solve(out$hessian)))
o <- confidenceIntervals.pow(bm_v_lambda)
print(ci)
print(paste("boot lambda:", o[["test"]][,"lambda"]))
print(paste("beta", o[["test"]][,"beta"]))


