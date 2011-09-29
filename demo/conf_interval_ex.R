# Compares confidence interval estimates from the Hessian approach, 
# (removed from more recent versions of fitContinuous) with the bootstrap
# apporoach.  Shows that the Hessian approximation can perform very poorly
# Won't run without custom compile of geiger with hessian option added back in  

require(pmc)
data(geospiza_lambda)
# fit the model with the hessian returned
lambda <- fitContinuous(data$phy, data$data, model="lambda", hessian=TRUE)
out <- lambda[[1]]
# Estimate confidence intervals from the hessian: can be highly unreliable 
ci <- sqrt(diag(solve(out$hessian)))
print(ci)

# Confidence intervals using the PMC package.  Compare to the above.
o <- confidenceIntervals.pow(bm_v_lambda)
print(paste("True CI, lambda:", o[["test"]][,"lambda"]))
print(paste("Trie CI, beta:", o[["test"]][,"beta"]))



# Compute confidence intervals of lambda and sigma for the large tree
data(simtree_lambda_dist)
lambda <- fitContinuous(data$phy, data$data, model="lambda", hessian=TRUE)
out <- lambda[[1]]
## Hessian method
ci <- sqrt(diag(solve(out$hessian)))
print(paste("Hessian-based approximation:", ci))

# PMC method
o <- confidenceIntervals.pow(bm_v_lambda)
print(paste("True CI, lambda:", o[["test"]][,"lambda"]))
print(paste("True CI, beta:", o[["test"]][,"beta"]))


