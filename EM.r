library("robustbase")


EM.cov= function(data, threshold = 0.00001, maxiter = 1000, verbose = F, use.robust = F)
{
	data = as.matrix(data)
	pattern = is.na(data)
	pattern.idx = integer(nrow(data))

	pattern.db = matrix(pattern[1, ], 1, ncol(pattern))
	pattern.idx[1] = 1
	for(i in 2:nrow(data)) {
		for(j in 1:nrow(pattern.db)) {
			if(all(pattern[i, ] == pattern.db[j, ])) {
				pattern.idx[i] = j
				break
			}
		}

		if(pattern.idx[i] == 0) {
			pattern.db = rbind(pattern.db, pattern[i, ])
			pattern.idx[i] = nrow(pattern.db)
		}
	}

	for(i in 1:ncol(data)) {
		if(!use.robust)
			data[is.na(data[, i]), i] = mean(data[, i], na.rm = T)
		else
			data[is.na(data[, i]), i] = median(data[, i], na.rm = T)
	}

	if(!use.robust) {
		S = cov(data)
		mu = colMeans(data, na.rm = T)
	}
	else {
		temp = rbcov(data)
		S = temp$cov
		mu = temp$center
	}
	ll = numeric(maxiter)

	covs = vector("list", nrow(pattern.db))
	means = vector("list", nrow(pattern.db))
	for(i in 1:nrow(pattern.db)) {
		var.obs = !pattern.db[i, ]
		data.idx = pattern.idx == i
			
		if(!use.robust) {
			covs[[i]] = cov(data[data.idx, var.obs])
			means[[i]] = colMeans(data[data.idx, var.obs])
		}
		
		else {
			temp = rbcov(data[data.idx, var.obs])
			covs[[i]] = temp$cov
			means[[i]] = temp$center
		}
	}	

	for(i in 1:maxiter) {
		next.S = matrix(0, ncol(data), ncol(data))
		
		for(j in 1:nrow(pattern.db)) {
			var.obs = !pattern.db[j, ]
			data.idx = pattern.idx == j
			num.data = sum(data.idx)

			project.mat = S[!var.obs, var.obs] %*% solve(S[var.obs, var.obs])
			obs.mu = num.data * (means[[j]] - mu[var.obs])
			mu[var.obs] = mu[var.obs] + obs.mu / nrow(data)
			mu[!var.obs] = mu[!var.obs] + project.mat %*% obs.mu / nrow(data)
		}

		for(j in 1:nrow(pattern.db)) {
			var.obs = !pattern.db[j, ]
			data.idx = pattern.idx == j
			num.data = sum(data.idx)

			inv.mat = solve(S[var.obs, var.obs])
			project.mat = S[!var.obs, var.obs] %*% inv.mat
			diff.mean = matrix(means[[j]] - mu[var.obs], sum(var.obs), 1)
			obs.S = num.data * (covs[[j]] + diff.mean %*% t(diff.mean))

			next.S[var.obs, var.obs] = next.S[var.obs, var.obs] + obs.S
			next.S[!var.obs, var.obs] = next.S[!var.obs, var.obs] + project.mat %*% obs.S
			next.S[var.obs, !var.obs] = t(next.S[!var.obs, var.obs])
			next.S[!var.obs, !var.obs] = next.S[!var.obs, !var.obs] + num.data * S[!var.obs, !var.obs] + 
				project.mat %*% (obs.S - num.data * S[var.obs, var.obs]) %*% t(project.mat)
			ll[i] = ll[i] + num.data * log(det(inv.mat)) - sum(eigen(inv.mat %*% obs.S)$value)
		}

		next.S = next.S / nrow(data)
		ll[i] = ll[i] / nrow(data)

		if(verbose) 
			cat(sum(abs(next.S - S)), "\t", ll[i], "\n")

		S = next.S
		if(i > 1)
			if(abs((ll[i] - ll[i-1]))< threshold)
				break	
	}

	return(list(S = S, mu = mu, ll = ll[1:i]))
}

rbcov = function(data, mu = NULL, ...) {
	if(is.null(mu))		
		return(covOGK(data, sigmamu = rbsd, ...))
	else {
		data = t(t(data) - mu)
		return(covOGK(data, sigmamu = rbsdcentered, ...))
	}
}


