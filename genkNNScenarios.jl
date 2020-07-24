# Generate scenarios of the demands for use within the kNN-SAA framework of Bertsimas and Kallus


# find the indices of the k minimum elements of an array
function mink(arr::Array{Float64}, k_val::Int64)

	v = (sortperm(arr,rev=false))[1:k_val]
	return v
end


# compute distances between covariate samples for kNN
function computeSampCovDist(inputs::Array{Float64,2})

	numSamples::Int64 = size(inputs,1)
	dist_cov = zeros(Float64,numSamples,numSamples)
	for s = 1:numSamples
		for t = s+1:numSamples
			dist_cov[s,t] = norm(inputs[s,:]-inputs[t,:],2)
			dist_cov[t,s] = dist_cov[s,t]
		end
	end

	return dist_cov
end


function generatekNNScenarios(outputs::Array{Float64,2},inputs::Array{Float64,2},numFolds::Int64,minExponent::Float64,maxExponent::Float64,covariate_obs::Array{Float64,1})	

	dist_cov = computeSampCovDist(inputs)


	numDependentVar::Int64 = size(outputs,2)
	numSamples::Int64 = size(outputs,1)
	

	min_error_est::Float64 = Inf
	numRepeatCV::Int64 = 1
	k_sweep = floor.(Int64,numSamples^minExponent):ceil.(Int64,numSamples^maxExponent)
	k_optimal::Int64 = k_sweep[1]
	error_k = zeros(Float64,length(k_sweep))


	# replicates of k-fold CV, if needed
	for rep = 1:numRepeatCV

		# first, randomly permute the given data points for k-fold CV
		rand_perm = randperm(numSamples)

		output_perm = zeros(Float64,numSamples,numDependentVar)
		for d = 1:numDependentVar
			output_perm[:,d] = outputs[rand_perm,d]
		end

		
		# next, compute distances for kNN
		dist_perm_tmp = zeros(Float64,numSamples,numSamples)
		dist_perm = zeros(Float64,numSamples,numSamples)
		for s = 1:numSamples
			dist_perm_tmp[:,s] = dist_cov[:,rand_perm[s]]
		end
		for s = 1:numSamples
			dist_perm[s,:] = dist_perm_tmp[rand_perm[s],:]
		end


		# then, compute the folds of data for k-fold CV
		numBaseElements::Int64 = floor.(Int64,numSamples*1.0/numFolds)
		numRemElements::Int64 = numSamples - numFolds*numBaseElements
		index_folds = zeros(Int64,numFolds,2)
		index_folds[1,1] = 1
		index_folds[1,2] = numBaseElements
		if(numRemElements > 0)
			index_folds[1,2] += 1
			numRemElements -= 1
		end
		for f = 2:numFolds
			index_folds[f,1] = index_folds[f-1,2] + 1
			index_folds[f,2] = index_folds[f,1] + (numBaseElements-1)
			if(numRemElements > 0)
				index_folds[f,2] += 1
				numRemElements -= 1
			end
		end
		

		# compute the k nearest neighbors
		ind_dist_sort = zeros(Int64,numSamples,k_sweep[end])
		for f = 1:numFolds
			for s = index_folds[f,1]:index_folds[f,2]
				dist_s_f = dist_perm[s,:]
				dist_s_f[index_folds[f,1]:index_folds[f,2]] = Inf
				ind_dist_sort[s,:] = mink(dist_s_f,k_sweep[end])
			end
		end
	
	
		# try different values of the parameter k and pick one that yields the smallest CV prediction error
		for i = 1:length(k_sweep)
			k_val::Int64 = k_sweep[i]
			for s = 1:numSamples
				est_output_s = mean(output_perm[ind_dist_sort[s,1:k_val],:],1)'
				error_k[i] += norm(est_output_s-output_perm[s,:],2)^2
			end
		end
	
	end

	
	for i = 1:length(k_sweep)
		if(error_k[i] < min_error_est)
			min_error_est = error_k[i]
			k_optimal = k_sweep[i]
		end
	end


	# finally, use kNN to predict the scenarios for the demands given a new realization of the covariates
	dist_obs = zeros(Float64,numSamples)
	for s = 1:numSamples
		dist_obs[s] = norm(inputs[s,:]-covariate_obs,2)
	end
	min_ind = mink(dist_obs,k_optimal)

	numSAAScenarios::Int64 = size(min_ind,1)
	demand_scen_kNN = zeros(Float64,numSAAScenarios,numDependentVar)
	for s = 1:numSAAScenarios
		demand_scen_kNN[s,:] = outputs[min_ind[s],:]
	end
	
	
	return demand_scen_kNN
end
