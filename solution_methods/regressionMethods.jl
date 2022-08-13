# Regression routines

# Routine for estimating Lasso coefficients using GLMNet
function getLassoCoeff(outputs::Array{Float64,2},inputs::Array{Float64,2})

	numDependentVar::Int64 = size(outputs,2)
	numIndependentVar::Int64 = size(inputs,2)
	
	if(size(outputs,1) != size(inputs,1))
		throw(ErrorException("Number of output samples != number of input samples"))
	end

	coeff = zeros(Float64,numIndependentVar,numDependentVar)
	
	for d = 1:numDependentVar
		path_cv = glmnetcv(inputs[:,2:numIndependentVar],outputs[:,d]);
		path_ind::Int64 = indmin(path_cv.meanloss)
		coeff[2:numIndependentVar,d] = path_cv.path.betas[:,path_ind]
	end

	tmpMat = outputs - inputs[:,2:numIndependentVar]*coeff[2:numIndependentVar,:]
	for d = 1:numDependentVar
		coeff[1,d] = mean(tmpMat[:,d])
	end
	
	return coeff
end



# Routines for kNN regression

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


# generate kNN regression scenarios depending on the SAA method
function generatekNNScenarios(outputs::Array{Float64,2},inputs::Array{Float64,2},numFolds::Int64,minExponent::Float64,maxExponent::Float64,covariate_obs::Array{Float64,1},solutionMethod::String)	

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

	soln_kNN = Dict()

	if solutionMethod == "reweighted-SAA"
		numSAAScenarios::Int64 = size(min_ind,1)
		demand_scen_kNN = zeros(Float64,numSAAScenarios,numDependentVar)
		for s = 1:numSAAScenarios
			demand_scen_kNN[s,:] = outputs[min_ind[s],:]
		end
		soln_kNN["demand_scen"] = demand_scen_kNN
	elseif solutionMethod == "ER-SAA"
		point_pred = zeros(Float64,numDependentVar)
		point_pred[:] = mean(outputs[min_ind,:],1)'
	
		demand_scen_kNN = zeros(Float64,numSamples,numDependentVar)
		additive_error_kNN = zeros(Float64,numSamples,numDependentVar)
		for s = 1:numSamples
			ind_dist_sort = mink(dist_cov[s,:],k_optimal)
			est_output = mean(outputs[ind_dist_sort,:],1)'
			additive_error_kNN[s,:] = (outputs[s,:] - est_output)
			demand_scen_kNN[s,:] = (outputs[s,:] - est_output) + point_pred
		end
		soln_kNN["demand_scen"] = demand_scen_kNN
		soln_kNN["demand_point_pred"] = point_pred
		soln_kNN["additive_error"] = additive_error_kNN
	elseif solutionMethod == "PointPrediction"
		demand_scen_kNN = zeros(Float64,1,numDependentVar)
		demand_scen_kNN[1,:] = mean(outputs[min_ind,:],1)'
		soln_kNN["demand_scen"] = demand_scen_kNN
	end
	
	
	return soln_kNN
end
