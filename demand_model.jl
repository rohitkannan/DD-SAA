# Generate random demand and covariate data for the resource allocation model


# Function to generate random correlation matrices. Based on https://stats.stackexchange.com/questions/124538/how-to-generate-a-large-full-rank-random-correlation-matrix-with-some-strong-cor
function generateRandomCorrMat(dim::Int64)

	betaparam::Float64 = 2.0

	partCorr = zeros(Float64,dim,dim)
	corrMat = eye(dim)

	for k = 1:dim-1
		for i = k+1:dim
			partCorr[k,i] = ((rand(Beta(betaparam,betaparam),1))[1] - 0.5)*2.0
			p::Float64 = partCorr[k,i]
			for j = (k-1):-1:1
				p = p*sqrt((1-partCorr[j,i]^2)*(1-partCorr[j,k]^2)) + partCorr[j,i]*partCorr[j,k]
			end
			corrMat[k,i] = p
			corrMat[i,k] = p
		end
	end

	permut = randperm(dim)
	corrMat = corrMat[permut, permut]

	return corrMat
end


# generate base parameters for the demand model given number of covariates
function generateBaseDemandModel(numCovariates::Int64)

	# generate parameters for the covariate distribution
	covariate_mean = zeros(Float64, maxNumCovariates)
	covariate_corrMat = generateRandomCorrMat(maxNumCovariates)
	covariate_covMat = covariate_corrMat + 1E-09*eye(maxNumCovariates)
	# covariate realizations for estimating scaling coefficients
	covariate_tmp = rand(Normal(),maxNumCovariates,maxNumDataSamples)

	# generate coefficients for the dependence of the demands on the covariates
	alpha = 50.0 + 5.0*rand(Normal(), numCustomers)
	beta = repeat([10.0,5.0,2.0],outer=[1,numCustomers])' + rand(Uniform(-4.0,4.0), numCustomers, 3)
	
	if(numCovariates < 3)
		throw(ErrorException("Expected numCovariates >= 3!"))
	end

	coeff_true = zeros(Float64,numCovariates+1,numCustomers)
	coeff_true[1,:] = alpha
	coeff_true[2,:] = beta[:,1]
	coeff_true[3,:] = beta[:,2]
	coeff_true[4,:] = beta[:,3]


	var_coeff_true = zeros(Float64,numCovariates+1,numCustomers)
	var_scaling = 2*(het_level - 1)^2
	var_coeff_true[2:4,:] = var_scaling*rand(Uniform(0.0,1.0), 3, numCustomers)


	# estimate scaling coefficients
	covariate_mean_case = covariate_mean[1:numCovariates]
	covariate_covMat_case = covariate_covMat[1:numCovariates,1:numCovariates]
	covariate_covMat_chol = (cholfact(Hermitian(covariate_covMat_case)))[:L]
	
	covariate_data = zeros(Float64, maxNumDataSamples, numCovariates+1)
	covariate_data[:,1] = ones(Float64,maxNumDataSamples)
	for s = 1:maxNumDataSamples
		covariate_data[s,2:numCovariates+1] = covariate_covMat_chol*covariate_tmp[1:numCovariates,s] + covariate_mean_case
	end
	covariate_data = abs.(covariate_data)

	var_coeff_scaling = ones(Float64,numCustomers)
	for j = 1:numCustomers
		var_coeff_scaling[j] = median(exp.(sum(var_coeff_true[i,j] .* log.(1 .+ covariate_data[:,i]) for i = 2:4)))
	end

	
	return covariate_mean, covariate_covMat, coeff_true, var_coeff_true, var_coeff_scaling
end



# generate demand and covariate data using a sparse (non)linear model
function generateDemandData(numCovariates::Int64,numSamples::Int64,degree::Float64,covariate_mean::Array{Float64},covariate_covMat::Array{Float64,2},coeff_true::Array{Float64,2},var_coeff_true::Array{Float64,2},var_coeff_scaling::Array{Float64})

	covariate_mean_case = covariate_mean[1:numCovariates]
	covariate_covMat_case = covariate_covMat[1:numCovariates,1:numCovariates]
	covariate_covMat_chol = (cholfact(Hermitian(covariate_covMat_case)))[:L]

	#*===========================================
	# first, construct the the covariates	
	# first covariate is simply the constant one (for the intercept)
	# the next set of covariates are iid samples from a multivariate folded normal distribution
	
	covariate_data = zeros(Float64, numSamples, numCovariates+1)

	covariate_data[:,1] = ones(Float64,numSamples)
	covariate_tmp = rand(Normal(),maxNumCovariates,maxNumDataSamples)
	for s = 1:numSamples
		covariate_data[s,2:numCovariates+1] = covariate_covMat_chol*covariate_tmp[1:numCovariates,s] + covariate_mean_case
	end
	covariate_data = abs.(covariate_data)
	
	#*===========================================
	# next, construct the random demands from the covariates using a sparse (non)linear model
	
	demand_errors_indep = demand_errors_scaling*rand(Normal(), maxNumDataSamples, numCustomers)
	demand_errors_var = ones(Float64, numSamples, numCustomers)
	for j = 1:numCustomers
		demand_errors_var[:,j] = exp.(sum(var_coeff_true[i,j] .* log.(1 .+ covariate_data[:,i]) for i = 1:numCovariates+1)) ./ var_coeff_scaling[j]
	end
	demand_data = ((covariate_data).^(degree))*coeff_true + (sqrt.(demand_errors_var) .* demand_errors_indep[1:numSamples,:])
		
	return demand_data, covariate_data, demand_errors_var
end



# generate a new realization of the covariates from a multivariate folded normal distribution
function generateCovariateReal(numCovariates::Int64,covariate_mean::Array{Float64},covariate_covMat::Array{Float64,2})

	covariate_mean_case = covariate_mean[1:numCovariates]
	covariate_covMat_case = covariate_covMat[1:numCovariates,1:numCovariates]
	covariate_covMat_chol = (cholfact(Hermitian(covariate_covMat_case)))[:L]

	covariate_tmp = rand(Normal(),maxNumCovariates)
	covariate_obs = ones(Float64, numCovariates+1)
	covariate_obs[2:numCovariates+1] = covariate_covMat_chol*(covariate_tmp[1:numCovariates]) + covariate_mean_case
	covariate_obs = abs.(covariate_obs)
	
	return covariate_obs
end



# generate samples from the true conditional distribution of the demands given a covariate realization
function generateTrueCondScenarios(numScenarios::Int64,covariate_obs::Array{Float64,1},degree::Float64,coeff_true::Array{Float64,2},var_coeff_true::Array{Float64,2},var_coeff_scaling::Array{Float64})

	demand_errors_indep = demand_errors_scaling*rand(Normal(), numScenarios, numCustomers)
	demand_errors_var = ones(Float64, numScenarios, numCustomers)
	for j = 1:numCustomers
		demand_errors_var[:,j] .= exp(sum(var_coeff_true[i,j] * log(1 + covariate_obs[i]) for i = 1:numCovariates+1)) / var_coeff_scaling[j]
	end
	demand_scen = repeat((((covariate_obs).^(degree))'*coeff_true)', outer = [1,numScenarios])' + (sqrt.(demand_errors_var) .* demand_errors_indep)
	
	return demand_scen
end
