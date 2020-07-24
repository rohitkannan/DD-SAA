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
# generates a common set of coefficients across all replicates
function generateBaseDemandModel(numCovariates::Int64)

	# generate parameters for the covariate distribution
	covariate_mean = zeros(Float64, maxNumCovariates)
	covariate_corrMat = generateRandomCorrMat(maxNumCovariates)
	covariate_covMat = covariate_corrMat + 1E-09*eye(maxNumCovariates)

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
	
	return covariate_mean, covariate_covMat, coeff_true
end



# generate demand and covariate data using a sparse (non)linear model
function generateDemandData(numCovariates::Int64,numSamples::Int64,degree::Float64,covariate_mean::Array{Float64},covariate_covMat::Array{Float64,2},coeff_true::Array{Float64,2})

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
	
	demand_errors = demand_errors_scaling*rand(Normal(), maxNumDataSamples, numCustomers)
	demand_data = ((covariate_data).^(degree))*coeff_true + demand_errors[1:numSamples,:]
		
	return demand_data, covariate_data
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
function generateTrueCondScenarios(numScenarios::Int64,covariate_obs::Array{Float64,1},degree::Float64,coeff_true::Array{Float64,2})

	demand_errors = demand_errors_scaling*rand(Normal(), numScenarios, numCustomers)
	demand_scen = repeat((((covariate_obs).^(degree))'*coeff_true)', outer = [1,numScenarios])' + demand_errors
	
	return demand_scen
end
