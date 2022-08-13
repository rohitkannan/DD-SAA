# Generate scenarios of the random variables

# generate a point prediction of the random variables for the given covariate observation
function generatePointPrediction(outputs::Array{Float64,2},inputs::Array{Float64,2},covariates_obs::Array{Float64,1},regressionMethod::String)

	numDependentVar::Int64 = size(outputs,2)
	numIndependentVar::Int64 = size(inputs,2)
	numSamples::Int64 = size(outputs,1)
	demand_scen_pointpred = zeros(Float64,1,numDependentVar)
		
	if(size(outputs,1) != size(inputs,1))
		throw(ErrorException("Number of output samples != number of input samples"))
	end
		
	if(regressionMethod == "ols")
	
		if(numSamples < numIndependentVar)
			throw(ErrorException("Expected number of samples > number of covariates in OLS!"))
		end


		XtX_fact = lufact((inputs')*inputs)
		XtY = (inputs')*outputs
		
		coeff_full = zeros(Float64,numIndependentVar,numDependentVar)
		for d = 1:numDependentVar
			coeff_full[:,d] = XtX_fact\(XtY[:,d])
		end
				
		demand_scen_pointpred[1,:] = (covariates_obs')*coeff_full

	elseif(regressionMethod == "lasso")
		
		coeff_full = getLassoCoeff(outputs,inputs)
				
		demand_scen_pointpred[1,:] = (covariates_obs')*coeff_full

	end
	
	return demand_scen_pointpred
end
