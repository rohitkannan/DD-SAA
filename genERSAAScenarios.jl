# Generate scenarios of the demands for use within the ER-SAA framework


function generateERSAAScenarios(outputs::Array{Float64,2},inputs::Array{Float64,2},covariates_obs::Array{Float64,1},regressionMethod::String,coeff_true::Array{Float64,2})

	numDependentVar::Int64 = size(outputs,2)
	numIndependentVar::Int64 = size(inputs,2)
	numSamples::Int64 = size(outputs,1)
	demand_scen_ERSAA = zeros(Float64,numSamples,numDependentVar)
	
	fitQuality::Float64 = 1.0
		
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
				
		demand_scen_ERSAA = repeat(((covariates_obs')*coeff_full)', outer=[1,numSamples])' + outputs - inputs*coeff_full
		
		fitQuality = norm(coeff_full-coeff_true,2)

	elseif(regressionMethod == "lasso")
		
		coeff_full = getLassoCoeff(outputs,inputs)
				
		demand_scen_ERSAA = repeat(((covariates_obs')*coeff_full)', outer=[1,numSamples])' + outputs - inputs*coeff_full
		
		fitQuality = norm(coeff_full-coeff_true,2)

	end
	
	return demand_scen_ERSAA, fitQuality
end
