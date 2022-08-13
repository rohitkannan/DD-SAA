# Generate scenarios of the random variables


# generate samples using the ER-SAA approach
function generateERSAAScenarios(outputs::Array{Float64,2},inputs::Array{Float64,2},covariates_obs::Array{Float64,1},regressionMethod::String,estimate_het::Bool)

	numDependentVar::Int64 = size(outputs,2)
	numIndependentVar::Int64 = size(inputs,2)
	numSamples::Int64 = size(outputs,1)
	demand_scen_ERSAA = zeros(Float64,numSamples,numDependentVar)
	err_tol_het::Float64 = 1E-02
		
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

		if estimate_het
			err_outputs = log.(max.(err_tol_het^2, (outputs - inputs*coeff_full).^2))
			err_inputs = log.(abs.(inputs .+ 1))
	        err_coeff = getLassoCoeff(err_outputs,err_inputs)

			add_err = zeros(Float64,numSamples,numDependentVar)
			for d = 1:numDependentVar
				add_err[:,d] = sqrt(exp(sum(err_coeff[i,d]*log(1 + covariates_obs[i]) for i = 1:numIndependentVar))) .* ((outputs[:,d] - inputs*coeff_full[:,d]) ./ sqrt.(exp.(sum(err_coeff[i,d]*log.(1 .+ inputs[:,i]) for i = 1:numIndependentVar))))
			end

			coeff_full_updated = zeros(Float64,numIndependentVar,numDependentVar)

			for d = 1:numDependentVar
				scaling_mat_inputs = zeros(Float64,numSamples,numIndependentVar)
        	    scaling_mat_outputs = zeros(Float64,numSamples,1)
				for samp = 1:numSamples
					scaling_mat_inputs[samp,:] = inputs[samp,:] ./ sqrt(exp(sum(err_coeff[i,d]*log(1 + inputs[samp,i]) for i = 1:numIndependentVar)))
                    scaling_mat_outputs[samp,1] = outputs[samp,d] / sqrt(exp(sum(err_coeff[i,d]*log(1 + inputs[samp,i]) for i = 1:numIndependentVar)))
				end

				XtX_fact_updated = lufact((inputs')*scaling_mat_inputs)
				XtY_updated = (inputs')*scaling_mat_outputs

	            coeff_full_updated[:,d] = XtX_fact_updated\(XtY_updated)
	        end

			for d = 1:numDependentVar
				add_err[:,d] = sqrt(exp(sum(err_coeff[i,d]*log(1 + covariates_obs[i]) for i = 1:numIndependentVar))) .* ((outputs[:,d] - inputs*coeff_full_updated[:,d]) ./ sqrt.(exp.(sum(err_coeff[i,d]*log.(1 .+ inputs[:,i]) for i = 1:numIndependentVar))))
			end
			
			demand_scen_ERSAA = repeat(((covariates_obs')*coeff_full_updated)', outer=[1,numSamples])' + add_err
		else
			demand_scen_ERSAA = repeat(((covariates_obs')*coeff_full)', outer=[1,numSamples])' + outputs - inputs*coeff_full
		end

	elseif(regressionMethod == "lasso")
		
		coeff_full = getLassoCoeff(outputs,inputs)

		if estimate_het
			throw(ErrorException("Heteroscedasticity estimation not implemented for Lasso regression"))
		else
			demand_scen_ERSAA = repeat(((covariates_obs')*coeff_full)', outer=[1,numSamples])' + outputs - inputs*coeff_full
		end

	elseif(regressionMethod == "kNN")
				
		soln_kNN = generatekNNScenarios(outputs,inputs,numFolds,minExponent,maxExponent,covariates_obs,"ER-SAA")
		if estimate_het
			err_outputs = log.(max.(err_tol_het^2, soln_kNN["additive_error"].^2))
			err_inputs = log.(abs.(inputs .+ 1))
            err_coeff = getLassoCoeff(err_outputs,err_inputs)

			add_err = zeros(Float64,numSamples,numDependentVar)
			for d = 1:numDependentVar
				add_err[:,d] = sqrt(exp(sum(err_coeff[i,d]*log(1 + covariates_obs[i]) for i = 1:numIndependentVar))) .* (soln_kNN["additive_error"][:,d] ./ sqrt.(exp.(sum(err_coeff[i,d]*log.(1 .+ inputs[:,i]) for i = 1:numIndependentVar))))
			end
			
			demand_scen_ERSAA = repeat(soln_kNN["demand_point_pred"], outer=[1,numSamples])' + add_err
		else
			demand_scen_ERSAA = soln_kNN["demand_scen"]
		end
	end
	
	return demand_scen_ERSAA
end
