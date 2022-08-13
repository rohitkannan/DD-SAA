# Generate scenarios of the demands for use within the kNN-based reweighted SAA framework of Bertsimas and Kallus
function generateBKScenarios(outputs::Array{Float64,2},inputs::Array{Float64,2},numFolds::Int64,minExponent::Float64,maxExponent::Float64,covariate_obs::Array{Float64,1})

	soln_kNN = generatekNNScenarios(outputs,inputs,numFolds,minExponent,maxExponent,covariate_obs,"reweighted-SAA")
	return soln_kNN["demand_scen"]
end