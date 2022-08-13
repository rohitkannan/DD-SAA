# Monte Carlo estimation of the out-of-sample cost of a given first-stage decision

# estimate the objective value of a given first-stage solution
function estimateSolnQuality(z_soln::Array{Float64},covariate_obs::Array{Float64},degree::Float64,numMCScenarios::Int64,numMCReplicates::Int64,coeff_true::Array{Float64,2},var_coeff_true::Array{Float64,2},var_coeff_scaling::Array{Float64})
	
	objEstimates = zeros(Float64,numMCReplicates)

	for mc = 1:numMCReplicates
		srand(randomSeeds_MC[mc])
		demand_scen_MC = generateTrueCondScenarios(numMCScenarios,covariate_obs,degree,coeff_true,var_coeff_true,var_coeff_scaling)		
		objEstimates[mc] = estimateCostOfSoln(z_soln,demand_scen_MC)
	end
	
	return objEstimates
end