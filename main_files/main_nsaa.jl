# Main file for solving the naive SAA problem

using Distributions, StatsFuns, StatsBase, JuMP, Gurobi


const gurobi_env = Gurobi.Env()


# read these as inputs from the corresponding "instance file"
modRepNum = parse(ARGS[1])
degNum = parse(ARGS[2])
dataRepNum = parse(ARGS[3])
sampSizeNum = parse(ARGS[4])



#*========= MODEL INFORMATION ==========
const caseNum = 1
const paramsFile = "params_case" * string(caseNum) * "_nsaa.jl"
#*======================================



#*========= INCLUDE FILES ====================
include("maxParameters.jl")
include("randomSeeds.jl")

include(paramsFile)
# the next two data files are defined in paramsFile
include(modelFile)
include(dataFile)

include("solveModel.jl")
include("estimateSolnCost.jl")
include("evaluateSoln.jl")
#*======================================



# set seed for reproducibility
srand(startingSeed)



# directory name for storing results
const baseDirName = "case" * string(caseNum) * "_nsaa/" * "mod_" * string(modRepNum) * "/" * "deg_" * string(degNum) * "/"
const subDirName = baseDirName * "rep_" * string(dataRepNum) * "/"
const subDirName2 = subDirName * "samp_" * string(sampSizeNum) * "/"
mkpath(subDirName2)



#*========= PRINT OPTIONS ====================
const storeResults = true

const infoFile = "ddsp.txt"
const modelDataFile = "model_data.txt"
const numSampFile = "num_samples.txt"

const nsaaObjFile = "nsaa_obj.txt"
const nsaaDDObjFile = "nsaa_ddobj.txt"
const nsaaTimeFile = "nsaa_time.txt"
#*======================================
	


#*========= GENERATE MODEL PARAMETERS ====================
	
srand(randomSeeds_models[modRepNum])	
numCovariates = 3

covariate_mean, covariate_covMat, coeff_true, var_coeff_true, var_coeff_scaling = generateBaseDemandModel(numCovariates)


#*========= STORE RESULTS ====================
if(storeResults)
	
	# write details to text file, including some key details about the test instance
	details_file = baseDirName * infoFile
	open(details_file, "w") do f
		write(f,"case number: $caseNum \n")
		write(f,"numResources: $numResources \n")
		write(f,"numCustomers: $numCustomers \n")
		write(f,"model replicate number: $modRepNum \n")
		write(f,"degree: $(degree[degNum]) \n")
		write(f,"sample sizes: $numDataSamples \n")
		write(f,"heteroscedasticity level: $het_level \n")
		write(f,"demand_errors_scaling: $demand_errors_scaling \n\n")
		
		write(f,"randomSeeds_MC: $randomSeeds_MC \n")
		write(f,"randomSeeds_models: $randomSeeds_models \n")
		write(f,"randomSeeds_data: $randomSeeds_data \n")
	end	
	
	mod_data_file = baseDirName * modelDataFile
	open(mod_data_file, "w") do f
		write(f,"covariate_mean = $covariate_mean \n")
		write(f,"covariate_covMat = $covariate_covMat \n")
		write(f,"coeff_true = $coeff_true \n")
		write(f,"var_coeff_true = $var_coeff_true \n")
		write(f,"var_coeff_scaling = $var_coeff_scaling \n")
	end
	
	samp_size_file = baseDirName * numSampFile
	open(samp_size_file, "w") do f
		for i = 1:length(numDataSamples)
			write(f,"$(numDataSamples[i]) \n")
		end
	end
	
end
#*============================================



#*========= GENERATE DATA REPLICATE ====================

srand(randomSeeds_data[dataRepNum])

covariate_obs = generateCovariateReal(numCovariates,covariate_mean,covariate_covMat)

demand_data, ~, ~ = generateDemandData(numCovariates,numDataSamples[sampSizeNum],degree[degNum],covariate_mean,covariate_covMat,coeff_true,var_coeff_true,var_coeff_scaling)



#*========= SOLVE N-SAA MODEL ====================

tic()


z_soln_NSAA, objDDNSAA = solveSAAModel(demand_data)

nsaaObjEstimates = estimateSolnQuality(z_soln_NSAA,covariate_obs,degree[degNum],numMCScenarios,numMCReplicates,coeff_true,var_coeff_true,var_coeff_scaling)


NSAATime = toq()


#*========= STORE RESULTS ====================
if(storeResults)
	
	# write data to text file
	obj_est_file = subDirName2 * nsaaObjFile
	open(obj_est_file, "w") do f
		for i = 1:length(nsaaObjEstimates)
			write(f,"$(nsaaObjEstimates[i]) \n")
		end
	end	
	
	ddobj_est_file = subDirName2 * nsaaDDObjFile
	open(ddobj_est_file, "w") do f
		write(f,"$objDDNSAA")
	end	
	
	time_obj_file = subDirName2 * nsaaTimeFile
	open(time_obj_file, "w") do f
		write(f,"$NSAATime")
	end	
	
end
#*============================================
