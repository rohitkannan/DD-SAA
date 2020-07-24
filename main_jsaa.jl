# Main file for solving the J-SAA and J+-SAA problems


using Distributions, StatsFuns, StatsBase, JuMP, Gurobi, GLMNet


const gurobi_env = Gurobi.Env()



# read these as inputs from the corresponding "instance file"
modRepNum = parse(ARGS[1])
degNum = parse(ARGS[2])
dataRepNum = parse(ARGS[3])
sampSizeNum = parse(ARGS[4])



#*========= MODEL INFORMATION ==========
const caseNum = 2
const paramsFile = "params_case" * string(caseNum) * "_jsaa.jl"
#*======================================



#*========= INCLUDE FILES ====================
include("maxParameters.jl")
include("randomSeeds.jl")

include(paramsFile)
# the next two data files are defined in paramsFile
include(modelFile)
include(dataFile)

include("regressionMethods.jl")
include("genJSAAScenarios.jl")
include("solveModel.jl")
include("estimateSolnCost.jl")
include("evaluateSoln.jl")
#*======================================



# set seed for reproducibility
srand(startingSeed)



# directory name for storing results
const baseDirName = "case" * string(caseNum) * "_jsaa/" * "mod_" * string(modRepNum) * "/" * "deg_" * string(degNum) * "/"
const subDirName = baseDirName * "rep_" * string(dataRepNum) * "/"
const subDirName2 = subDirName * regressionMethod * "/" * "samp_" * string(sampSizeNum) * "/"
mkpath(subDirName2)



#*========= PRINT OPTIONS ====================
const storeResults = true

const infoFile = "ddsp.txt"
const modelDataFile = "model_data.txt"
const numSampFile = "num_samples_" * regressionMethod * ".txt"

const jsaaObjFile = "jsaa_obj.txt"
const jsaaDDObjFile = "jsaa_ddobj.txt"
const jpsaaObjFile = "jpsaa_obj.txt"
const jpsaaDDObjFile = "jpsaa_ddobj.txt"
const jsaaTimeFile = "jsaa_time.txt"

const covRealFile = "covariate_obs.txt"
#*======================================
	


#*========= GENERATE MODEL PARAMETERS ====================
	
srand(randomSeeds_models[modRepNum])

covariate_mean, covariate_covMat, coeff_true = generateBaseDemandModel(numCovariates)


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
		write(f,"demand_errors_scaling: $demand_errors_scaling \n")
		write(f,"regressionMethod: $regressionMethod \n\n")
		
		write(f,"randomSeeds_MC: $randomSeeds_MC \n")
		write(f,"randomSeeds_models: $randomSeeds_models \n")
		write(f,"randomSeeds_data: $randomSeeds_data \n")
	end	
	
	mod_data_file = baseDirName * modelDataFile
	open(mod_data_file, "w") do f
		write(f,"covariate_mean = $covariate_mean \n")
		write(f,"covariate_covMat = $covariate_covMat \n")
		write(f,"coeff_true = $coeff_true \n")
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

demand_data, covariate_data = generateDemandData(numCovariates,numDataSamples[sampSizeNum],degree[degNum],covariate_mean,covariate_covMat,coeff_true)


#*========= STORE RESULTS ====================
if(storeResults)
	
	cov_real_file = subDirName * covRealFile
	open(cov_real_file, "w") do f
		write(f,"covariate_obs = $covariate_obs")
	end	
	
end
#*============================================



#*========= SOLVE THE J-SAA and J+-SAA MODELS ====================

tic()


demand_scen_jsaa, demand_scen_jpsaa = generateJSAAScenarios(demand_data,covariate_data,covariate_obs,regressionMethod)

z_soln_jsaa, objDDJSAA = solveSAAModel(demand_scen_jsaa)

jsaaObjEstimates = estimateSolnQuality(z_soln_jsaa,covariate_obs,degree[degNum],numMCScenarios,numMCReplicates,coeff_true)


z_soln_jpsaa, objDDJpSAA = solveSAAModel(demand_scen_jpsaa)

jpsaaObjEstimates = estimateSolnQuality(z_soln_jpsaa,covariate_obs,degree[degNum],numMCScenarios,numMCReplicates,coeff_true)


jsaaTime = toq()


#*========= STORE RESULTS ====================
if(storeResults)
	
	# write data to text file
	obj_est_file1 = subDirName2 * jsaaObjFile
	open(obj_est_file1, "w") do f
		for i = 1:length(jsaaObjEstimates)
			write(f,"$(jsaaObjEstimates[i]) \n")
		end
	end	
	
	obj_est_file2 = subDirName2 * jpsaaObjFile
	open(obj_est_file2, "w") do f
		for i = 1:length(jpsaaObjEstimates)
			write(f,"$(jpsaaObjEstimates[i]) \n")
		end
	end	
	
	ddobj_est_file1 = subDirName2 * jsaaDDObjFile
	open(ddobj_est_file1, "w") do f
		write(f,"$objDDJSAA")
	end	
	
	ddobj_est_file2 = subDirName2 * jpsaaDDObjFile
	open(ddobj_est_file2, "w") do f
		write(f,"$objDDJpSAA")
	end	
	
	time_obj_file = subDirName2 * jsaaTimeFile
	open(time_obj_file, "w") do f
		write(f,"$jsaaTime")
	end	
	
end
#*============================================
