# Main file for solving the kNN-based reweighted SAA problem of Bertsimas and Kallus


using Distributions, StatsFuns, StatsBase, JuMP, Gurobi


const gurobi_env = Gurobi.Env()



# read these as inputs from the corresponding "instance file"
modRepNum = parse(ARGS[1])
degNum = parse(ARGS[2])
dataRepNum = parse(ARGS[3])
sampSizeNum = parse(ARGS[4])



#*========= MODEL INFORMATION ==========
const caseNum = 1
const paramsFile = "params_case" * string(caseNum) * "_knn.jl"
#*======================================



#*========= INCLUDE FILES ====================
include("maxParameters.jl")
include("randomSeeds.jl")

include(paramsFile)
# the next two data files are defined in paramsFile
include(modelFile)
include(dataFile)

include("genkNNScenarios.jl")
include("solveModel.jl")
include("estimateSolnCost.jl")
include("evaluateSoln.jl")
#*======================================



# set seed for reproducibility
srand(startingSeed)



# directory name for storing results
const baseDirName = "case" * string(caseNum) * "_knn/" * "mod_" * string(modRepNum) * "/" * "deg_" * string(degNum) * "/"
const subDirName = baseDirName * "rep_" * string(dataRepNum) * "/"
const subDirName2 = subDirName * "samp_" * string(sampSizeNum) * "/"
mkpath(subDirName2)



#*========= PRINT OPTIONS ====================
const storeResults = true

const infoFile = "ddsp.txt"
const modelDataFile = "model_data.txt"
const numSampFile = "num_samples.txt"

const knnObjFile = "knn_obj.txt"
const knnDDObjFile = "knn_ddobj.txt"
const knnTimeFile = "knn_time.txt"
const koptFile = "k_opt.txt"

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
		write(f,"numFolds: $numFolds \n")
		write(f,"minExponent: $minExponent \n")
		write(f,"maxExponent: $maxExponent \n\n")
		
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



#*========= SOLVE kNN-based SAA MODEL ====================

tic()


demand_scen_knn = generatekNNScenarios(demand_data,covariate_data,numFolds,minExponent,maxExponent,covariate_obs)

z_soln_knn, objDDknn = solveSAAModel(demand_scen_knn)

knnObjEstimates = estimateSolnQuality(z_soln_knn,covariate_obs,degree[degNum],numMCScenarios,numMCReplicates,coeff_true)


knnTime = toq()


#*========= STORE RESULTS ====================
if(storeResults)
	
	# write data to text file
	obj_est_file = subDirName2 * knnObjFile
	open(obj_est_file, "w") do f
		for i = 1:length(knnObjEstimates)
			write(f,"$(knnObjEstimates[i]) \n")
		end
	end	
	
	ddobj_est_file = subDirName2 * knnDDObjFile
	open(ddobj_est_file, "w") do f
		write(f,"$objDDknn")
	end	
	
	time_obj_file = subDirName2 * knnTimeFile
	open(time_obj_file, "w") do f
		write(f,"$knnTime")
	end	
		
	kopt_file = subDirName2 * koptFile
	open(kopt_file, "w") do f
		write(f,"$(size(demand_scen_knn,1)) \n")
	end
	
end
#*============================================
