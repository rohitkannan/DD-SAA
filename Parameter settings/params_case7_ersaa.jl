# parameters for the case study

# model and data generation files
const modelFile = "resource_20_30.jl"
const dataFile = "demand_model_iid.jl"

# parameters for estimating the true optimal objective value using MC sampling
const numMCScenarios = 1000
const numMCReplicates = 30

# number of covariates
const numCovariates = 10

# regression method
const regressionMethod = "ols"

# determine number of samples depending on the number of covariates
const numDataSamples = Int64[16,22,55,220,1100]

# degree of nonlinearity in demand model
const degree = Float64[1]

# scaling factor for the additive demand errors
const demand_errors_scaling = 20.0

# starting seed
const startingSeed = 1280
