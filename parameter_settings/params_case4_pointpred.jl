# parameters for point prediction + Lasso with d_x = 10 and omega = 1

# model and data generation files
const modelFile = "resource_20_30.jl"
const dataFile = "demand_model.jl"

# parameters for estimating the true optimal objective value using MC sampling
const numMCScenarios = 1000
const numMCReplicates = 30

# number of covariates
const numCovariates = 10

# regression method
const regressionMethod = "lasso"

# determine number of samples depending on the number of covariates
const numDataSamples = Int64[13,14,16,22,33]

# degree of nonlinearity in demand model
const degree = Float64[1,0.5,2]

# heteroscedasticity level
const het_level = 1

# scaling factor for the additive demand errors
const demand_errors_scaling = 5.0

# starting seed
const startingSeed = 1280
