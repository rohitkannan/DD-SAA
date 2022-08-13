# parameters for the kNN-based reweighted SAA approach of Bertsimas and Kallus with d_x = 3 and omega = 1

# model and data generation files
const modelFile = "resource_20_30.jl"
const dataFile = "demand_model.jl"

# parameters for estimating the true optimal objective value using MC sampling
const numMCScenarios = 1000
const numMCReplicates = 30

# number of covariates
const numCovariates = 3

# determine number of samples depending on the number of covariates
const numDataSamples = Int64[6,8,20,80,400]

# degree of nonlinearity in demand model
const degree = Float64[1,0.5,2]

# heteroscedasticity level
const het_level = 1

# scaling factor for the additive demand errors
const demand_errors_scaling = 5.0

# number of folds for K-fold CV to determine optimal k in kNN
const numFolds = 5

# minimum and maximum exponents to determine optimal k in kNN
const minExponent = 0.1
const maxExponent = 0.9

# starting seed
const startingSeed = 1280
