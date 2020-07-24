# parameters for the case study

# model and data generation files
const modelFile = "resource_20_30.jl"
const dataFile = "demand_model_iid.jl"

# parameters for estimating the true optimal objective value using SAA
const numMCScenarios = 1000

# degree of nonlinearity in demand model
const degree = Float64[1,0.5,2]

# scaling factor for the additive demand errors
const demand_errors_scaling = 5.0

# starting seed
const startingSeed = 1280
