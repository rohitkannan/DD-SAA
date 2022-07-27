# DD-SAA

**NOTE:** Codes correspond to Version 1 of the paper found <a href = "http://www.optimization-online.org/DB_FILE/2020/07/7932.pdf" target="_blank">here</a>. Codes corresponding to Version 2 of the paper <a href = "https://github.com/rohitkannan/DD-SAA/blob/master/Data-Driven%20SAA%20with%20Covariate%20Information%20(R1).pdf" target="_blank">here</a> will be updated soon.

Julia codes for the paper "Data-Driven Sample Average Approximation with Covariate Information". Results in the paper were generated using Julia 0.6.4, JuMP 0.18.5, Gurobi 8.1.0, and GLMNet 0.3.0.

To generate results for the FI-SAA, N-SAA, ER-SAA, J-SAA & J+-SAA, or kNN-SAA methods, run the "main_fullinf_saa.jl", "main_nsaa.jl", "main_ersaa.jl", "main_jsaa.jl", or "main_knn.jl" file, respectively.

Parameter settings for each of the case studies are provided in the "params_case*" files within the _Parameter settings_ folder, and the instances for each case study are provided in the "instances_case*" files within the _Instances_ folder. To run each of the "main" files, set the variable "caseNum" to the appropriate value, and write a script to read the instances from the corresponding "instances" file.
