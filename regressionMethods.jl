# Routine for estimating Lasso coefficients using GLMNet


function getLassoCoeff(outputs::Array{Float64,2},inputs::Array{Float64,2})

	numDependentVar::Int64 = size(outputs,2)
	numIndependentVar::Int64 = size(inputs,2)
	
	if(size(outputs,1) != size(inputs,1))
		throw(ErrorException("Number of output samples != number of input samples"))
	end

	coeff = zeros(Float64,numIndependentVar,numDependentVar)
	
	for d = 1:numDependentVar
		path_cv = glmnetcv(inputs[:,2:numIndependentVar],outputs[:,d]);
		path_ind::Int64 = indmin(path_cv.meanloss)
		coeff[2:numIndependentVar,d] = path_cv.path.betas[:,path_ind]
	end

	tmpMat = outputs - inputs[:,2:numIndependentVar]*coeff[2:numIndependentVar,:]
	for d = 1:numDependentVar
		coeff[1,d] = mean(tmpMat[:,d])
	end
	
	return coeff
end
