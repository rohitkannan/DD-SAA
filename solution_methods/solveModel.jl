# Solve a SAA of the resource allocation model given scenarios of the demands
function solveSAAModel(demand_scen::Array{Float64,2})

	numScenarios::Int64 = size(demand_scen,1)
	
	mod = Model(solver=GurobiSolver(gurobi_env,Presolve=0,OutputFlag=0,Threads=maxNumThreads))

	@variable(mod, z[1:numResources] >= 0)
	@variable(mod, v[1:numResources,1:numCustomers,1:numScenarios] >= 0)
	@variable(mod, w[1:numCustomers,1:numScenarios] >= 0)

	
	@objective(mod, Min, sum(costVector[i]*z[i] for i=1:numResources) + (1.0/numScenarios)*sum(recourseCostCoeff[j]*w[j,s] for j=1:numCustomers for s=1:numScenarios))
	
	@constraint(mod, [i=1:numResources, s=1:numScenarios], sum(v[i,j,s] for j = 1:numCustomers) <= rho[i]*z[i])
	@constraint(mod, [j=1:numCustomers, s=1:numScenarios], sum(mu[i,j]*v[i,j,s] for i = 1:numResources) + w[j,s] >= demand_scen[s,j])

	
	status = solve(mod)
	z_soln = getvalue(z)
	objValue_soln::Float64 = getobjectivevalue(mod)
	
	return z_soln, objValue_soln
end