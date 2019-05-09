# GFFJ Interface

#=
Set up MILP model for Gap-Find
Input:
  isRev: true if corresponding reaction is reversible;
  isCyt: true if corresponding metabolite is in cytosol;
  isExt: true if corresponding metabolite is in extracellular compartment;
  stoiMatrix: stoichiometric matrix, #compounds * #reactions
  fluxLB: flux lower bound;
  fluxUB: flux upper bound;
  solver: CPLEX, Gurobi, GLPK;
  epsilon: minimum amount to be considered active;
  bigM: constant used in MILP model;
  nonZero: minimum stoichiometric coefficient to be considered non-zero.
Output:
  m: the JuMP model;
  objVal: objective value, i.e. # non-blocked compounds;
  status: termination status;
  binX: 1 if corresponding compound is non-blocked;
=#
function find_gaps(isRev::Array{Bool}, isCyt::Array{Bool}, isExt::Array{Bool},
    stoiMatrix::Array{Float64}, fluxLB::Array{Float64}, fluxUB::Array{Float64};
    epsilon::Float64 =0.001, bigM::Float64 =1000.0, nonZero::Float64 =1e-8,
    solver::Module=Gurobi)
  t1 = time()
  noCompound, noRxn = size(stoiMatrix)
  m = Model(with_optimizer(solver.Optimizer))
  # define variables using @variable
  @variable(m, fluxLB[i] <= v[i=1:noRxn] <= fluxUB[i])  # fluxes
  @variable(m, x[1:noCompound], Bin)  # 1 if metabolite i could be produced
  @variable(m, w[1:noCompound, 1:noRxn], Bin)  # 1 if metabolite i is produced in reaction j
  # define objective
  @objective(m, Max, sum(x))
  # define constraints
  # Mass Balance
  for i = 1:noCompound
    if (isCyt[i] && !isExt[i])  # cytosl but not extracellular
      @constraint(m, sum(stoiMatrix[i,j]*v[j] for j=1:noRxn) >= 0)
    elseif (!isCyt[i] && !isExt[i])  # other compartments but not extracellular
      @constraint(m, sum(stoiMatrix[i,j]*v[j] for j=1:noRxn) == 0)
    end
  end
  # Production Constrains for (ir)reversible reactions
  for i = 1:noCompound
    for j = 1:noRxn
      if (abs(stoiMatrix[i,j]) > nonZero)  # non-zero coefficient
        @constraint(m, stoiMatrix[i,j]*v[j] >= epsilon - bigM*(1-w[i,j]))
        @constraint(m, stoiMatrix[i,j]*v[j] <= bigM*w[i,j])
      end
    end
  end
  # Binary Constraints
  for i = 1:noCompound
    @constraint(m, sum(w[i,j] for j=1:noRxn if (
      (!isRev[j] && (stoiMatrix[i,j] > nonZero)) ||
      (isRev[j] && (abs(stoiMatrix[i,j]) > nonZero)))) >= x[i])
  end
  # Calling solver
  JuMP.optimize!(m)
  t2 = time()
  # get some results
  objVal = objective_value(m)
  status = termination_status(m)
  solveTime = t2 - t1
  binX = value.(x)
  return m, objVal, status, solveTime, binX
end




#=
Set up MILP model for Gap-Fill
Input:
  isMd: true if corresponding reaction is in the model;
  isDb: true if corresponding reaction is in the database,
        i.e. not in the model;
  isRev: true if corresponding reaction is reversible;
  isCyt: true if corresponding metabolite is in cytosol;
  isExt: true if corresponding metabolite is in extracellular compartment;
  noProdID: indices of no-production-metabolites
  stoiMatrix: stoichiometric matrix, #compounds * #reactions
  fluxLB: flux lower bound;
  fluxUB: flux upper bound;
  epsilon: minimum amount to be considered active;
  bigM: constant used in MILP model;
  nonZero: minimum stoichiometric coefficient to be considered non-zero.
Output:
  results: a dictionary, each key is an index of no-production-mebanolites,
    each value is a Set containing all possible gap filling ways with minimum
    number of added reactions.
=#
function fill_gaps_min(isMd::Array{Bool}, isDb::Array{Bool}, isRev::Array{Bool},
    isCyt::Array{Bool}, isExt::Array{Bool}, noProdID::Array{Integer},
    stoiMatrix::Array{Float64}, fluxLB::Array{Float64}, fluxUB::Array{Float64};
    epsilon::Float64 =0.001, bigM::Float64 =1000.0, nonZero::Float64 =1e-9)
  solver=Gurobi
  results = Dict{Int,Any}()
  for i in noProdID
    m, objVal, status, solveTime, yresult = fix_metabolite(
        isMd, isDb, isRev, isCyt, isExt, i,
        stoiMatrix, fluxLB, fluxUB,
        solver, epsilon, bigM, nonZero)
    println(objVal)
    println(yresult)
    # store result anyway (maybe no result)
    if haskey(results, i)
      push!(results[i], yresult)
    else
      results[i] = Set()
      push!(results[i], yresult)
    end
    if objVal > 0.5
      # apply integer cuts
      Z = objVal
      while (objVal <= Z) && (objVal > 0.5)
        # println("=============current Z value $Z")
        m, objVal, status, solveTime, yresult = fix_metabolite(
            isMd, isDb, isRev, isCyt, isExt, i,
            stoiMatrix, fluxLB, fluxUB,
            solver, epsilon, bigM, nonZero, Z, results[i])
        println(objVal)
        println(yresult)
        # store valid result only
        if (objVal <= Z) && (objVal > 0.5)
          push!(results[i], yresult)
        end
      end
    end
  end
  return results
end


function fix_metabolite(isMd, isDb, isRev, isCyt, isExt, ID,
    stoiMatrix, fluxLB, fluxUB, solver, epsilon, bigM, nonZero, intCuts...)
  t1 = time()
  println("---tackling no-production metabolite \"$ID\"..")
  noCompound, noRxn = size(stoiMatrix)
  m = Model(with_optimizer(solver.Optimizer))
  # define variables using @variables blocks
  @variable(m, v[1:noRxn])  # fluxes
  # # 1 if the reversibility of rnx i in model is relaxed
  # @variable(m, x[i=1:noRxn; (isMd[i] && !isRev[i])], Bin)
  # 1 if rxn i is relaxed or added
  @variable(m, y[i=1:noRxn; (isDb[i] || (isMd[i] && !isRev[i]))], Bin)
  # 1 if metabolite ID is produced in reaction j
  @variable(m, w[1:noRxn], Bin)
  # @variable(m, w[1:noCompound, 1:noRxn], Bin)

  # define objective
  @objective(m, Min, sum(y))
  # define constraints
  # Mass Balance (Eq. 6.9 & 6.10 on P124)
  for i = 1:noCompound
    if (isCyt[i] && !isExt[i])
      # cytosol but not extracellular Eq. 6.9
      @constraint(m, sum(stoiMatrix[i,j]*v[j] for j=1:noRxn) >= 0)
    elseif (!isCyt[i] && !isExt[i])
      # other compartments but not extracellular, Eq. 6.10
      @constraint(m, sum(stoiMatrix[i,j]*v[j] for j=1:noRxn) == 0)
    end
  end
  # Production Constrains for (ir)reversible reactions (Eq. 6.11 on P124)
  for j in 1:noRxn
    # keep this pre-condition to reduce # of constraints
    if (abs(stoiMatrix[ID,j]) > nonZero)  # non-zero coefficient
      @constraint(m, stoiMatrix[ID,j]*v[j] >= epsilon - bigM*(1-w[j]))
      @constraint(m, stoiMatrix[ID,j]*v[j] <= bigM*w[j])
    end
  end
  # Binary Constraints (Eq. 6.12 on P124)
  @constraint(m, sum(w[j] for j=1:noRxn if
    (abs(stoiMatrix[ID,j]) > nonZero)) >= 1)
  # Bound constraints (Eq. 6.13-15 on P124)
  for j in 1:noRxn
    if isMd[j] # Eq. 6.13-14
      @constraint(m, v[j] <= fluxUB[j])
      if isRev[j]  # Eq. 6.13
        @constraint(m, v[j] >= fluxLB[j])
      else # Eq. 6.14
        @constraint(m, v[j] >= -bigM*y[j])
      end
    else  # Eq. 6.15
      @constraint(m, v[j] >= fluxLB[j]*y[j])
      @constraint(m, v[j] <= fluxUB[j]*y[j])
    end
  end
  # add integer cuts
  if length(intCuts) != 0
    for cut in intCuts[2]
      @constraint(m, sum(y[j] for j in cut) <= intCuts[1]-1)
    end
  end
  # Calling solver
  JuMP.optimize!(m)
  t2 = time()
  # get some results
  status = termination_status(m)
  solveTime = t2 - t1
  yresult = Set{Int64}()
  if (status == MOI.OPTIMAL) || (status == MOI.TIME_LIMIT && has_values(m))
    objVal = objective_value(m)
    # extract y
    for i in 1:noRxn
      if isDb[i] || (isMd[i] && !isRev[i])
        if value(y[i]) == 1
          push!(yresult, i)
        end
      end
    end
  else
    println("The model was not solved correctly.")
    objVal = -1
  end
  return m, objVal, status, solveTime, yresult
end
