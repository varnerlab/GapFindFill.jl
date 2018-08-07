# JuGapFill.jl
# redo GapFill in Julia
using JuMP, Gurobi
include("loadGAMS.jl")

#*************SETS********************************
reactionList = read_GAMS_format_file("gapfill/reactions.txt")
@assert(length(reactionList) == 2888)
modelRxnList = read_GAMS_format_file("gapfill/model_reactions.txt")
@assert(length(modelRxnList) == 1578)
revRxnList = read_GAMS_format_file("gapfill/reversible_rxns.txt")
@assert(length(revRxnList) == 1959)
databaseRxnList = read_GAMS_format_file("gapfill/database_reactions.txt")
@assert(length(databaseRxnList) == 1310)

compoundList = read_GAMS_format_file("gapfill/compounds.txt")
@assert(length(compoundList) == 1822)
cytCompoundList = read_GAMS_format_file("gapfill/cytosolic_metabolites.txt")
@assert(length(cytCompoundList) == 714)
extCompoundList = read_GAMS_format_file("gapfill/extracellular_metabolites.txt")
@assert(length(extCompoundList) == 678)
noProdCompList = read_GAMS_format_file("gapfill/no_production_metabolites.txt")
@assert(length(noProdCompList) == 8)
#*****************PARAMETERS*****************************
fluxUBList = read_GAMS_format_file("gapfill/upperbounds_on_fluxes.txt")
@assert(length(fluxUBList) == 2888)
fluxLBList = read_GAMS_format_file("gapfill/lowerbounds_on_fluxes.txt")
@assert(length(fluxLBList) == 2888)
stoiList = read_GAMS_format_file("gapfill/S_matrix.txt")
@assert(length(stoiList) > 8000)

#****************PARAMETERS PREPROCESSING***************************
# map rxn_name to index
rxnDict = Dict{String, Integer}()
for (id, name) in enumerate(reactionList)
  rxnDict[name] = id
end
# map compound_name to index
compoundDict = Dict{String, Integer}()
for (id, name) in enumerate(compoundList)
  compoundDict[name] = id
end
noRxn = length(reactionList)
noCompound = length(compoundList)

# isRev -- isRev[i] is true if rxn in rxnDict(rxn, i) is reversible
isRev = zeros(Bool, noRxn)
for el in revRxnList
  isRev[rxnDict[el]] = true
end
@assert(length(revRxnList) == sum(isRev))
# isMd -- isMd[i] is true if rxn in rxnDict(rxn, i) is in model
isMd = zeros(Bool, noRxn)
for el in modelRxnList
  isMd[rxnDict[el]] = true
end
@assert(length(modelRxnList) == sum(isMd))
# isDb -- isDb[i] is true if rxn in rxnDict(rxn, i) is in database
isDb = zeros(Bool, noRxn)
for el in databaseRxnList
  isDb[rxnDict[el]] = true
end
@assert(length(databaseRxnList) == sum(isDb))

# isCyt -- isCyt[i] is true if cp in compoundDict(cp, i) is in cytosol
isCyt = zeros(Bool, noCompound)
for el in cytCompoundList
  isCyt[compoundDict[el]] = true
end
@assert(length(cytCompoundList) == sum(isCyt))
# isExt -- isExt[i] is true if cp in compoundDict(cp, i) is in extracellular
isExt = zeros(Bool, noCompound)
for el in extCompoundList
  isExt[compoundDict[el]] = true
end
@assert(length(extCompoundList) == sum(isExt))
# no production metabolites list
noProdID = Array{Integer, 1}()
for el in noProdCompList
  push!(noProdID, compoundDict[el])
end

# get S matrix from stoiList
stoiMatrix = zeros(Float64, noCompound, noRxn)
for el in stoiList
  stoiMatrix[compoundDict[el[1]],rxnDict[el[2]]] = el[3]
end
println("size of S matrix: ", size(stoiMatrix))
# get upperbound on fluxes
fluxUB = zeros(Float64, noRxn)
for (name, val) in fluxUBList
  fluxUB[rxnDict[name]] = val
end
println("big-M in upperbound: ", maximum(fluxUB))
# get lowerbound on fluxes
fluxLB = zeros(Float64, noRxn)
for (name, val) in fluxLBList
  fluxLB[rxnDict[name]] = val
end
println("big-M in lowerbousnd: ", minimum(fluxLB))

# constants
epsilon = 0.001;
bigM = 1000;
nonZero = 1e-8;

#***************JuMP MODEL**********************************
function fixNoProductionMetabolite(Target::Integer)
  println("=================================")
  println("tackling no-production metabolite \"$(compoundList[Target])\" ......")
  solver = GurobiSolver()
  m = Model(solver=solver)
  # define variables using @variables blocks
  @variable(m, v[1:noRxn])  # fluxes
  # 1 if the reversibility of rnx i in model is relaxed
  @variable(m, x[i=1:noRxn; (isMd[i] && !isRev[i])], Bin)
  # 1 if rxn i from database is added
  @variable(m, y[i=1:noRxn; isDb[i]], Bin)
  # 1 if metabolite i is produced in reaction j
  @variable(m, w[1:noCompound, 1:noRxn], Bin)
  # define objective
  @objective(m, Min, sum(x)+sum(y))
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
  ID = Target
  for j in 1:noRxn
    if (abs(stoiMatrix[ID,j]) > nonZero)  # non-zero coefficient
      @constraint(m, stoiMatrix[ID,j]*v[j] >= epsilon - bigM*(1-w[ID,j]))
      @constraint(m, stoiMatrix[ID,j]*v[j] <= bigM*w[ID,j])
    end
  end
  # Binary Constraints (Eq. 6.12 on P124)
  @constraint(m, sum(w[ID,j] for j=1:noRxn if (
    abs(stoiMatrix[ID,j]) > nonZero)) >= 1)
  # Bound constraints (Eq. 6.13-15 on P124)
  for j in 1:noRxn
    if isMd[j] # Eq. 6.13-14
      @constraint(m, v[j] <= fluxUB[j])
      if isRev[j]  # Eq. 6.13
        @constraint(m, v[j] >= fluxLB[j])
      else # Eq. 6.14
        @constraint(m, v[j] >= -bigM*x[j])
      end
    else  # Eq. 6.15
      @constraint(m, v[j] >= fluxLB[j]*y[j])
      @constraint(m, v[j] <= fluxUB[j]*y[j])
    end
  end
  # Calling solver
  println("Calling Solver.............")
  t1 = time()
  status = solve(m)
  t2 = time()
  println("*****************Result Analysis")
  objVal = getobjectivevalue(m)
  println("Best objective: $(objVal)")
  # println("Solve Time by \"getsolvetime()\": $(getsolvetime(m))")
  println("Solve Time by \"time()\": $(t2-t1)")
  println("Solve Status: $(status)")
end

#*************SOLVING*********************
foreach(fixNoProductionMetabolite, noProdID)
println("=============\nDONE fixing")
