# JuGapFind.jl
# redo GapFind in Julia
using JuMP, Gurobi
include("loadGAMS.jl")

#*************SETS********************************
reactionList = read_GAMS_format_file("gapfind/reactions.txt")
@assert(length(reactionList) == 2383)
# foreach(println, reactionList)
revRxnList = read_GAMS_format_file("gapfind/reversible_reactions.txt")
@assert(length(revRxnList) == 852)
# foreach(println, revRxnList)
compoundList = read_GAMS_format_file("gapfind/compounds.txt")
@assert(length(compoundList) == 1668)
# foreach(println, compoundList)
cytCompoundList = read_GAMS_format_file("gapfind/cytosolic_compounds.txt")
@assert(length(cytCompoundList) == 951)
# foreach(println, cytCompoundList)
extCompoundList = read_GAMS_format_file("gapfind/extracellular_compounds.txt")
@assert(length(extCompoundList) == 299)
# foreach(println, extCompoundList)

#*****************PARAMETERS*****************************
fluxUBList = read_GAMS_format_file("gapfind/upperbound_on_fluxes.txt")
@assert(length(fluxUBList) == 2383)
# foreach(println, fluxUBList)
fluxLBList = read_GAMS_format_file("gapfind/lowerbound_on_fluxes.txt")
@assert(length(fluxLBList) == 2383)
# foreach(println, fluxLBList)
stoiList = read_GAMS_format_file("gapfind/S_matrix.txt")
@assert(length(stoiList) == 9329)
# foreach(println, stoiList)
# constants
epsilon = 0.001;
bigM = 1000;
nonZero = 1e-8;

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
# get S matrix from stoiList
stoiMatrix = zeros(Float64, noCompound, noRxn)
for elm in stoiList
  stoiMatrix[compoundDict[elm[1]],rxnDict[elm[2]]] = elm[3]
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
println("big-M in lowerbound: ", minimum(fluxLB))

#***************JuMP MODEL**********************************
solver = GurobiSolver()
m = Model(solver=solver)
# define variables using @variables blocks
# @variable(m, fluxLB[i] <= v[i=1:noRxn] <= fluxUB[i])  # fluxes
# @variable(m, x[1:noCompound], Bin)  # 1 if metabolite i could be produced
# @variable(m, xp[1:noCompound, 1:noRxn], Bin)  # 1 if metabolite i is produced in reaction j
@variables m begin
  v[i=1:noRxn], (lowerbound = fluxLB[i], upperbound = fluxUB[i])  # fluxes
  x[1:noCompound], Bin  # 1 if metabolite i could be produced
  w[1:noCompound, 1:noRxn], Bin  # 1 if metabolite i is produced in reaction j
end
# define objective
@objective(m, Max, sum(x))
# define constraints
# Mass Balance (Eq. 6.4 & 6.5 on P121)
for i = 1:noCompound
  if (isCyt[i] && !isExt[i])  # cytosl but not extracellular
    @constraint(m, sum(stoiMatrix[i,j]*v[j] for j=1:noRxn) >= 0)
  elseif (!isCyt[i] && !isExt[i])  # other compartments but not extracellular
    @constraint(m, sum(stoiMatrix[i,j]*v[j] for j=1:noRxn) == 0)
  end
end
# Production Constrains for (ir)reversible reactions (Eq. 6.1 on P121)
for i = 1:noCompound
  for j = 1:noRxn
    if (abs(stoiMatrix[i,j]) > nonZero)  # non-zero coefficient
      @constraint(m, stoiMatrix[i,j]*v[j] >= epsilon - bigM*(1-w[i,j]))
      @constraint(m, stoiMatrix[i,j]*v[j] <= bigM*w[i,j])
    end
  end
end
# Binary Constraints (Eq. 6.2 on P121)
for i = 1:noCompound
  @constraint(m, sum(w[i,j] for j=1:noRxn if (
    (!isRev[j] && (stoiMatrix[i,j] > nonZero)) ||
    (isRev[j] && (abs(stoiMatrix[i,j]) > nonZero)))) >= x[i])
end
# Calling solver
println("===========Calling Solver=============")
t1 = time()
status = solve(m)
t2 = time()
println("===========Result Analysis=============")
objVal = getobjectivevalue(m)
println("Best objective: $(objVal)")
println("Solve Time by \"getsolvetime()\": $(getsolvetime(m))")
println("Solve Time by \"time()\": $(t2-t1)")
println("Solve Status: $(status)")
#*************VERIFICATION*********************
println("==========VERIFICATION====================")
stdResultList = read_GAMS_format_file("gapfind/sample_results.txt")
stdResult = Set(stdResultList)
@assert(length(stdResult) == (noCompound - objVal))
xValue = getvalue(x)
println(sum(xValue))
@assert(objVal == sum(xValue))
for (id, val) in enumerate(xValue)
  if val == 0
    @assert(in(compoundList[id], stdResult))
  end
end
println("Consistent results")

#******Identify Root no-production metabolites****************
println("====Root no-production metabolites=======")
rootNoProd = Array{Any, 1}()
for i in 1:noCompound
  isRootNoP = true
  for j in 1:noRxn
    c1 = stoiMatrix[i,j] > nonZero
    c2 = isRev[j] && abs(stoiMatrix[i,j]) > nonZero
    if (c1 || c2)
      isRootNoP = false
      break
    end
  end
  if isRootNoP
    name = compoundList[i]
    push!(rootNoProd, name)
    println(name)
  end
end
println("total number of root no-production metabolites: $(length(rootNoProd))")
println("total number of no-production metabolites: $(noCompound - objVal)")

# generator expression
