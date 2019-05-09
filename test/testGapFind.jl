using GFFJ, Test, Gurobi
include("loadGAMS.jl")

println("======Tesing Gap Finding===")
println("======loading data===")
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

println("======calling find_gaps===")
#*********Calling function*********************
model, objVal, status, solveTime, xValue = find_gaps(isRev, isCyt, isExt,
  stoiMatrix, fluxLB, fluxUB, solver=Gurobi)
println("======results report===")
println("Best objective: $(objVal)")
println("Solve Time by \"time()\": $(solveTime)")
println("Solve Status: $(status)")

#***************VERIFICATION*********************
println("======VERIFICATION===")
stdResultList = read_GAMS_format_file("gapfind/sample_results.txt")
stdResult = Set(stdResultList)
@assert(length(stdResult) == (noCompound - objVal))
println(sum(xValue))
@assert(objVal == sum(xValue))
try
  for (id, val) in enumerate(xValue)
    if val == 0
      @assert(in(compoundList[id], stdResult))
    end
  end
  println("Consistent results")
  println("======Gap Finding Test Passed===")
catch e
  println("Inconsistent results")
end
