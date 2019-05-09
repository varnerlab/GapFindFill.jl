using GFFJ, Test, Gurobi
include("loadGAMS.jl")

println("======Tesing Gap Filling===")
println("======loading data===")
#*************SETS********************************
compoundList = read_GAMS_format_file("gapfill/compounds.txt")
@assert(length(compoundList) == 1822)
reactionList = read_GAMS_format_file("gapfill/reactions.txt")
# reactionList = modelRxnList + databaseRxnList
@assert(length(reactionList) == 2888)
modelRxnList = read_GAMS_format_file("gapfill/model_reactions.txt")
@assert(length(modelRxnList) == 1578)
databaseRxnList = read_GAMS_format_file("gapfill/database_reactions.txt")
@assert(length(databaseRxnList) == 1310)
revRxnList = read_GAMS_format_file("gapfill/reversible_rxns.txt")
@assert(length(revRxnList) == 1959)
noProdCompList = read_GAMS_format_file("gapfill/no_production_metabolites.txt")
@assert(length(noProdCompList) <= 8)
cytCompoundList = read_GAMS_format_file("gapfill/cytosolic_metabolites.txt")
@assert(length(cytCompoundList) == 714)
extCompoundList = read_GAMS_format_file("gapfill/extracellular_metabolites.txt")
@assert(length(extCompoundList) == 678)

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
# number of reactions & compounds
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


isExt = zeros(Bool, noCompound)

println("======calling fill_gaps_min===")
# test fill_gaps_min
results = fill_gaps_min(isMd, isDb, isRev, isCyt, isExt, noProdID,
    stoiMatrix, fluxLB, fluxUB)
println("======results report===")
for (key, val) in results
  println("To fix metabolite $(compoundList[key])")
  for el in val
    println(join([reactionList[el2] for el2 in el], " "))
  end
  println("\n")
end

println("======Gap Filling Test Passed===")
