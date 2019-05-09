module GFFJ

using JuMP, Gurobi, CPLEX, GLPK
include("GFFJInterface.jl")

# Export the interface function -
export find_gaps, fill_gaps_min

end # module
