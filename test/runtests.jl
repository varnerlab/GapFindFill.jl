using GapFindFill, Test
# using JuMP, Gurobi, CPLEX, GLPK

println("Starting tests....................")
@time begin include("testGapFind.jl") end
@time begin include("testGapFill.jl") end
