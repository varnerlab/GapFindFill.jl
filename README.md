# GFFJ: Gap Finding and Filling in Julia
An implementation of [Gap Finding & Filling](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-212) in Julia. 

## Requirement
In order to use SEML, the user needs to [install Julia](https://julialang.org/downloads/platform.html) first. This version is compatible with [Julia v1.1.0](https://julialang.org/downloads/index.html).

[Gurobi.jl](https://github.com/JuliaOpt/Gurobi.jl) is used in GFFJ as the default solver for both gap finding and filling. [Gurobi](http://www.gurobi.com/) provides free academic license for non-commercial use. 
Users can also choose to use [GLPK.jl](https://github.com/JuliaOpt/GLPK.jl) or [CPLEX.jl](https://github.com/JuliaOpt/CPLEX.jl) for gap finding. GLPK is free for all users, while CPLEX provides free academic license. But note that comparing to Gurobi and CPLEX, it takes GLPK much more time and memory to solve the same problem. 

For gap filling, Gurobi is set as the only solver since it worked way better than the other two in our tests. 

## Usage
Within [Julia](http://http://julialang.org), press "__]__" to enter __pkg>__ mode. 
To install GFFJ, issue 

```julia
add https://github.com/varnerlab/Julia_GapFill_Repository.git
```
To test the GFFJ installation use:

```julia
test GFFJ 
```
which runs two test examples from the __test__ directory. 

To delete GFFJ package use the command:

```julia
rm GFFJ
```

To use GFFJ in your project simply issue the command:

```julia
using GFFJ
```

Two interfaces are provided, "__find_gaps()__" and "__fill_gaps_min()__". 

The "__find_gaps()__" interface:
```julia
function find_gaps(isRev::Array{Bool}, isCyt::Array{Bool}, isExt::Array{Bool},
    stoiMatrix::Array{Float64}, fluxLB::Array{Float64}, fluxUB::Array{Float64};
    epsilon::Float64 =0.001, bigM::Float64 =1000.0, nonZero::Float64 =1e-8,
    solver::Module=Gurobi)
```
Inputs description: 

Argument | Required | Description 
:--- | :--- | :---
isRev | yes | true if corresponding reaction is reversible
isCyt | yes | true if corresponding metabolite is in cytosol
isExt | yes | true if corresponding metabolite is in extracellular compartment;
stoiMatrix | yes | stoichiometric matrix, \|compounds\| * \|reactions\|
fluxLB | yes | flux lower bound;
fluxUB | yes | flux upper bound;
epsilon | optional | minimum amount to be considered active;
bigM | optional | constant used in MILP model;
nonZero | optional | minimum stoichiometric coefficient to be considered non-zero;
solver | optional | CPLEX, Gurobi or GLPK.

Outputs description: 

Argument | Description 
:--- | :--- 
m | the JuMP model;
objVal | objective value, i.e. number of non-blocked compounds;
status | termination status;
binX | 1 if corresponding compound is non-blocked;


The "__fill_gaps_min()__" interface: 
```julia 
function fill_gaps_min(isMd::Array{Bool}, isDb::Array{Bool}, isRev::Array{Bool},
    isCyt::Array{Bool}, isExt::Array{Bool}, noProdID::Array{Integer},
    stoiMatrix::Array{Float64}, fluxLB::Array{Float64}, fluxUB::Array{Float64};
    epsilon::Float64 =0.001, bigM::Float64 =1000.0, nonZero::Float64 =1e-9)
```
Inputs description:  

Argument | Required | Description 
:--- | :--- | :---
isMd | yes | true if corresponding reaction is in the model;
isDb | yes | true if corresponding reaction is in the database, i.e., not in the model;
isRev | yes | true if corresponding reaction is reversible
isCyt | yes | true if corresponding metabolite is in cytosol
isExt | yes | true if corresponding metabolite is in extracellular compartment;
noProdID | yes | indices of no-production-metabolites;
stoiMatrix | yes | stoichiometric matrix, \|compounds\| * \|reactions\|
fluxLB | yes | flux lower bound;
fluxUB | yes | flux upper bound;
epsilon | optional | minimum amount to be considered active;
bigM | optional | constant used in MILP model;
nonZero | optional | minimum stoichiometric coefficient to be considered non-zero.

Outputs description: 

Argument | Description 
:--- | :--- 
results | a dictionary, each key is an index of no-production-mebanolites, each value is a Set containing all possible gap filling ways with minimum number of added reactions.


## Reference:
- Maranas, Costas D., and Ali R. Zomorrodi. Optimization methods in metabolic networks. John Wiley & Sons, 2016.
- Kumar, Vinay Satish, Madhukar S. Dasika, and Costas D. Maranas. "Optimization based automated curation of metabolic reconstructions." BMC bioinformatics 8.1 (2007): 212.
- [GapFind/GapFill](http://www.maranasgroup.com/software.htm) in [GAMS](https://www.gams.com/)


## Support or Contact
Having trouble at installation or function? Feel free to contact the [authors](https://github.com/varnerlab) or [authors](https://www.cheme.cornell.edu/faculty-directory/jeffrey-d-varner).
