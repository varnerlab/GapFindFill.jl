# GFFJ.jl Documentation 
- [Statement of need](https://github.com/varnerlab/GFFJ/blob/master/README.md#statement-of-need)
- [Installation instruction](https://github.com/varnerlab/GFFJ/blob/master/README.md#installation-instruction)
- [Example](https://github.com/varnerlab/GFFJ/blob/master/README.md#example)
- [API documentation](https://github.com/varnerlab/GFFJ/blob/master/README.md#api-documentation)
- [Reference](https://github.com/varnerlab/GFFJ/blob/master/README.md#reference)
- [Support or Contact](https://github.com/varnerlab/GFFJ/blob/master/README.md#support-or-contact)


## Statement of need 
The current implementation of *GapFind* and *GapFill* is in GAMS, which charges a significant amount of license fee from each single user, even though many solvers are free for academic purpose.
To promote the usage of this computational tool, we developed this open-source Julia package, GFFJ.jl, to enable researchers to use *GapFind* and *GapFill* for free by harnessing the power of academic free solvers provided by [Gurobi](https://www.gurobi.com/) and [IBM](https://www.ibm.com/analytics/cplex-optimizer).
GFFJ.jl is implemented in Julia and makes use of the high-level interface [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl).
JuMP is a domain-specific modeling language for mathematical optimization embedded in Julia. 
With JuMP, it it easier for users to specify and call different optimizers to solve optimization problems in GFFJ.jl than using interfaces provided by solvers directly. 
Built upon the generic high-level programming language Julia, users can embed GFFJ.jl in their complex work flows to simplify task processing. While GAMS, as a specific optimization tool, does not provide support of processing other tasks, nor being able to be integrated with other programming languages. 


## Installation instruction 
**Requirement**.
In order to use SEML, the user needs to [install Julia](https://julialang.org/downloads/platform.html) first. This version is compatible with [Julia v1.1](https://julialang.org/downloads/oldreleases.html).
[Gurobi.jl](https://github.com/JuliaOpt/Gurobi.jl) is used in GFFJ as the default solver for both gap finding and filling. [Gurobi](http://www.gurobi.com/) provides free academic license for non-commercial use. 
Users can also choose to use [GLPK.jl](https://github.com/JuliaOpt/GLPK.jl) or [CPLEX.jl](https://github.com/JuliaOpt/CPLEX.jl) for gap finding. GLPK is free for all users, while CPLEX provides free academic license. But note that comparing to Gurobi and CPLEX, it takes GLPK much more time and memory to solve the same problem. 
For gap filling, Gurobi is set as the only solver since it worked way better than the other two in our tests. 

For [Julia v1.0](https://julialang.org/downloads/) users, any package compatibility issues while testing GFFJ.jl can be resolved by pinning CPLEX.jl and Gurobi.jl to a specific version by running following commands in `pkg>` mode: 
```julia
pin Gurobi@0.6.0
pin CPLEX@0.5.0 
```

For [Julia v1.2](https://julialang.org/downloads/) users, GFFJ.jl will also be compatible with Julia v1.2 once [CPLEX.jl](https://github.com/JuliaOpt/CPLEX.jl) is updated for Julia v1.2. 

**Installation**.
Within [Julia](http://http://julialang.org), press `]` to enter `pkg>` mode. 
To install GFFJ, issue 
```julia
add https://github.com/varnerlab/GFFJ.git
```
To use GFFJ in your project simply issue the command:
```julia
using GFFJ
```

**Test**. 
To test the GFFJ installation use the following command in `pkg>` mode:
```julia
test GFFJ 
```
which runs two test examples under [test](https://github.com/varnerlab/GFFJ/tree/master/test) directory.
The expected outcomes from this test are illustrated in [Example](https://github.com/varnerlab/GFFJ/blob/master/README.md#example) section. 

**Uninstallation**.
To delete GFFJ package use the following command in `pkg>` mode:
```julia
rm GFFJ
```


## Example 
Two examples are provided under [test](https://github.com/varnerlab/GFFJ/tree/master/test) showing how to use GFFJ.jl to solve two problems in [Manaras paper]((https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-212)). 
[testGapFind.jl](https://github.com/varnerlab/GFFJ/blob/master/test/testGapFind.jl) and 
[testGapFill.jl](https://github.com/varnerlab/GFFJ/blob/master/test/testGapFill.jl) demonstrate how to set up *GapFind* and *GapFill* models, respectively. 
We reported our experimental results here for users reference. All experiments were run on an Intel Core i7-6700 CPU with Ubuntu 10.04.

For *find_gaps*, the testing example contains 1668 compounds and 2383 reactions, which is of size of real problems. 
The expected outcome is to find 115 blocked metabollites in the network, namely, 1553 non-blocked metabolites. The following table shows running time comparison between GAMS and GFFJ.jl on gap finding.

Software | Solver | Running time (s) 
:--- | :--- | :---
GAMS | CPLEX | 3.3 
GAMS | Gurobi | 0.6 
GFFJ.jl | CPLEX | 68.4
GFFJ.jl | Gurobi | 64.5

GAMS is pretty fast in solving large size problems like the testing example, but GFFJ.jl finished the job in less than 2 mins, which is also acceptable for real problems as large as this example. 
Using either Gurobi or CPLEX in GFFJ.jl did not make much difference in running time, while GLPK was unable to solve the testing example within 1 hr.
Thus, although *find_gaps* allows users to specify GLPK as the solver, it is recommended only for small scale problems.

For *fill_gaps_min*, the testing example contains 1822 compounds and 2888 reactions, which is also of size of real problems. Five no-production-metabolites were chosen as testing cases. The expected outcomes are summarized in the following table. 

No-production-metabolite | Solution 1 | Solution 2
:--- | :--- | :---
2doxg6p[c] | EX_2DOXG6P(e) & 2DOXG6P_t | 
2dglc[c] | EX_2DGLC(e) & 2DGLC_t | EX_2DOXG6P(e) & 2DOXG6P_t
alatrna[c] | TRNAALA_t & EX_TRNAALA(e) | ALATRNA_t & EX_ALATRNA(e)
2dr5p[c] | DRIB_t & EX_DRIB(e) | 2DR5P_t & EX_2DR5P(e)
4gudbutn[c] | 4GUDBD_t & EX_4GUDBD(e) | 4GUDBUTN_t & EX_4GUDBUTN(e)

The following table shows running time comparison between GAMS and GFFJ.jl on gap filling. 

Software | Solver | Running time (s) 
:--- | :--- | :---
GAMS | Gurobi | 6.2
GFFJ.jl | Gurobi | 28.7

The difference between GAMS and GFFJ.jl is smaller than solving gap finding tasks.
This is largely because GFFJ.jl got a speed-up from just-in-time feature of Julia as similar problems were solved repeatedly.


## API documentation 
Two interfaces are provided, `find_gaps()` and `fill_gaps_min()` for gap finding and filling respectively. 

The `find_gaps()` interface:
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

The `fill_gaps_min()` interface: 
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
Having trouble at installation or function? Feel free to contact [VarnerLab](https://github.com/varnerlab) or [authors](https://www.cheme.cornell.edu/faculty-directory/jeffrey-d-varner).
