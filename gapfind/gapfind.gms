*****************************************************************************
*                               GAPFIND                                     *
*                               -------                                     *
* Related publication: PMID:17584497                                        *
*    - This code identifies the no-production metabolites for the           *
*      iAF1260 model of E. coli as an example                               *
*    - The code can be easily modified for any other metabolic model        *
*    - The results of this code are stored in a file called:                *
*      no_production_metabolites.txt. A sample of results is given in the   *
*      in the file sample_results.txt                                       *
*                                                                           *
* This code is for demonstration purposes only and the results may not      *
* have any biological implications                                          *
*                                                                           *
* -- Ali R. Zomorrodi - Chemical & Biological Systems Optimization Lab      *
*****************************************************************************

$INLINECOM /*  */

OPTION decimals = 8
       sysout = off
       solprint = on
       reslim = 100000
       iterlim = 10000000
       domlim = 10
       limcol = 10
       limrow = 10
       optca = 0.0
       optcr = 0.0
       work = 10000000
       mip = cplex;



SETS

	i   Name of the metabolites in the model
$include "compounds.txt"
/* or, include here your own file */

	j   Name of the reactions
$include "reactions.txt"
/* or, include here your own file */

	rev(j)  Reversible reactions (including exchange rxns)
$include "reversible_reactions.txt"
* use string as index? and sparse array?
/* or, include here your own file */

	cytosol(i)	Cytosolic metabolites
$include "cytosolic_compounds.txt"
/* or, include here your own file */

	extracellular(i)  Extracellular metabolites
$include "extracellular_compounds.txt"
/* or, include here your own file */
;

PARAMETERS
  UB(j)     Upper bound on fluxes
$include "upperbound_on_fluxes.txt"
/* or, include here your own file */

  LB(j)     Lower bound on fluxes
$include "lowerbound_on_fluxes.txt"
/* or, include here your own file */

  S(i,j)    Stoichiometric matrix
$include "S_matrix.txt"
* in form of sparse matrix
/* or, include here your own file */

  epsilon   A small number
;

epsilon = 0.001;


*********************** Define the variables *********************
VARIABLES
	v(j)    fluxes
	z       Objective function
;

v.lo(j)=LB(j);
v.up(j)=UB(j);

biNARY VARIABLE
	xp(i)       1 if metabolite i is produced and zero otherwise
	w(i,j)      1 if metabolite i is produced by reaction j and zero otherwise
;

EQUATIONS
   obj
   massbalance_cytosol(i)
   massbalance_other(i)
   prodconsirrev1
   prodconsirrev2
   prodconsrev1
   prodconsrev2
   binarycons
;

* Objective Function maximizing number of metabolites that can be produced in the network*
obj..                                                  z =e= sum(i,xp(i));

* Mass balance constraints for cytosolic  and non-cytosolic metabolites
massbalance_cytosol(i)$(cytosol(i) and (not extracellular(i))).. sum(j,S(i,j)*v(j)) =g= 0;
massbalance_other(i)$(not cytosol(i) and not extracellular(i)).. sum(j,S(i,j)*v(j)) =e= 0;

* Production constraints for irreversible and reversible reactions
prodconsirrev1(i,j)$([S(i,j) gt 0] and [not rev(j)])..  v(j) =g= epsilon*w(i,j);
prodconsirrev2(i,j)$([S(i,j) gt 0] and [not rev(j)])..  v(j) =l= 1000*w(i,j);
prodconsrev1(i,j)$([S(i,j) ne 0] and rev(j))..          S(i,j)*v(j) =g= epsilon-1000*(1-w(i,j));
prodconsrev2(i,j)$([S(i,j) ne 0] and rev(j))..          S(i,j)*v(j) =l= 1000*w(i,j);

binarycons(i)..    sum(j$(([S(i,j) ne 0] and rev(j)) or ([S(i,j) gt 0] and [not rev(j)])), w(i,j)) =g= xp(i);


******************** Definition of the model *********************
MODEL GAPFIND
/
	obj
	massbalance_cytosol
	massbalance_other
	prodconsirrev1
	prodconsirrev2
	prodconsrev1
	prodconsrev2
	binarycons
/;



******************** Solve the model and store the results  ******************
SOLVE GAPFIND USING MIP MAXIMIZING z;

FILE res /no_production_metabolites.txt/;
PUT res;

if(GAPFIND.modelstat=1,
  PUT "/"/;
   LOOP(i$(xp.l(i) = 0),
      PUT "'",i.tl:0,"'"/;
  );
  PUT "/"/;
else
  PUT "No optimal solution! model status = ",GAPFIND.modelstat:0:0/;
);
