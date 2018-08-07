*****************************************************************************
*                               GAPFILL                                     *
*                               -------                                     *
* Related publication: PMID:17584497                                        *
*    - This code identifies determines the minimal number of modifications  *
*      required to restore connectivity of each problem metabolite          * 
*      in the model                                                         *
*    - A sample list of problem metabolites is given in the file            *
*      no_production_metabolites.txt, which is imported into GAMS as a      *
*      set called problme_met(i)                                            *
*    - The code was written for the iMM904 model of yeast as an example     *
*    - The code can be easily modified for any other metabolic model        *
*    - The results of this code are stored in a file called results.txt. A  *
*      sample of results from this code is available in sample_results.txt  *
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
       domlim = 1000
       limcol = 1000
       limrow = 1000
       optca = 0.0
       optcr = 0.0
       work = 10000000
       nlp = minos5
       mip = cplex; 
       
SETS
	i     Name of the metabolites in the model
$include "compounds.txt"
/* or, include your own file here */

	j     Name of the rxns in the original model and those in the set database
$include "reactions.txt"
/* or, include your own file here */

        model_rxn(j)  Reactions in the original model
$include "model_reactions.txt"
/* or, include your own file here */

	database(j)  Customozed database of reactions
* This database contains the transport rxns between compartments and cytosol and 
* between cytosl and extracellular environment as well as the exchange rxns that are
* mmissing from the original mode as well as the exchange rxns that are
* mmissing from the original model. The equation for the reactions in this database 
* is given in the file database.txt 
$include "database_reactions.txt"
/* or, include your own file here */
	
	rev(j)   Reversible reactions (including exchange rxns)	
$include "reversible_rxns.txt"
/* or, include your own file here */

	problem_met(i)   Problem metabolites
$include "no_production_metabolites.txt"
/* or, include your own file here */
	
	cytosolic(i)    Cytosolic metabolites
$include "cytosolic_metabolites.txt"
/* or, include your own file here */

	extracellular(i)  Extracellular metabolites
$include "extracellular_metabolites.txt"
/* or, include your own file here */

       current_metab(i)  A dynamic set containing the problem metaoblite is solved for

       iter /1*10000/

       eqncounter(iter) Dynamic set containing the current number of runs of the loop
;

* At first there is nothing in the set eqncounter
eqncounter(iter)=no;

********************* Define parameters *******************
PARAMETERS
  UB(j)  Upper bound on fluxes
$include "upperbounds_on_fluxes.txt"
/* or, include your own file here */

  LB(j)  Lower bound on fluxes
$include "lowerbounds_on_fluxes.txt"
/* or, include your own file here */

  S(i,j)   Stoichiomteric matrix 
$include "S_matrix.txt"

  prev_y(j,iter) Stores the optimal value of binary variables at each iteration

  epsilon  A small value

  n_min  The minimum number of added rxns

  counter

  done  Represents when we are done with finding the rxns to resol
;

epsilon = 0.001;

********************* Define variables *******************
VARIABLES
	v(j)      Fluxes
	Z         objective	
;

biNARY VARIABLE
	y(j)          1 if reaction j from the database is added and zero otherwise	
	w(i,j)        1 if reaction j (either database or model) produces metabolite i and zero otherwise
;

********************* Define equations *******************
EQUATIONS
	obj
	massbalance_cytosol(i)
	massbalance_other(i)
	boundcon1(j)
	boundcon2(j)
	boundcon3(j)
	boundcon4(database)
	boundcon5(database)
	prodconst1(i,j)
	prodconst2(i,j)
	binarycons(i)

* Equations for deleting the previous solutions
       integercut(iter)  This eqn is defined over a dynamic set containing the current number of iterations
;


* Minimize the number of reactions added from the database
obj ..              z =e= sum(j, y(j));

* Mass balance constraints on cytosolic and non cytosolic metabolites
massbalance_cytosol(i)$(cytosolic(i)) ..   sum(j$(S(i,j)),S(i,j)*v(j)) =g= 0;
massbalance_other(i)$(not cytosolic(i))..  sum(j$(S(i,j)),S(i,j)*v(j)) =e= 0;

* Production constraints that ensure production of problem metabolite
prodconst1(i,j)$(current_metab(i) and [S(i,j) ne 0]).. S(i,j)*v(j) =g= epsilon-1000*(1-w(i,j));
prodconst2(i,j)$(current_metab(i) and [S(i,j) ne 0]).. S(i,j)*v(j) =l= 1000*w(i,j);

* Bound constraint on reactions in the database
boundcon1(j)$(model_rxn(j))..                   v(j) =g= LB(j);
boundcon2(j)$(model_rxn(j) and rev(j))..        v(j) =l= UB(j);
boundcon3(j)$(model_rxn(j) and [not rev(j)])..  v(j) =g= -1000*y(j);
boundcon4(database)..                           v(database) =g= LB(database)*y(database);
boundcon5(database)..                           v(database) =l= UB(database)*y(database);

binarycons(i)$(current_metab(i))..   sum(j$(s(i,j) ne 0), w(i,j)) =g= 1;

* Integer cuts to preclude the previously found solutions
integercut(eqncounter)..   sum(j$(prev_y(j,eqncounter) eq 1),y(j)) =l= sum(j,prev_y(j,eqncounter)) - 1;


********** Define and solve the model. Store the results in output files  *************
model GAPFILL 
/
        obj
        massbalance_cytosol
        massbalance_other
        boundcon1
        boundcon2
        boundcon3
        boundcon4
        boundcon5
        prodconst1
        prodconst2
        binarycons
        integercut
/
;

* Include th option file
GAPFILL.optfile=1;

* Initializing prev_y for the first run of the loop
prev_y(j,iter)=no;
n_min=0;
done=0;
counter=0;

alias(i,i1);

FILE resultfile /results.txt/;
resultfile.pw=500;
PUT resultfile;

LOOP(i1$(problem_met(i1)),
   current_metab(i)=no;
   current_metab(i1)=yes;   /* current_metab contains only i1 */

   eqncounter(iter)=no;

   SOLVE GAPFILL USING MIP MINIMIZING z;

   PUT /"*********** ",i1.tl:0," *************"/;
   if(GAPFILL.modelstat = 10,       /*If the model is infeasible */
     PUT resultfile;
     PUT " Reactions in the set Dataset cannot fix this problem metabolite! "/;
   elseif (GAPFILL.modelstat ne 1),  /* If a non-optimum solution is found */
     PUT resultfile;
     PUT "ERROR: modelstat= "GAPFILL.modelstat/; /* Refer to GAMS manual for the modelstat*/
   elseif (z.l=0),   /* If it turns out that no rxns is needed */
     PUT resultfile;
     PUT "Objective = 0, No rxns need to be added! This is not a problem metabolite!"/;
   else     /* If at least one optimal solution exists */
      counter=0;
      n_min=z.l;   /* store the minimal number of rxns needed */
      PUT "This problem metabolite can be solved by adding a minimum of ",n_min:0:1," rxns to the model"//; 
      done=0;
      LOOP(iter$(not done),  /* Now find all alternative solutions */
         counter=counter+1;
         SOLVE GAPFILL using mip minimizing z;
         if(GAPFILL.modelstat eq 10,
            done=1;
            PUT resultfile;
            PUT "ERROR: modelstat= "GAPFILL.modelstat:0," (Integer infeasible): No more solutions with the reactions in this database!"/;
         elseif (GAPFILL.modelstat ne 1),
            done=1;
            PUT resultfile;
            PUT "ERROR: modelstat= "GAPFILL.modelstat:0,"  No optimal solution!"/;
         elseif (z.l > n_min+1),   /* If number of rxns needed is greater than n_min+1 */ 
            done=1;           /* Do not continue */ 
            PUT "Ended with obj = ",z.l:0:1/;
         else                /* Continue finding other alternative solutions */
            eqncounter(iter)=yes;  /* Add one more element to our dynamic set */
            prev_y(j,iter)=0;  /* Initiate prev_y for this iteration */ 
            PUT resultfile;         
            PUT counter:0:0,") ";
            /* Note that if a rxn written in the output is a rxn of the origianl model */
            /* then it is an irreversible rxn, i.e., it means that the irreversability */
            /* on the model rxn need to be relaxed*/ 
            LOOP(j$(y.l(j) = 1),
                PUT " '",j.tl:0,"'";
                prev_y(j,iter)=1;  /* Store the soln at this iteration in prev_y */
            );
            PUT /;
         );
         if(counter = 10,      /* Find a maximum of 10 alternative solutions */
             done=1;
         );
      );
   );

);
