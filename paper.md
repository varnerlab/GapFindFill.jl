---
title: 'GFFJ.jl: optimization-based gap finding and filling in Julia'
tags:
  - Julia
  - metabolic network 
  - optimization
authors:
  - name: Zhiping Zhang
    orcid: 0000-0002-9099-4357
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Jeffrey D. Varner
    affiliation: 1
affiliations:
 - name: Robert Frederick Smith School of Chemical and Biomolecular Engineering, Cornell University, Ithaca NY, 14853 USA
   index: 1
 - name: Institution 2
   index: 2
date: 04 October 2019
bibliography: paper.bibtex

# # Optional fields if submitting to a AAS journal too, see this blog post:
# # https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

``describing the high-level functionality and purpose of the software for a diverse, non-specialist audience.`` 

**Abstract**
Metabolic reconstructions are widely used to study metabolism of
biological systems. Various computational methods have been proposed to
automatically curate large-scale metabolic networks. The application of
network curration tools is currently hampered by the scalability of
algorithms and the availability of software supports. `GFFJ.jl` is a
high-level and open-source implementation of optimization-based gap
finding and filling algorithm in Julia. It has two advantages over the
current implementation in `GAMS`. On one hand, it harnesses the power of
free academic license solvers, such as Gurobi, to solve large-scale
metabolic network curation problems as mixed integer linear programming
problems. On the other hand, it is built upon high-performance general
purpose programming language, Julia, thus, endows users the capability
of embedding network curation optimization into other computational
tasks. The code is freely available on
<https://github.com/varnerlab/GFFJ.git>.

**Introduction**
Metabolic reconstructions of different organisms from experimental evidence and bioinformatics based knowledge are widely used to facilitate the study of biological systems [@feist2009reconstruction; @schellenberger2011quantitative; @thiele2010protocol; @henry2010high]. 
All of these reconstructions are inherently incomplete due to our lack of complete experimental and/or homology information \cite{kumar2007optimization}. 
A number of computational approaches have been proposed to detect gaps in a metabolic network and subsequently generate hypotheses to fix these gaps \cite{vlassis2014fast, becker2008context, jerby2010computational, agren2012reconstruction, wang2012reconstruction, zur2010imat}.  
The optimization based \textit{GapFind} and \textit{GapFill} approach proposed by Maranas and coworkers is one of the most widely used \cite{kumar2007optimization, maranas2016optimization}.

The \textit{GapFind} identifies all no-production metabolites, by solving the following mixed integer linear programming problem \cite{maranas2016optimization}: 


where $I$ and $J$ are the set of compounds and reactions in the network, respectively,
$x_i$ is $1$ if compound $i$ can be produced in the network, otherwise 0,
$\epsilon$ denotes a minimum threshold for a reaction to be treated as active,
$S_{ij}$ denotes the stoichiometric coefficient for species $i$ in reaction $j$,
$v_j$ denotes the flux through reaction $j$,
$w_{ij}$ is $1$ if reaction $j$ is producing compound $i$ actively, 0 otherwise.  
$J^{ir}$ denotes the set of irreversible reactions. 
$LB_j$ and $UB_j$ are lower and upper bounds on flux $j$, respectively. 
$I^{cyt}$ denotes the set of cytosolic compounds.  

The \textit{GapFill} tries to propose ways of bridging each gap independently by solving a new mixed integer linear programming problem repeatedly \cite{kumar2007optimization}: 

# Statement of Need 

`` illustrates the research purpose of the software.``

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$


# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this: ![Example figure.](figure.png)

# Acknowledgements

We acknowledged the financial support to
JV from the departmental funding of Robert Frederick Smith
School of Chemical and Biomolecular Engineering, Cornell University.

# References
