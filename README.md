# gembox

In-house toolbox for genome-scale metabolic modeling (GEM) for the Ruppin lab. This is an R package that implements GEM functionalities.

Provided as is without warranty of any kind.

### Current Functions

* Query of metabolic model/network data, and basic manipulations of metabolic models
* FBA, FVA, FCA
* Simulating reaction knockouts; MOMA
* Sampling: ACHR
* Metabolic-EP
* Incorporating gene expression data: iMAT, PRIME, GIMME, E-Fmin, the method of Lee and Smallbone, and their variations
* MTA, rMTA and their variations
* Differential flux analysis and metabolic pathway enrichment analysis
* Visualization of metabolic network and differential flux results

### Dependencies

* ILOG CPLEX Optimization Studio or Gurobi (free academic licenses available)
* R packages
  - Depends: Matrix, data.table  
  - Imports: stringr, parallel, pbmcapply  
  - LinkingTo: Rcpp, RcppArmadillo  
  - Suggests: Rcplex2 (ruppinlab/Rcplex2), gurobi, R.matlab, sybilSBML, rlist, fgsea, igraph, ggplot2, RColorBrewer, visNetwork, hypergraph, hyperdraw
  - At least one of Rcplex2 (ruppinlab/Rcplex2) and gurobi is required for running optimizations

This package for now only works on Linux and MacOS.
