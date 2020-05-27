# gembox

In-house toolbox for genome-scale metabolic modeling (GEM) for the Ruppin lab.

This is an R package that implements basic functionality of GEM, as well as several advanced GEM algorithms by our lab including iMAT, MTA and their variations.

### Dependencies

* ILOG CPLEX Optimization Studio or Gurobi (free academic licenses available)
* R packages
  - Depends: Matrix, data.table  
  - Imports: stringr, parallel  
  - LinkingTo: Rcpp, RcppArmadillo  
  - Suggests: Rcplex2 (ruppinlab/Rcplex2), gurobi, R.matlab, sybilSBML, rlist, fgsea, igraph, hypergraph, hyperdraw, RColorBrewer, visNetwork
  - At least one of Rcplex2 (ruppinlab/Rcplex2) and gurobi is required for running optimizations

This package for now only works on Linux and MacOS.
