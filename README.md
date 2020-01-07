# gembox

In-house toolbox for genome-scale metabolic modeling (GEM) for the Ruppin lab.

This is an R package that implements basic functionality of GEM, as well as several advanced GEM algorithms by our lab including iMAT, MTA and their variations.

### Dependencies

* ILOG CPLEX Optimization Studio (free academic version available from IBM official website)
* R packages
  - Depends: Matrix, data.table  
  - Imports: Rcplex2 (ruppinlab/Rcplex2), stringr, parallel, igraph  
  - LinkingTo: Rcpp, RcppArmadillo  
  - Suggests: R.matlab, sybilSBML, rlist, fgsea, hypergraph, hyperdraw, RColorBrewer, visNetwork

This package has so far only been tested on Linux systems.
