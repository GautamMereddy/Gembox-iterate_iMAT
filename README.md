# Gembox - Iterative iMAT Pipeline

Jason H. Yang Lab @ Rutgers New Jersey Medical School

Author: Gautam Mereddy

Date: March 5, 2024



What is Gembox?
===============

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
  - For iterate_imat: tools, reticulate

This package for now only works on Linux and MacOS.

What is iterate iMAT?
=====================

In Gembox's R folder, an additional script called iterate_imat.R has been added which includes additional functions that allow for the automatic, simultaenous creation of multiple context-specific genome-scale metabolic models from multiple transcriptomic datasets using Gembox's iMAT algorithm. This is done through a single function, iterate_imat(). The iterate_imat() function uses parallel processing to simultaneously apply the iMAT gene integration algorithm to multiple samples in a gene by sample expression matrix. The function applies the same iMAT parameters (upper/lower bounds) to each gene expression matrix sample and returns a context-specific metabolic model for each sample (in Matlab's .mat format). This function was designed to be used in conjunction with our other Python-based, in-house pipeline iterate_sampling (based off the opencobra/cobrapy repository). The iterate_imat() function is capable of automatically flux sampling each context-specific iMAT model made through COBRApy's OptGpSampler using the Reticulate package. 

All functions in iterate_imat.R have been tested with the following specifications:

- R 3.6.0
- Gembox 0.1.0
- Rcplex2 0.3-3
- CPLEX 12.6.3
- Python 3.7.11
- COBRApy 0.21.0
- Gurobi 9.5.2

### Using iterate iMAT

After installing Gembox and it's dependencies, the functions in iterate_imat.R can be called by first sourcing the iterate_imat R file. To be able to automatically flux sample the iMAT-derived metabolic models, a Python environment with COBRApy must be speciifed in Reticulate first. 



