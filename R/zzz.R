.onLoad <- function(libname, pkgname) {

  if (requireNamespace("Rcplex2", quietly=TRUE)) {
    solver("rcplex")
    message("Solver set to rcplex.")
  } else if (requireNamespace("gurobi", quietly=TRUE)) {
    solver("gurobi")
    message("Solver set to gurobi.")
  } else {
    warning("No solver found! At least one of the packages \"Rcplex2\", \"gurobi\" needed for running optimization.")
  }
  
} 
