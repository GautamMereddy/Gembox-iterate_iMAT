###### global variables used by package functions ######


### --- global variables --- ###

# an environment object used to store the global variables, such that the variables can be changed by user if needed
.pkg.var <- new.env(parent=emptyenv())
# current solver; todo: set default solver based on solver availability
.pkg.var$solver  <- "rcplex"


### --- global constants --- ###

.pkg.const <- list()

## various default parameters
# solver parameters (dependent on the solver used)
.pkg.const$lp <- list()
.pkg.const$lp$rcplex <- list(trace=0, maxcalls=2e4, tilim=120, threads=1, parallel.mode=1, nsol=1)
.pkg.const$qp <- list()
.pkg.const$qp$rcplex <- list(trace=0, maxcalls=2e4, tilim=120, threads=1, parallel.mode=1, nsol=1)
.pkg.const$mip <- list()
.pkg.const$mip$rcplex <- list(trace=1, maxcalls=2e4, tilim=1200, threads=1, parallel.mode=1, nodesel=0, epagap=1e-6, epgap=1e-4, nsol=1, solnpoolagap=0, solnpoolgap=0, solnpoolintensity=2)
# sampling parameters
.pkg.const$samp <- list(method="achr", n.sample=5e3, steps.per.pnt=400, n.warmup=5000, nc=1L)
# parameters for various algorithms
.pkg.const$imat <- list(mode=0, nc=1L, nsteps=1, sol=0, sol.major.cutoff=0.667, flux.act=1, flux.inact=0.1, flux.delta=0.1, flux.delta.rel=0, flux.bound=1000)
.pkg.const$mta <- list(v.min=-50, v.max=50, v.min.c=-1000, v.max.c=1000, alpha=0.9, epsil=0.01)

# cplex solver status codes
.pkg.const$cpx.stat.code = c(
  `1`="CPX_STAT_OPTIMAL",
  `2`="CPX_STAT_UNBOUNDED",
  `3`="CPX_STAT_INFEASIBLE",
  `4`="CPX_STAT_INForUNBD",
  `5`="CPX_STAT_OPTIMAL_INFEAS",
  `6`="CPX_STAT_NUM_BEST",
  `10`="CPX_STAT_ABORT_IT_LIM",
  `11`="CPX_STAT_ABORT_TIME_LIM",
  `12`="CPX_STAT_ABORT_OBJ_LIM",
  `13`="CPX_STAT_ABORT_USER",
  `14`="CPX_STAT_FEASIBLE_RELAXED_SUM",
  `15`="CPX_STAT_OPTIMAL_RELAXED_SUM",
  `16`="CPX_STAT_FEASIBLE_RELAXED_INF",
  `17`="CPX_STAT_OPTIMAL_RELAXED_INF",
  `18`="CPX_STAT_FEASIBLE_RELAXED_QUAD",
  `19`="CPX_STAT_OPTIMAL_RELAXED_QUAD",
  `20`="CPX_STAT_OPTIMAL_FACE_UNBOUNDED",
  `21`="CPX_STAT_ABORT_PRIM_OBJ_LIM",
  `22`="CPX_STAT_ABORT_DUAL_OBJ_LIM",
  `23`="CPX_STAT_FEASIBLE",
  `24`="CPX_STAT_FIRSTORDER",
  `25`="CPX_STAT_ABORT_DETTIME_LIM",
  `30`="CPX_STAT_CONFLICT_FEASIBLE",
  `31`="CPX_STAT_CONFLICT_MINIMAL",
  `32`="CPX_STAT_CONFLICT_ABORT_CONTRADICTION",
  `33`="CPX_STAT_CONFLICT_ABORT_TIME_LIM",
  `34`="CPX_STAT_CONFLICT_ABORT_IT_LIM",
  `35`="CPX_STAT_CONFLICT_ABORT_NODE_LIM",
  `36`="CPX_STAT_CONFLICT_ABORT_OBJ_LIM",
  `37`="CPX_STAT_CONFLICT_ABORT_MEM_LIM",
  `38`="CPX_STAT_CONFLICT_ABORT_USER",
  `39`="CPX_STAT_CONFLICT_ABORT_DETTIME_LIM",
  `101`="CPXMIP_OPTIMAL",
  `102`="CPXMIP_OPTIMAL_TOL",
  `103`="CPXMIP_INFEASIBLE",
  `104`="CPXMIP_SOL_LIM",
  `105`="CPXMIP_NODE_LIM_FEAS",
  `106`="CPXMIP_NODE_LIM_INFEAS",
  `107`="CPXMIP_TIME_LIM_FEAS",
  `108`="CPXMIP_TIME_LIM_INFEAS",
  `109`="CPXMIP_FAIL_FEAS",
  `110`="CPXMIP_FAIL_INFEAS",
  `111`="CPXMIP_MEM_LIM_FEAS",
  `112`="CPXMIP_MEM_LIM_INFEAS",
  `113`="CPXMIP_ABORT_FEAS",
  `114`="CPXMIP_ABORT_INFEAS",
  `115`="CPXMIP_OPTIMAL_INFEAS",
  `116`="CPXMIP_FAIL_FEAS_NO_TREE",
  `117`="CPXMIP_FAIL_INFEAS_NO_TREE",
  `118`="CPXMIP_UNBOUNDED",
  `119`="CPXMIP_INForUNBD",
  `120`="CPXMIP_FEASIBLE_RELAXED_SUM",
  `121`="CPXMIP_OPTIMAL_RELAXED_SUM",
  `122`="CPXMIP_FEASIBLE_RELAXED_INF",
  `123`="CPXMIP_OPTIMAL_RELAXED_INF",
  `124`="CPXMIP_FEASIBLE_RELAXED_QUAD",
  `125`="CPXMIP_OPTIMAL_RELAXED_QUAD",
  `126`="CPXMIP_ABORT_RELAXED",
  `127`="CPXMIP_FEASIBLE",
  `128`="CPXMIP_POPULATESOL_LIM",
  `129`="CPXMIP_OPTIMAL_POPULATED",
  `130`="CPXMIP_OPTIMAL_POPULATED_TOL",
  `131`="CPXMIP_DETTIME_LIM_FEAS",
  `132`="CPXMIP_DETTIME_LIM_INFEAS"
)


### --- showing/setting global variables --- ###

# only solver is meant to be changeable; all other global variables are meant to be constants

solver <- function(x=c("","rcplex")) {
  # show the current solver if no argument is passed, or set solver
  # for now, only Rcplex is available
  # todo: add cplexAPI and gurobi; check solver availability
  x <- match.arg(x)
  cur.solver <- .pkg.var$solver
  if (x=="") {
    message("Current solver is: ", cur.solver, ".")
  } else if (x!=cur.solver) {
    # set solver
    assign("solver", x, envir=.pkg.var)
    cur.solver <- x
  }
  invisible(cur.solver)
}

get.default.pars <- function(x=c("lp","qp","mip","samp","imat","mta"), print=TRUE) {
  # print and return the default values of parameters
  x <- match.arg(x)
  if (x %in% c("lp","qp","mip")) {
    res <- .pkg.const[[x]][[.pkg.var$solver]]
  } else res <- .pkg.const[[x]]
  if (print) print(res)
  invisible(res)
}

get.pars <- function(x=c("lp","qp","mip","samp","imat","mta"), pars) {
  # a helper function used in other functions to set parameters based on the default parameters and a list of pars
  x <- match.arg(x)
  default.pars <- get.default.pars(x, print=FALSE)
  new.pars <- c(default.pars[!names(default.pars) %in% names(pars)], pars)
}
