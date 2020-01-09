###### functions for optimization solver APIs ######


solve.model <- function(model, ..., x0=NULL, pars=list()) {
  # solve the optimization problem as coded in model, unless explicitly overriding model parameters by passing arguments to "..."
  # arguments that can be overrided by passing to "..." are listed below (otherwise will throw an error):
  valid.args <- c("csense","c","lb","ub","Q")
  # x0: initial solution (warm start)
  # pars: a list of solver parameters
  pars <- get.pars("mip", pars) # by default use the "mip" default parameters

  args <- list(...)
  if (any(!names(args) %in% valid.args)) stop("Arguments passed to ... invalid.")
  if ("csense" %in% names(args)) csense <- args$csense else csense <- model$csense
  if ("c" %in% names(args)) c <- args$c else c <- model$c
  if ("Q" %in% names(args)) Q <- args$Q else Q <- model$Q
  A <- rbind(model$S, model$S)
  b <- c(model$rowlb, model$rowub)
  sense <- rep(c("G","L"), c(length(model$rowlb), length(model$rowub)))
  if ("lb" %in% names(args)) lb <- args$lb else lb <- model$lb
  if ("ub" %in% names(args)) ub <- args$ub else ub <- model$ub
  vtype <- model$vtype

  if (.pkg.var$solver=="rcplex") {
    if ("nsol" %in% names(pars)) {
      nsol <- pars$nsol
      pars$nsol <- NULL
    } else nsol <- 1
    res <- Rcplex2::Rcplex(c, A, b, Q, lb, ub, x0, pars, csense, sense, vtype, nsol)
    tmpf <- function(x) {
      res <- x[c("xopt", "obj")]
      res$stat <- x$status
      res$stat.str <- .pkg.const$cpx.stat.code[as.character(x$status)]
      if (res$stat %in% c(2,118,12)) {
        # unbounded (2,118) or possibly unbounded (12, maybe others but I'm not sure) problem
        if (objsense=="min") res$obj <- -Inf
        if (objsense=="max") res$obj <- Inf
      } else if (!res$stat %in% c(1,101,102,128,129,130)) {
        warning("In solve.model(): Potential issue, solver status: ", res$stat.str, ", returning NA for xopt and objective value.", call.=FALSE)
        res$obj <- NA
        res$xopt <- NA
      }
      res
    }
    if (!is.null(names(res))) res <- list(res)
    res <- lapply(res, tmpf)
  }
  # todo: add cplexAPI and gurobi
  res
}
