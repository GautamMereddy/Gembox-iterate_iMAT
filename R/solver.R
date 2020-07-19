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

    tmpf <- function(x) {
      res <- list(stat=.pkg.const$cpx.stat.code[as.character(x$status)], obj=x$obj, xopt=x$xopt)
      if (x$status %in% c(2,118)) {
        # unbounded (2,118) problem
        if (objsense=="min") res$obj <- -Inf
        if (objsense=="max") res$obj <- Inf
      }
      res
    }

    res <- tryCatch({
      res <- Rcplex2::Rcplex(c, A, b, Q, lb, ub, x0, pars, csense, sense, vtype, nsol)
      if (!is.null(names(res))) res <- list(res)
      res <- lapply(res, tmpf)
      if (!res[[1]]$stat %in% .pkg.const$ok.stat) warning("Potential issue, solver status: ", res[[1]]$stat, ".")
      res
    }, error=function(e) {
      warning("In solve.model(): Error while trying to solve model (see info below), NA returned.\n", e)
      res <- list(list(stat=NA, obj=NA, xopt=NA))
    })

  } else if (.pkg.var$solver=="gurobi") {

    m <- list(obj=c, A=A, rhs=b, Q=Q, lb=lb, ub=ub, modelsense=csense, sense=ifelse(sense=="G",">","<"), vtype=vtype)
    if (!is.null(x0)) {
      if (any(vtype!="C")) m$start <- x0 else m$pstart <- x0
    }

    res <- tryCatch({
      res <- gurobi::gurobi(m, pars)
      if (!res$status %in% .pkg.const$ok.stat) warning("Potential issue, solver status: ", res$status, ".")
      if ("pool" %in% names(res)) {
        res <- lapply(res$pool, function(x) list(stat=res$status, obj=x$objval, xopt=x$x))
      } else res <- list(list(stat=res$status, obj=res$objval, xopt=res$x))
      res
    }, error=function(e) {
      warning("In solve.model(): Error while trying to solve model (see info below), NA returned.\n", e)
      res <- list(list(stat=NA, obj=NA, xopt=NA))
    })

  }

  res
}
