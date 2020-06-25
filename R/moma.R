###### minimization of metabolic adjustment (MOMA) algorithm ######


minimize.dflux <- function(model, flux0, solv.pars=list()) {
  # solve the QP problem of min ||(flux-flux0)||_2
  solv.pars <- get.pars("qp", solv.pars)

  solve.model(model, csense="min", c=-2*flux0, Q=.sparseDiagonal(ncol(model$S), x=2), pars=solv.pars) # need to use .sparseDiagonal() instead of Diagonal() since Rcplex requires Q to be of the class dsparseMatrix
}

moma <- function(model, rxns="all", nc=1L, flux0="biomass", obj="biomass", biomass.rgx="biomass", xopt=FALSE, solv.pars=list()) {
  # run moma, the original formulation: first obtain a reference flux0 corresponding to the maximum biomass in the wildtype model (or provide a specific flux0); then minimize the euclidian distance to flux0 in the KO model and obtain the resulting biomass (or other objective value)
  # flux0: should either be "biomass" (then take the first solution of maximizing biomass as flux0), or "max.c" or "min.c" (then will maximize/minimize model$c to obtain flux0), or a vector for the actual flux0
  # obj: should either be "biomass" (then will get the biomass value as the reault) or "c" (then will get the model$c value) or a function (then will use this function on the output of solve.model() to get the result; the return of this function should be a 1-row data.table)
  # rxns are the indices or IDs or rxns to run moma on, by default all rxns; or a list, each element containing multiple reactions to KO at the same time
  # nc: number of cores
  # biomass.rgx is the regex used to find the biomass rxn
  # xopt: whether to include xopt in the result

  if (length(flux0)==1 && is.character(flux0)) {
    if (flux0=="biomass") {
      bm.idx <- get.biomass.idx(model, biomass.rgx)
      flux0 <- get.opt.flux(model, rxns=bm.idx, xopt=TRUE)$xopt
    } else if (flux0=="max.c") {
      flux0 <- get.opt.flux(model, 1:ncol(model$S), model$c, xopt=TRUE)$xopt
    } else if (flux0=="min.c") {
      flux0 <- get.opt.flux(model, 1:ncol(model$S), model$c, dir="min", xopt=TRUE)$xopt
    } else stop("Invalid value for the flux0 argument.")
  } else if (!is.numeric(flux0)) {
    stop("Invalid value for the flux0 argument.")
  }

  if (length(obj)==1 && is.character(obj)) {
    if (obj=="biomass") {
      bm.idx <- get.biomass.idx(model, biomass.rgx)
      obj <- function(x) data.table(obj=x[[1]]$xopt[bm.idx])
    } else if (obj=="c") obj <- function(x) data.table(obj=sum(x[[1]]$xopt*model$c))
  } else if (!is.function(obj)) stop("Invalid value for the obj argument.")
  
  moma0 <- function(model, flux0, obj, pars, xopt) {
    tmp <- minimize.dflux(model, flux0, pars)
    if (xopt) cbind(obj(tmp), data.table(xopt=list(tmp[[1]]$xopt))) else obj(tmp)
  }

  run.ko.screen(model, rxns, moma0, flux0=flux0, obj=obj, pars=solv.pars, xopt=xopt, nc=nc, simplify=TRUE)
}

moma2 <- function(model, rxns="all", nc=1L, obj0=c("biomass","max.c","min.c"), obj="biomass", biomass.rgx="biomass", solv.pars=list()) {
  # run moma (the alternative formulation: first obtain the maximum biomass or other objective in the wildtype model; then model the wildtype and KO as a single model, minimize the euclidian distance between flux_wt and flux_ko while fixing the biomass/other objective of the wildtype at the optimal value; finally obtain the resulting biomass/other objective of KO)
  # obj0: the objective for the wildtype model; either "default" (biomass) or "max.c"/"min.c" (maximize/minimize model$c)
  # obj: should either be "biomass" (then will get the biomass value as the reault) or "c" (then will get the model$c value) or a function (then will use this function on the output of solve.model() to get the result; the return of this function should be a 1-row data.table)
  # rxns are the indices or IDs or rxns to run moma on, by default all rxns; or a list, each element containing multiple reactions to KO at the same time
  # nc: number of cores
  # biomass.rgx: the regex used to find the biomass rxn

  bm.idx <- get.biomass.idx(model, biomass.rgx)
  model.wt <- model
  obj0 <- match.arg(obj0)
  if (obj0=="biomass") {
    obj0.opt.wt <- get.opt.flux(model, bm.idx)
    model.wt <- set.rxn.bounds(model.wt, bm.idx, obj0.opt.wt, obj0.opt.wt)
  } else {
    if (obj0=="max.c") obj0.opt.wt <- get.opt.flux(model, 1:ncol(model$S), model$c)
    if (obj0=="min.c") obj0.opt.wt <- get.opt.flux(model, 1:ncol(model$S), model$c, dir="min")
    model.wt <- add.constraint(model.wt, 1:ncol(model$S), model$c, obj0.opt.wt, obj0.opt.wt)
  }
  # create a combined model with the optimal-biomass-constrained wildtype model placed in the second place
  model.comb <- c.model(model, model.wt)
  model.comb$rxns[1:length(model$rxns)] <- stringr::str_sub(model.comb$rxns[1:length(model$rxns)], 1, -3)

  if (length(obj)==1 && is.character(obj)) {
    if (obj=="biomass") {
      obj <- function(x) data.table(obj=x[[1]]$xopt[bm.idx])
    } else if (obj=="c") {
      obj <- function(x) data.table(obj=sum(x[[1]]$xopt*model$c))
    } else stop("Invalid value for the obj argument.")
  } else if (!is.function(obj)) stop("Invalid value for the obj argument.")
  moma0 <- function(model, obj, pars) {
    pars <- get.pars("qp", pars)
    c <- rep(0, ncol(model$S))
    tmp <- .sparseDiagonal(ncol(model$S)/2, x=2) # need to use .sparseDiagonal() instead of Diagonal() since Rcplex requires Q to be of the class dsparseMatrix
    Q <- rbind(cbind(tmp, -1*tmp), cbind(-1*tmp, tmp)) # needs -1*tmp, -tmp doesn't work for the class somehow
    obj(solve.model(model, csense="min", c=c, Q=Q, pars=pars))
  }

  # need to parse rxns here so that they are for the KO model only (run.ko.screen will do it across all rxns in the combined model)
  if (length(rxns)==1 && rxns=="all") rxns <- 1:length(model$rxns) else rxns <- all2idx(model, rxns)
  run.ko.screen(model.comb, rxns, moma0, obj=obj, pars=solv.pars, nc=nc, simplify=TRUE)
}



