###### flux analysis: functions involving running simple optimization problems ######


get.opt.flux <- function(model, rxns="biomass", coefs=1, dir="max", ko=NULL, xopt=FALSE, lp.pars=list()) {
  # get the max or min flux of a single reaction in the model, given in rxns as indices or IDs as in model$rxns; if rxns="biomass" or is a single string containing the substing "biomass" (the default), optimize the biomass reaction as found by get.biomass.idx using rxns as the regex
  # or, get the max or min flux of a linear combination of reactions, with the corresponding linear coefficients given in coefs
  # to knockout reaction(s), pass KO as their indices or IDs; the lb's and ub's of these reactions will be set to 0
  # if xopt=TRUE, the optimal xopt vector will also be returned, then will return a list of list(obj, xopt); otherwise just the optimal obj value
  # lp.pars: a list of control parameters passed to solver

  lp.pars <- get.pars("lp", lp.pars)
  c <- rep(0, ncol(model$S))
  if (length(rxns)==1 && grepl("biomass", rxns, ignore.case=TRUE)) {
    i <- get.biomass.idx(model, rgx=rxns)
  } else i <- all2idx(model, rxns)
  c[i] <- coefs
  lb <- model$lb
  ub <- model$ub
  ko <- all2idx(model, ko)
  lb[ko] <- 0
  ub[ko] <- 0

  res <- solve.model(model, csense=dir, c=c, lb=lb, ub=ub, pars=lp.pars)[[1]]
  if (res$stat!=1) {
    # if solver status not 1 (optimal), to be safe, set result to NA
    res$xopt <- NA
    if (!res$stat %in% c(2,12)) res$obj <- NA # exclude the case of solver status 2 (unbounded) and 12, where obj will already by +/-Inf
  }
  if (xopt) {
    return(list(obj=res$obj, xopt=res$xopt))
  } else return(res$obj)
}

fba <- get.opt.flux # flux balance analysis as a synonym of get.opt.flux

get.opt.fluxes <- function(model, rxns="all", dir="max", nc=1L, lp.pars=list()) {
  # a wrapper around get.opt.flux for getting the min or max fluxes of multiple single reactions at once; return a named vector
  # thus this function is not for optimizing linear combination of reactions; and the direction will be either min or max for all rxns; no xopt is returned

  lp.pars <- get.pars("lp", lp.pars)
  if (length(rxns)==1 && rxns=="all") {
    rxns <- 1:length(model$rxns) # I use model$rxns instead of ncol(S) since S can contain extra columns
  } else rxns <- all2idx(model, rxns)
  res <- unlist(parallel::mclapply(rxns, get.opt.flux, model=model, dir=dir, lp.pars=lp.pars, mc.cores=nc))
  names(res) <- model$rxns[rxns]
  res
}

fva <- function(model, rxns="all", nc=1L, max.biomass=FALSE, biomass.rgx="biomass", lp.pars=list()) {
  # flux variability analysis, for rxns given as indices of IDs as in model$rxns; "all" for all rxns
  # max.biomass: whether to fix biomass at the maximum; biomass.rgx: regex used to find the biomass reaction
  # return a named matrix with two columns "min" and "max", rxns in the rows

  lp.pars <- get.pars("lp", lp.pars)
  if (max.biomass) {
    bm.idx <- get.biomass.idx(model, biomass.rgx)
    model <- set.rxn.bounds(model, bm.idx, lbs=1, ubs=1, relative=TRUE, nc=1L, lp.pars)
  }
  min <- get.opt.fluxes(model, rxns, "min", nc, lp.pars)
  max <- get.opt.fluxes(model, rxns, "max", nc, lp.pars)
  cbind(min, max)
}

run.ko.screen <- function(model, rxns="all+ctrl", f, ..., nc=1L, simplify=TRUE) {
  # a wrapper function to run algorithm `f` under the original wildtype model as well as models where each rxn in rxns is knocked out (i.e. both lb and ub set to 0)
  # model is the formulated model for running the algorithm `f`; it should be derived from the original metabolic model by adding columns and rows to the S matrix, such that indices of the rxns and mets do not change
  # rxns are the indices or IDs or rxns to include in the KO screen, by default all rxns plus the control wild-type model; rxns="all" to run for all rxns w/o the ctrl; if using indices, 0 means ctrl; if using IDs, "ctrl", means ctrl
  # ... are further parameters passed to f()
  # nc: the number of cores for running across rxns
  # simplify: if FALSE, return a list from mclapply across rxns; if TRUE, then assume the run for each rxn yields a summary data.table, and will try to combine them with rbindlist

  if (rxns=="all+ctrl") {
    rxns <- 0:length(model$rxns) # I use model$rxns instead of ncol(S) since S can contain extra columns
    names(rxns) <- c("ctrl", model$rxns)
  } else if (rxns=="all") {
    rxns <- 1:length(model$rxns)
    names(rxns) <- model$rxns
  } else if (rxns=="ctrl") {
    rxns <- 0
    names(rxns) <- "ctrl"
  } else if (is.character(rxns)) {
    tmp <- rxns
    rxns <- match(rxns, model$rxns)
    names(rxns) <- tmp
    if ("ctrl" %in% names(rxns)) rxns["ctrl"] <- 0
  } else if (is.numeric(rxns)) {
    names(rxns)[rxns!=0] <- model$rxns[rxns[rxns!=0]]
    names(rxns)[rxns==0] <- "ctrl"
  }

  res <- parallel::mclapply(rxns, function(i) {
    m <- model
    m$lb[i] <- 0 # if i==0, lb will not be changed (i.e. the control)
    m$ub[i] <- 0 # if i==0, ub will not be changed (i.e. the control)
    f(model, ...)
  }, mc.cores=nc)

  if (simplify) {
    if (any(!sapply(res, is.data.table))) {
      warning("Outputs contain non-data.table objects, cannot simplify.")
    } else res <- cbind(data.table(id=rxns), rbindlist(res, idcol="rxn"))
  }
  res
}


