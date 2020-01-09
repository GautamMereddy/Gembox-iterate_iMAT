###### flux coupling analysis ######


prep.model.for.fca <- function(model, nc=1L, lp.pars=list()) {
  # prepare metabolic model for flux coupling analysis: first convert.rev.rxns() then rm.blocked.rxns()

  lp.pars <- get.pars("lp", lp.pars)
  model <- convert.rev.rxns(model)
  rm.blocked.rxns(model, nc=nc, solv.pars=lp.pars)
}

fca <- function(model, rxn0, rxn1, prep.model=FALSE, nc=1L, lp.pars=list()) {
  # flux coupling analysis for a single pair of reactions; the formalization is rxn1/rxn0 (i.e. normalized by rxn0)
  # prep.model: whether to run prep.model.for.fca(); by default the model should already be properly prepared

  lp.pars <- get.pars("lp", lp.pars)
  if (prep.model) model <- prep.model.for.fca(model, nc, lp.pars)
  # formalize fca model
  model$lb <- rep(0, length(model$lb))
  ub <- model$ub
  model$ub <- rep(1e20, length(model$ub))
  model$lb[rxn0] <- 1
  model$ub[rxn0] <- 1
  trxns <- which(apply(model$S, 2, function(x) all(x>=0) || all(x<=0))) # transport reactions
  n.trxns <- length(trxns)
  nr <- nrow(model$S)
  nc <- ncol(model$S)
  model$S <- rbind(cbind(model$S,                                             rep(0, nr)),
                   cbind(sparseMatrix(1:n.trxns, trxns, dims=c(n.trxns, nc)), -ub[trxns]))
  model$lb <- c(model$lb, 0)
  model$ub <- c(model$ub, 1e20)
  model$rowlb <- c(model$rowlb, rep(-1e20, n.trxns))
  model$rowub <- c(model$rowub, rep(0, n.trxns))

  # run
  rmax <- get.opt.flux(model, rxn1, solv.pars=lp.pars)
  rmin <- get.opt.flux(model, rxn1, dir="min", solv.pars=lp.pars)
  c(rmin=rmin, rmax=rmax)
}


