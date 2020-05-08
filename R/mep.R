###### the Metabolic-EP algorithm (Braunstein et al. Nat Commun 2017) ######


### --- prepare data for MEP --- ###

rm.fixed.rxns <- function(model, nc=1L) {
  fva.res <- fva(model, nc=nc)
  vmin <- pmax(model$lb, fva.res$vmin)
  vmax <- pmin(model$ub, fva.res$vmax)
  fixed <- vmin==vmax
  notfixed <- !fixed
  if (any(fixed)) {
    S <- model$S
    model$S = S[,notfixed]
    model$b = model$b - as.vector(S[,fixed,drop=FALSE] %*% cbind(vmin[fixed]))
    model$lb = vmin[notfixed]
    model$ub = vmax[notfixed]
    model$rxns = model$rxns[notfixed]
    model$rxnNames = model$rxnNames[notfixed]
  }
  v <- vmin
  v[notfixed] <- NA
  list(model=model, v=v, idx=which(notfixed))
}


### --- MEP --- ###

mep <- function(model, nc=1L, fix.flux=FALSE, fflux.id=0, fflux.mean=0, fflux.var=0, mep.pars=list()) {
  # the main function for running MEP
  # nc: number of cores to use for FVA in pre-processing model
  # fix.flux: whether to incorporate experimentally measured fluxes, if so, pass in the relevant info via fflux.id, fflux.mean and fflux.var
  # return: list(model, result.mep)
  pars <- get.pars("mep", mep.pars)

  # pre-process model by removing fixed reactions (i.e. lb==ub via FVA)
  message("Pre-processing model...")
  model1 <- rm.fixed.rxns(model, nc=nc)
  message("Done pre-processing model.")

  # run MEP
  mep.res <- mepc(model1$model, pars$beta, pars$damp, pars$max.iter, pars$dlb, pars$dub, pars$epsil, fix.flux, fflux.id, fflux.mean, fflux.var)

  # post-processing: incorporate MEP results into the original metabolic model
  tmpf <- function(x, a=0) {
    if (a==0) res <- model1$v else if (a==1) res <- rep(0, length(model1$v))
    res[model1$idx] <- mep.res[[x]]
    mep.res[[x]] <<- res
  }
  for (i in c("means.post","means.prior","means.trunc")) tmpf(i)
  for (i in c("vars.post","vars.prior","vars.trunc")) tmpf(i, 1)
  mep.res$cov.idx <- model1$idx

  list(model=model, result.mep=mep.res)
}


