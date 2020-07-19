imat.dv0 <- function(model, expr=NULL, v0, df, w=NULL, imat.pars=list(), solv.pars=list(), samp.pars=NULL) {
  # the main function for running iMAT diff flux relative to a given v0
  # set samp.pars to NULL to skip the sampling
  # return both the imat.model, and the result model (the updated metabolic model), in a list(imat.model, result.model)

  # formulate iMAT model
  if (!is.null(expr)) imat.model <- form.imat(model, expr, imat.pars) else imat.model <- model
  imat.model <- form.imat.dv0(imat.model, v0, df, w, imat.pars)

  # solve the iMAT model
  imat.model <- run.imat.dv0(imat.model, imat.pars, solv.pars)

  # update the original metabolic model based on iMAT result
  res.model <- update.model.imat.dv0(model, imat.model, imat.pars)

  if (!is.null(samp.pars)) {
    # sample the result model
    res.model <- sample.model(res.model, samp.pars)
  }

  list(imat.model=imat.model, result.model=res.model)
}



form.imat.dv0 <- function(model, v0, df, w, imat.pars) {
  # formulate iMAT model for diff flux relative to v0
  # expr is the output from exprs2int()
  # imat.pars: the parameters for iMAT
  pars <- get.pars("imat", imat.pars)

  n.mets <- nrow(model$S)
  n.rxns <- ncol(model$S)
  S <- model$S

  # 1. df>0 (v>v0) reactions: specify the z+ indicator variables
  rxns.act <- which(df>0)
  n.act <- length(rxns.act)
  v.act <- (1+pars$flux.delta.rel)*abs(v0)[rxns.act]+pars$flux.delta
  if (n.act!=0) {
    m1 <- sparseMatrix(1:n.act, rxns.act, dims=c(n.act, n.rxns))
    m2 <- Diagonal(n.act, x=(-v.act-pars$flux.bound))
    S <- rbind(cbind(S, sparseMatrix(NULL, NULL, dims=c(n.mets, n.act))), cbind(m1, m2))
  }

  # 2. Reversible df>0 (|v|>|v0|) reactions: specify the additional z- indicator variables
  rxns.act.rev <- which(df>0 & model$lb[1:length(model$rxns)]<0)
  n.act.rev <- length(rxns.act.rev)
  #v.act <- (1+pars$flux.delta.rel)*abs(v0)[rxns.act.rev]+pars$flux.delta
  v.act <- pars$flux.delta.rel*abs(v0)[rxns.act.rev]+pars$flux.delta
  if (n.act.rev!=0) {
    m1 <- sparseMatrix(1:n.act.rev, rxns.act.rev, dims=c(n.act.rev, ncol(S)))
    m2 <- Diagonal(n.act.rev, x=v.act+pars$flux.bound)
    S <- rbind(cbind(S, sparseMatrix(NULL, NULL, dims=c(nrow(S), n.act.rev))), cbind(m1, m2))
  }

  # 3. df<0 (|v|<|v0|): specify the z+ indicator variables
  # 3a. specify v<|v0|
  rxns.inact <- which(df<0)
  n.inact <- length(rxns.inact)
  #v.inact <- (1-pars$flux.delta.rel)*abs(v0)[rxns.inact]-pars$flux.delta
  v.inact <- (1/pars$flux.delta.rel)*abs(v0)[rxns.inact]-pars$flux.delta
  if (n.inact!=0) {
    m1 <- sparseMatrix(1:n.inact, rxns.inact, dims=c(n.inact, ncol(S)))
    m2 <- Diagonal(n.inact, x=pars$flux.bound-v.inact)
    S <- rbind(cbind(S, sparseMatrix(NULL, NULL, dims=c(nrow(S), n.inact))), cbind(m1, m2))
  }
  # 3b. for those reversible inactive reactions, need to further specify v>-|v0|
  rxns.inact.rev <- which(df<0 & model$lb[1:length(model$rxns)]<0)
  n.inact.rev <- length(rxns.inact.rev)
  #v.inact <- (1-pars$flux.delta.rel)*abs(v0)[rxns.inact.rev]-pars$flux.delta
  v.inact <- (1/pars$flux.delta.rel)*abs(v0)[rxns.inact.rev]-pars$flux.delta
  if (n.inact.rev!=0) {
    m3 <- sparseMatrix(1:n.inact.rev, rxns.inact.rev, dims=c(n.inact.rev, ncol(S)-n.inact))
    m4 <- sparseMatrix(1:n.inact.rev, match(rxns.inact.rev, rxns.inact), x=v.inact-pars$flux.bound, dims=c(n.inact.rev, n.inact))
    S <- rbind(S, cbind(m3, m4))
  }

  # 4. df==0, z0
  rxns0 <- which(df==0)
  n.rxns0 <- length(rxns0)
  #v01 <- ifelse(v0[rxns0]>0, 1+pars$flux.delta.rel, 1-pars$flux.delta.rel)*v0[rxns0]+pars$flux.delta
  #v02 <- ifelse(v0[rxns0]>0, 1-pars$flux.delta.rel, 1+pars$flux.delta.rel)*v0[rxns0]-pars$flux.delta
  v01 <- v0[rxns0]+pars$flux.delta
  v02 <- v0[rxns0]-pars$flux.delta
  if (n.rxns0!=0) {
    m1 <- sparseMatrix(1:n.rxns0, rxns0, dims=c(n.rxns0, ncol(S)))
    m2 <- Diagonal(n.rxns0, x=pars$flux.bound-v01)
    S <- rbind(cbind(S, sparseMatrix(NULL, NULL, dims=c(nrow(S), n.rxns0))), cbind(m1, m2))
    m3 <- sparseMatrix(1:n.rxns0, rxns0, dims=c(n.rxns0, ncol(S)-n.rxns0))
    m4 <- Diagonal(n.rxns0, x=-v02-pars$flux.bound)
    S <- rbind(S, cbind(m3, m4))
  }
  model$S <- S

  # other parameters
  n <- n.act + n.act.rev + n.inact + n.inact.rev + n.rxns0*2
  model$rowlb <- c(model$rowlb, rep(-pars$flux.bound, n))
  model$rowub <- c(model$rowub, rep(pars$flux.bound, n))
  n <- ncol(S) - n.rxns
  model$lb <- c(model$lb, rep(0, n))
  model$ub <- c(model$ub, rep(1, n))
  #if (is.null(w)) w <- c(sum(df!=0,na.rm=TRUE), sum(df==0,na.rm=TRUE))
  if (is.null(w)) w <- c(1,1)
  if ("exprs.int" %in% names(model)) {
    model$c <- c(model$c, rep(c(w[1]/sum(df!=0,na.rm=TRUE), w[2]/sum(df==0,na.rm=TRUE)), c(n-n.rxns0, n.rxns0)))
    model$vtype <- c(model$vtype, rep("I",n))
    model$var.ind <- c(model$var.ind, rep(c("z+","z-","z+1","z0"), c(n.act, n.act.rev, n.inact, n.rxns0))) # iMAT variable type indicators (v: fluxex; z+/-/0: indicator variables)
  } else {
    model$c <- rep(c(0, w[1]/sum(df!=0,na.rm=TRUE), w[2]/sum(df==0,na.rm=TRUE)), c(n.rxns, n-n.rxns0, n.rxns0))
    model$vtype <- ifelse(model$c==0, "C", "I")
    model$var.ind <- rep(c("v","z+","z-","z+1","z0"), c(n.rxns, n.act, n.act.rev, n.inact, n.rxns0)) # iMAT variable type indicators (v: fluxex; z+/-/0: indicator variables)
    model$csense="max"
  }

  model$flux0 <- v0
  model$dfluxes.int <- df
  model$rxns.pos <- rxns.act
  model$rxns.pos.rev <- rxns.act.rev
  model$rxns.neg <- rxns.inact
  model$rxns.neg.rev <- rxns.inact.rev
  model$rxns.steady <- rxns0

  model
}

run.imat.dv0 <- function(imat.model, imat.pars, solv.pars) {
  # solve the iMAT differential flux MILP model
  # return the imat.model with the results from running iMAT added:
  # if the part for the original iMAT formulation is in imat.model, then will contain imat.model$fluxes.int.imat: a vector in the order of the rxns for the two consecutive models, with values 0/9/1/-1 representing a rxn being inactive, activity level not enforced, active in the forward direction, and active in the backward direction, as that returned by run.imat()
  # imat.model$dfluxes.int.imat: a vector in the order of the model rxns, with values 1/2/-1/-2/0/9, where: 1/-1 means positive or negative dflux, 2/-2 means pos/neg dflux of the alternative form (corresponding to (z-) = 1, for reversible rxns only), 0 means steady flux, 9 means rxn is not forced to have either differential or steady flux, as determined by iMAT
  # imat.model$solver.out: the output of MILP solver
  imat.pars <- get.pars("imat", imat.pars)
  solv.pars <- get.pars("mip", solv.pars)

  milp.out <- solve.model(imat.model, pars=solv.pars)
  xopt <- get.xopt(imat.model, milp.out, imat.pars)
  if (any(c("y+","y-","y0") %in% imat.model$var.ind)) imat.model$fluxes.int.imat <- get.imat.opt.flux.int(imat.model, xopt)
  imat.model$dfluxes.int.imat <- get.imat.opt.dv0.int(imat.model, xopt)
  imat.model$solver.out <- milp.out

  imat.model
}

get.imat.opt.dv0.int <- function(imat.model, xopt) {
  # given the formulated imat.model involving differential fluxes, and the xopt of solving that model, return a vector in the order of the model rxns, with values 1/2/-1/0/9, where:
  # 1/-1 means positive or negative dflux, 2 means pos dflux of the alternative form (corresponding to (z-) = 1, for reversible rxns only), 0 means steady flux, 9 means rxn is not forced to have either differential or steady flux

  zp <- xopt[imat.model$var.ind=="z+"]
  zm <- xopt[imat.model$var.ind=="z-"]
  zp1 <- xopt[imat.model$var.ind=="z+1"]
  z0 <- xopt[imat.model$var.ind=="z0"]
  
  pos <- imat.model$rxns.pos[zp==1]
  pos.alt <- imat.model$rxns.pos.rev[zm==1]
  neg <- imat.model$rxns.neg[zp1==1]
  z01 <- imat.model$rxns.steady[z0==1]
  
  zres <- rep(9L, length(imat.model$dfluxes.int))
  zres[pos] <- 1L
  zres[pos.alt] <- 2L
  zres[neg] <- -1L
  zres[z01] <- 0L
  zres
}

update.model.imat.dv0 <- function(model, imat.res, imat.pars) {
  # update metabolic model (model) based on the result of iMAT (differential flux formulation) (imat.res), using parameters given in imat.pars
  # this function will update the parts specifying differential flux constraints; if results of the original iMAT formulation is present (i.e. imat.res$fluxes.int.imat), it will also perform the corresponding updates
  # return the updated model
  pars <- get.pars("imat", imat.pars)

  if ("fluxes.int.imat" %in% names(imat.res)) model <- update.model.imat(model, imat.res, imat.pars)
  dfint <- imat.res$dfluxes.int.imat
  v0 <- imat.res$flux0
  #model$lb[dfint==1] <- (1+pars$flux.delta.rel)*abs(v0)[dfint==1] + pars$flux.delta
  #model$ub[dfint==2] <- -(1+pars$flux.delta.rel)*abs(v0)[dfint==2] - pars$flux.delta
  #model$ub[dfint==-1] <- (1-pars$flux.delta.rel)*abs(v0)[dfint==-1] - pars$flux.delta
  #model$lb[dfint==-1 & model$lb<0] <- -(1-pars$flux.delta.rel)*abs(v0)[dfint==-1 & model$lb<0] + pars$flux.delta
  model$lb[dfint==1] <- (pars$flux.delta.rel)*abs(v0)[dfint==1] + pars$flux.delta
  model$ub[dfint==2] <- -(pars$flux.delta.rel)*abs(v0)[dfint==2] - pars$flux.delta
  model$ub[dfint==-1] <- (1/pars$flux.delta.rel)*abs(v0)[dfint==-1] - pars$flux.delta
  model$lb[dfint==-1 & model$lb<0] <- -(1/pars$flux.delta.rel)*abs(v0)[dfint==-1 & model$lb<0] + pars$flux.delta
  model$lb[dfint==0] <- v0[dfint==0] - pars$flux.delta
  model$ub[dfint==0] <- v0[dfint==0] + pars$flux.delta
  #model$ub[dfint==0] <- ifelse(v0[dfint==0]>0, 1+pars$flux.delta.rel, 1-pars$flux.delta.rel)*v0[dfint==0]+pars$flux.delta
  #model$lb[dfint==0] <- ifelse(v0[dfint==0]>0, 1-pars$flux.delta.rel, 1+pars$flux.delta.rel)*v0[dfint==0]-pars$flux.delta
  
  model
}

