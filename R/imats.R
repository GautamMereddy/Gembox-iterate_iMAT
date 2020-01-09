###### variations of the iMAT algorithm for multiple samples or cells (iMATs) ######


### --- iMAT for two samples --- ###

imat2 <- function(model, expr1, expr2, dflux, imat.pars=list(), solv.pars=list(), solv.pars1=list(), samp.pars=list()) {
  # the main function for running iMAT for two samples in a single step
  # dflux is the discretized (-1/0/1) vector in the order or model rxns representing the intended flux changes of model2 compared to model1
  # if imat.pars$nsteps==1, return a list(imat.model, result.model)
  # if imat.pars$nsteps==2, return a list(imat.model1, imat.model2, imat.model.dflux, result.model), where imat.model1/2 are the original imat models from the first step for expr1/2, and imat.model.dflux is the diff flux part of the imat model in the second step

  imat.pars <- get.pars("imat", imat.pars)
  if (imat.pars$nsteps==1) {
    # if running imat2 in a single step
    # original imat model for expr1 and expr2
    model1 <- form.imat(model, expr1, imat.pars)
    model2 <- form.imat(model, expr2, imat.pars)
  } else if (imat.pars$nsteps==2) {
    # if running imat2 in two steps, as the 1st step, run the original imat, obtain the updated metabolic models for both samples
    tmp.model1 <- imat(model, expr1, imat.pars, solv.pars1, samp.pars=NULL)
    tmp.model2 <- imat(model, expr2, imat.pars, solv.pars1, samp.pars=NULL)
    model1 <- tmp.model1$result.model
    model2 <- tmp.model2$result.model
  }
  # formulate diff flux part of the model
  imat.model <- form.imat.dflux(model1, model2, dflux, imat.pars)
  # solve the iMAT model
  imat.model <- run.imat.dflux(imat.model, imat.pars, solv.pars)
  # update the original metabolic model based on iMAT result
  modelx2 <- c.model(model, model)
  res.model <- update.model.imat.dflux(modelx2, imat.model, imat.pars)

  if (!is.null(samp.pars)) {
    # sample the updated metabolic model
    res.model <- sample.model(res.model, samp.pars)
  }

  # if two steps, include both the original iMAT models for expr1 and expr2 and the diff flux iMAT model in result
  if (imat.pars$nsteps==1) {
    return(list(imat.model=imat.model, result.model=res.model))
  } else if (imat.pars$nsteps==2) {
    return(list(imat.model1=tmp.model1$imat.model, imat.model2=tmp.model2$imat.model, imat.model.dflux=imat.model, result.model=res.model))
  }

}


### --- helper functions for iMAT for two samples --- ###

form.imat.dflux0 <- function(model, i1, i2, df, rr, pars) {
  # a helper function to formulate the differential flux part of iMATs
  # this function adds model constrains for the df of a single pair of rxns; i1 and i2 are their indices, df is the diff flux; rr is T/F to mark whether the pair or rxns are reversible
  # pars: a list of parameters (imat.pars)

  S <- model$S
  if (df==0) { # z0
    S <- rbind(cbind(S, sparseMatrix(NULL,NULL,dims=c(nrow(S),1))),
               sparseMatrix(rep(1,3), c(i1, i2, ncol(S)+1), x=c(1, -1, 2*pars$flux.bound-pars$flux.delta), dims=c(1, ncol(S)+1)),
               sparseMatrix(rep(1,3), c(i1, i2, ncol(S)+1), x=c(1, -1, -2*pars$flux.bound+pars$flux.delta), dims=c(1, ncol(S)+1)))
    model$rowlb <- c(model$rowlb, rep(-2*pars$flux.bound, 2))
    model$rowub <- c(model$rowub, rep(2*pars$flux.bound, 2))
    model$lb <- c(model$lb, 0)
    model$ub <- c(model$ub, 1)
    model$c <- c(model$c, 1)
    model$vtype <- c(model$vtype, "I")
    model$var.ind <- c(model$var.ind, "z0")
  } else {
    # z+
    if (df>0) {
      S <- rbind(cbind(S, sparseMatrix(NULL,NULL,dims=c(nrow(S),1))),
                 sparseMatrix(rep(1,3), c(i1, i2, ncol(S)+1), x=c(1, pars$flux.delta.rel-1, (2-pars$flux.delta.rel)*pars$flux.bound+pars$flux.delta), dims=c(1,ncol(S)+1)))
    } else if (df<0) {
      S <- rbind(cbind(S, sparseMatrix(NULL,NULL,dims=c(nrow(S),1))),
                 sparseMatrix(rep(1,3), c(i1, i2, ncol(S)+1), x=c(pars$flux.delta.rel-1, 1, (2-pars$flux.delta.rel)*pars$flux.bound+pars$flux.delta), dims=c(1,ncol(S)+1)))
    }
    model$rowlb <- c(model$rowlb, (pars$flux.delta.rel-2)*pars$flux.bound)
    model$rowub <- c(model$rowub, (2-pars$flux.delta.rel)*pars$flux.bound)
    model$lb <- c(model$lb, 0)
    model$ub <- c(model$ub, 1)
    model$c <- c(model$c, 1)
    model$vtype <- c(model$vtype, "I")
    model$var.ind <- c(model$var.ind, "z+")
    # if reversible rxn, need to add an extra constraint on z+, and similarly a pair of constraints on z-, and the constraint that (z+) + (z-) <= 1
    if (rr) {
      if (df>0) {
        # additional z+
        S <- rbind(S, sparseMatrix(rep(1,3), c(i1, i2, ncol(S)), x=c(1, 1-pars$flux.delta.rel, (pars$flux.delta.rel-2)*pars$flux.bound-pars$flux.delta), dims=c(1,ncol(S))))
        # z-
        S <- rbind(cbind(S, sparseMatrix(NULL,NULL,dims=c(nrow(S),1))),
                   sparseMatrix(rep(1,3), c(i1, i2, ncol(S)+1), x=c(1, 1-pars$flux.delta.rel, (2-pars$flux.delta.rel)*pars$flux.bound+pars$flux.delta), dims=c(1,ncol(S)+1)),
                   sparseMatrix(rep(1,3), c(i1, i2, ncol(S)+1), x=c(1, pars$flux.delta.rel-1, (pars$flux.delta.rel-2)*pars$flux.bound-pars$flux.delta), dims=c(1,ncol(S)+1)))
      } else if (df<0) {
        # additional z+
        S <- rbind(S, sparseMatrix(rep(1,3), c(i1, i2, ncol(S)), x=c(1-pars$flux.delta.rel, 1, (pars$flux.delta.rel-2)*pars$flux.bound-pars$flux.delta), dims=c(1,ncol(S))))
        # z-
        S <- rbind(cbind(S, sparseMatrix(NULL,NULL,dims=c(nrow(S),1))),
                   sparseMatrix(rep(1,3), c(i1, i2, ncol(S)+1), x=c(1-pars$flux.delta.rel, 1, (2-pars$flux.delta.rel)*pars$flux.bound+pars$flux.delta), dims=c(1,ncol(S)+1)),
                   sparseMatrix(rep(1,3), c(i1, i2, ncol(S)+1), x=c(pars$flux.delta.rel-1, 1, (pars$flux.delta.rel-2)*pars$flux.bound-pars$flux.delta), dims=c(1,ncol(S)+1)))
      }
      # (z+) + (z-) <= 1
      S <- rbind(S, sparseMatrix(rep(1,2), c(ncol(S)-1, ncol(S)), dims=c(1,ncol(S))))
      model$rowlb <- c(model$rowlb, rep((pars$flux.delta.rel-2)*pars$flux.bound, 3), 0)
      model$rowub <- c(model$rowub, rep((2-pars$flux.delta.rel)*pars$flux.bound, 3), 1)
      model$lb <- c(model$lb, 0)
      model$ub <- c(model$ub, 1)
      model$c <- c(model$c, 1)
      model$vtype <- c(model$vtype, "I")
      model$var.ind <- c(model$var.ind, "z-")
    }
  }
  model$S <- S
  model
}

form.imat.dflux <- function(model1, model2=model1, dflux, imat.pars) {
  # formulate the differential flux part of iMATs between two models (i.e. two samples)
  # model1 and model2 can be metabolic models (either "raw" or updated models from imat()), or iMAT models (i.e. output from form.imat())
  # dflux is the discretized (-1/0/1) vector in the order or model rxns representing the intended flux changes of model2 compared to model1
  pars <- get.pars("imat", imat.pars)

  model <- c.model(model1, model2, c.vars=c("c","vtype"), cx.vars="var.ind")
  # need to add a few things to model
  if ("exprs.int" %in% names(model1)) {
    # if model1 (and assume model2 as well) is an imat model from form.imat()
    model$exprs.int1 <- model1$exprs.int
    model$fluxes.int1 <- model1$fluxes.int
    model$exprs.int2 <- model2$exprs.int
    model$fluxes.int2 <- model2$fluxes.int
    n.rxns1 <- length(model1$rxns)
    model$rxns.act <- c(model1$rxns.act, model2$rxns.act+n.rxns1)
    model$rxns.act.rev <- c(model1$rxns.act.rev, model2$rxns.act.rev+n.rxns1)
    model$rxns.inact <- c(model1$rxns.inact, model2$rxns.inact+n.rxns1)
    model$rxns.inact.rev <- c(model1$rxns.inact.rev, model2$rxns.inact.rev+n.rxns1)
  } else {
    # if model1 (and assume model2 as well) is a metabolic model
    model$c <- rep(0, ncol(model$S))
    model$vtype <- rep("C", ncol(model$S))
    model$var.ind <- c(rep("v_1", ncol(model1$S)), rep("v_2", ncol(model2$S)))
  }
  model$dfluxes.int <- dflux
  model$rxns.pos <- which(dflux>0)
  model$rxns.neg <- which(dflux<0)
  model$rxns.change <- which(dflux!=0)
  model$rxns.change.rev <- which(dflux!=0 & model1$lb[1:length(model1$rxns)]<0)
  model$rxns.steady <- which(dflux==0)
  model$csense <- "max"

  nc1 <- ncol(model1$S)
  for (i in 1:length(dflux)) {
    rr <- model1$lb[i]<0 # use model1 is OK since even if it's updated by iMAT the reversibilities of reactions don't change
    model <- form.imat.dflux0(model, i, nc1+i, dflux[i], rr, pars)
  }
  model
}

run.imat.dflux <- function(imat.model, imat.pars, solv.pars) {
  # solve the iMAT differential flux MILP model
  # return the imat.model with the results from running iMAT added:
  # if the part for the original iMAT formulation is in imat.model, then will contain imat.model$fluxes.int.imat: a vector in the order of the rxns for the two consecutive models, with values 0/9/1/-1 representing a rxn being inactive, activity level not enforced, active in the forward direction, and active in the backward direction, as that returned by run.imat()
  # imat.model$dfluxes.int.imat: a vector in the order of the model rxns, with values 1/2/-1/-2/0/9, where: 1/-1 means positive or negative dflux, 2/-2 means pos/neg dflux of the alternative form (corresponding to (z-) = 1, for reversible rxns only), 0 means steady flux, 9 means rxn is not forced to have either differential or steady flux, as determined by iMAT
  # imat.model$solver.out: the output of MILP solver
  imat.pars <- get.pars("imat", imat.pars)
  solv.pars <- get.pars("mip", solv.pars)

  milp.out <- solve.model(imat.model, pars=solv.pars)
  xopt <- get.xopt(imat.model, milp.out, imat.pars)
  if (any(c("y+_1","y-_1","y0_1") %in% imat.model$var.ind)) imat.model$fluxes.int.imat <- get.imat.opt.flux.int(imat.model, xopt)
  imat.model$dfluxes.int.imat <- get.imat.opt.dflux.int(imat.model, xopt)
  imat.model$solver.out <- milp.out

  imat.model
}

get.imat.opt.dflux.int <- function(imat.model, xopt) {
  # given the formulated imat.model involving differential fluxes, and the xopt of solving that model, return a vector in the order of the model rxns, with values 1/2/-1/-2/0/9, where:
  # 1/-1 means positive or negative dflux, 2/-2 means pos/neg dflux of the alternative form (corresponding to (z-) = 1, for reversible rxns only), 0 means steady flux, 9 means rxn is not forced to have either differential or steady flux

  zp <- xopt[imat.model$var.ind=="z+"]
  zm <- xopt[imat.model$var.ind=="z-"]
  z0 <- xopt[imat.model$var.ind=="z0"]
  zp1 <- imat.model$rxns.change[zp==1]
  zm1 <- imat.model$rxns.change.rev[zm==1]
  z01 <- imat.model$rxns.steady[z0==1]
  pos <- intersect(zp1, imat.model$rxns.pos) # rxns whose fluxes 2>1 (the (z+)=1 case)
  neg <- intersect(zp1, imat.model$rxns.neg) # rxns whose fluxes 2<1 (the (z+)=1 case)
  pos.alt <- intersect(zm1, imat.model$rxns.pos) # rxns whose fluxes 2>1, alternative form (i.e. (z-)=1 case)
  neg.alt <- intersect(zm1, imat.model$rxns.neg) # rxns whose fluxes 2<1, alternative form (i.e. (z-)=1 case)

  zres <- rep(9L, length(imat.model$dfluxes.int))
  zres[pos] <- 1L
  zres[pos.alt] <- 2L
  zres[neg] <- -1L
  zres[neg.alt] <- -2L
  zres[z01] <- 0L
  zres
}

update.model.imat.dflux <- function(model, imat.res, imat.pars) {
  # update metabolic model (model) based on the result of iMAT (differential flux formulation) (imat.res), using parameters given in imat.pars
  # this function will update the parts specifying differential flux constraints; if results of the original iMAT formulation is present (i.e. imat.res$fluxes.int.imat), it will also perform the corresponding updates
  # return the updated model
  pars <- get.pars("imat", imat.pars)

  if ("fluxes.int.imat" %in% names(imat.res)) model <- update.model.imat(model, imat.res, imat.pars)
  dfint <- imat.res$dfluxes.int.imat
  n.rxns <- length(dfint)
  for (i in 1:n.rxns) {
    if (dfint[i]==1) {
      model <- add.constraint(model, c(i,i+n.rxns), c(1,pars$flux.delta.rel-1), (pars$flux.delta.rel-2)*pars$flux.bound, -pars$flux.delta)
      if (model$lb[i]<0) model <- add.constraint(model, c(i,i+n.rxns), c(1,1-pars$flux.delta.rel), pars$flux.delta, (2-pars$flux.delta.rel)*pars$flux.bound)
    } else if (dfint[i]==-1) {
      model <- add.constraint(model, c(i,i+n.rxns), c(pars$flux.delta.rel-1,1), (pars$flux.delta.rel-2)*pars$flux.bound, -pars$flux.delta)
      if (model$lb[i]<0) model <- add.constraint(model, c(i,i+n.rxns), c(1-pars$flux.delta.rel,1), pars$flux.delta, (2-pars$flux.delta.rel)*pars$flux.bound)
    } else if (dfint[i]==2) {
      model <- add.constraint(model, c(i,i+n.rxns), c(1,1-pars$flux.delta.rel), (pars$flux.delta.rel-2)*pars$flux.bound, -pars$flux.delta)
      model <- add.constraint(model, c(i,i+n.rxns), c(1,pars$flux.delta.rel-1), pars$flux.delta, (2-pars$flux.delta.rel)*pars$flux.bound)
    } else if (dfint[i]==-2) {
      model <- add.constraint(model, c(i,i+n.rxns), c(1-pars$flux.delta.rel,1), (pars$flux.delta.rel-2)*pars$flux.bound, -pars$flux.delta)
      model <- add.constraint(model, c(i,i+n.rxns), c(pars$flux.delta.rel-1,1), pars$flux.delta, (2-pars$flux.delta.rel)*pars$flux.bound)
    } else if (dfint[i]==0) {
      model <- add.constraint(model, c(i,i+n.rxns), c(1,-1), -2*pars$flux.bound, pars$flux.delta)
      model <- add.constraint(model, c(i,i+n.rxns), c(1,-1), -pars$flux.delta, 2*pars$flux.bound)
    }
  }

  model
}


### --- iMAT for multicelllar model of two or three cells --- ###

imat.mc <- function(model, expr, dflux, imat.pars=list(), solv.pars=list(), solv.pars1=list(), samp.pars=list()) {
  # the main function for running iMAT for multicellular model of two or three cells
  # if imat.pars$nsteps==1, return a list(imat.model, result.model)
  # if imat.pars$nsteps==2, return a list(imat.model, imat.model.dflux, result.model), where imat.model is the original imat model from the first step, and imat.model.dflux is the diff flux part of the imat model in the second step

  imat.pars <- get.pars("imat", imat.pars)
  if (imat.pars$nsteps==1) {
    # formulate iMAT model including the diff flux part
    imat.model <- form.imat(model, expr, imat.pars)
    imat.model <- form.imat.dflux.mc(imat.model, dflux, imat.pars)
  } else if (imat.pars$nsteps==2) {
    # first run the original iMAT and get the updated metabolic model
    tmp.model <- imat(model, expr, imat.pars, solv.pars1, samp.pars=NULL)
    # then formulate the diff flux part of the iMAT model from the above updated metabolic model
    imat.model <- form.imat.dflux.mc(tmp.model$result.model, dflux, imat.pars)
  }

  # solve the iMAT model and also obtain the updated metabolic model
  res <- run.imat.mc(imat.model, imat.pars, solv.pars)

  # if two steps, include both the original iMAT model and the diff flux iMAT model in result
  if (imat.pars$nsteps==2) {
    res$imat.model.dflux <- res$imat.model
    res$imat.model <- tmp.model$imat.model
  }

  if (!is.null(samp.pars)) {
    # sample the updated metabolic model
    res$result.model <- sample.model(res$result.model, samp.pars)
  }

  res
}


### --- helper functions for iMAT for multicellar model of two or three cells --- ###

form.imat.dflux.mc <- function(model, dflux, imat.pars=list()) {
  # formulate the differential flux part of iMATs among cells for a multicellular model (two or three cells)
  # model can be a multicellular metabolic model (either "raw" or updated model from imat()), or an iMAT model (i.e. output from form.imat())
  # for model of two cells, dflux is the discretized (-1/0/1) vector in the order or model rxns representing the intended flux changes of cell2 compared to cell1
  # for model of three cells, dflux is a list of length 3 named by "d12", "d23", "d13", containing the dflux vectors as above for cell2-vs-cell1, 3-vs-2, and 3-vs-1, respectively
  pars <- get.pars("imat", imat.pars)

  res.model <- model
  # if model is not an imat model from form.imat(), then need to add a few things
  if (!"exprs.int" %in% names(model)) {
    res.model$c <- rep(0, ncol(model$S))
    res.model$vtype <- rep("C", ncol(model$S))
    res.model$var.ind <- rep("v", ncol(model$S))
    res.model$csense <- "max"
  }
  irxns <- res.model$irxn.ids # irxns.ids should be contained as a part of a multicellular metabolic model
  if (is.list(dflux)) dflux <- lapply(dflux, function(x) x[irxns]) else dflux <- dflux[irxns]
  n <- length(irxns)
  nc <- sum(res.model$vtype=="C") - n

  # model of two cells
  if (!is.list(dflux)) {
    for (i in 1:n) {
      rr <- model$lb[irxns[i]]<0
      i1 <- irxns[i]
      i2 <- nc + i
      df <- dflux[i]
      res.model <- form.imat.dflux0(res.model, i1, i2, df, rr, pars)
    }
  }

  # model of three cells
  if (is.list(dflux) && length(dflux)==3) {
    for (i in 1:n) {
      rr <- model$lb[irxns[i]]<0
      # 1--2
      i1 <- irxns[i]
      i2 <- nc - n + i
      df <- dflux$d12[i]
      res.model <- form.imat.dflux0(res.model, i1, i2, df, rr, pars)
      # 2--3
      i1 <- nc - n + i
      i2 <- nc + i
      df <- dflux$d23[i]
      res.model <- form.imat.dflux0(res.model, i1, i2, df, rr, pars)
      # 1--3
      i1 <- irxns[i]
      i2 <- nc + i
      df <- dflux$d13[i]
      res.model <- form.imat.dflux0(res.model, i1, i2, df, rr, pars)
    }
  }

  res.model
}

run.imat.mc <- function(imat.model, imat.pars, solv.pars) {
  # solve the multicellular iMAT MILP model with differential flux formulation
  # return a list(imat.model, result.model), where imat.model has the the output of MILP solver added as imat.model$solver.out, and result.model is the updated metabolic model
  imat.pars <- get.pars("imat", imat.pars)
  solv.pars <- get.pars("mip", solv.pars)

  imat.model$solver.out <- solve.model(imat.model, pars=solv.pars)
  xopt <- get.xopt(imat.model, imat.model$solver.out, imat.pars)

  # the de part of the model
  # rows and cols to keep from imat.model$S: the "v" part and the z==1 part.
  tmp <- imat.model$S[, (imat.model$var.ind %in% c("z+","z-","z0") & xopt==0) | imat.model$var.ind %in% c("y+","y-","y0")]
  rind <- apply(tmp, 1, function(x) all(x==0))
  cind <- imat.model$var.ind=="v"
  x <- rowSums(imat.model$S[rind, imat.model$var.ind %in% c("z+","z-","z0") & xopt==1])
  res.model <- subset.model(imat.model, rind, cind)
  res.model$rowlb <- res.model$rowlb - x
  res.model$rowub <- res.model$rowub - x

  # the original imat part of the model, if present
  if ("y+" %in% imat.model$var.ind) {
    yp <- xopt[imat.model$var.ind=="y+"]
    ym <- xopt[imat.model$var.ind=="y-"]
    y0 <- xopt[imat.model$var.ind=="y0"]

    fw <- imat.model$rxns.act[yp==1]
    bk <- imat.model$rxns.act.rev[ym==1]
    inact <- imat.model$rxns.inact[y0==1]
    inact.rev <- intersect(inact, imat.model$rxns.inact.rev)

    # update model
    res.model$lb[fw] <- imat.pars$flux.act
    res.model$ub[bk] <- -imat.pars$flux.act
    res.model$ub[inact] <- imat.pars$flux.inact
    res.model$lb[inact.rev] <- -imat.pars$flux.inact
  }

  list(imat.model=imat.model, result.model=res.model)
}


