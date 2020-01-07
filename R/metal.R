###### linear programming version of the metabolic transformation algorithm (METAL) ######


### --- metal --- ###

metal <- function(model, flux0, dflux, rxns="all+ctrl", nc=1L, detail=TRUE, solv.pars=list()) {
  # the main function for running metal
  # flux0 is the reference flux vector from sampling an iMAT output model
  # dflux is the flux change, i.e. output from de.dt2dflux(); I make it separate as usually we need to try different parameters in de.dt2dflux()
  # rxns are the indices or IDs or rxns to run metal on, by default all rxns plus the control wild-type model; rxns="all" to run for all rxns w/o the ctrl; if using indices, 0 means ctrl; if using IDs, "ctrl", means ctrl
  # nc: number of cores to use for rxns
  # detail: whether to return more details or only the MTA scores
  # return both the metal.model, and the summary data.table of MTA scores for the rxns, in a list(metal.model, result)

  # formulate metal model
  metal.model <- form.metal(model, flux0, dflux)

  # solve the metal models
  res <- run.ko.screen(metal.model, rxns, run.metal, solv.pars=solv.pars, detail=detail, nc=nc)

  list(metal.model=metal.model, result=res)
}


### --- helper functions for the individual internal steps of metal --- ###

form.metal <- function(model, flux0, dflux) {
  # formulate metal model (no extra parameters needed)

  n.mets <- nrow(model$S)
  n.rxns <- ncol(model$S)

  # constraint matrix
  ## for reactions meant to remain steady
  st <- which(dflux==0 & model$c!=1 & !grepl("biomass", model$rxns, ignore.case=TRUE))
  n.st <- length(st)
  S <- rbind(
    cbind( model$S,                                              sparseMatrix(NULL, NULL, dims=c(n.mets, 2*n.st)) ),
    cbind( sparseMatrix(1:n.st, st, x=1, dims=c(n.st, n.rxns)),  Diagonal(n.st),  Diagonal(n.st, -1) )
  )
  ## for **reversible** reactions that is meant to have reduced fluxes (i.e. these have the potential to "overshoot", thus need to be treated specially)
  rvdn.b <- model$lb<0 & dflux<0
  rvdn <- which(rvdn.b)
  n.rvdn <- length(rvdn)
  if (n.rvdn>0) {
    S <- rbind(
      cbind( S,                                                                                                             sparseMatrix(NULL, NULL, dims=c(n.mets+n.st, 2*n.rvdn)) ),
      cbind( sparseMatrix(1:n.rvdn, rvdn, x=1, dims=c(n.rvdn, n.rxns)),  sparseMatrix(NULL, NULL, dims=c(n.rvdn, 2*n.st)),  Diagonal(n.rvdn),  Diagonal(n.rvdn, -1) )
    )
  }

  # constraints
  rowlb <- c(model$rowlb, flux0[st], rep(0, n.rvdn))
  rowub <- c(model$rowub, flux0[st], rep(0, n.rvdn))
  lb <- c(model$lb, rep(0, 2*n.st+2*n.rvdn))
  ub <- c(model$ub, rep(Inf, 2*n.st+2*n.rvdn))

  # objective function and others
  ## reactions meant to change in the forward direction, excluding those in rvdn (i.e. w/o the potential to "overshoot")
  fw0.b <- flux0>0 & dflux>0 | flux0==0 & dflux>0 & model$lb>=0
  ## reactions meant to change in the backward direction, excluding those in rvdn (i.e. w/o the potential to "overshoot")
  bk0.b <- flux0<0 & dflux>0 | flux0>0 & dflux<0 & model$lb>=0
  w <- abs(dflux) / sum(abs(dflux), na.rm=TRUE) # weight
  c <- c(ifelse(fw0.b, -w, ifelse(bk0.b, w, 0)), rep(1/n.st, 2*n.st), w[rvdn], w[rvdn])
  c[is.na(c)] <- 0

  # things to keep for downstream analysis of MTA score
  ## reactions meant to change in the forward direction
  fw <- which(flux0>0 & dflux>0 | flux0<0 & dflux<0)
  ## reactions meant to change in the backward direction
  bk <- which(flux0<0 & dflux>0 | flux0>0 & dflux<0)
  ## **reversible** reactions with flux0==0 and dflux==1: it's fine for them to change either forward or backward; Note that I did not include these in the objective function: they would require maximize |v|, which is slightly tricker; I simply neglect such cases and adjust for them later in the scores
  fw.or.bk <- which(flux0==0 & dflux>0 & model$lb<0)

  # return MTA model
  list(flux0=flux0, dflux=dflux, w=w,
       fw0.b=fw0.b, bk0.b=bk0.b, rvdn.b=rvdn.b, fw=fw, bk=bk, fw.or.bk=fw.or.bk, st=st, n.st=n.st,
       rxns=model$rxns, mets=model$mets, csense="min", c=c, S=S, rowlb=rowlb, rowub=rowub, lb=lb, ub=ub)
}

get.metal.score <- function(model, lp.out, detail) {
  # given one metal model and the output from solving the model, compute the MTA score and other associated info
  # return a 1-row data.table

  lp.out <- lp.out[[1]]
  if (is.na(lp.out$xopt)) {
    if (detail) {
      return(data.table(solv.stat=lp.out$stat.str, v.opt=NA, rxns.change.yes=NA, rxns.change.no=NA, advs.change.yes=NA, advs.change.no=NA, advs.steady=NA, score.change=NA, score.steady=NA, score.mta=NA))
    } else {
      return(data.table(solv.stat=lp.out$stat.str, score.change=NA, score.steady=NA, score.mta=NA))
    }
  }

  v0 <- model$flux0
  v <- lp.out$xopt[1:length(v0)]
  fw0 <- model$fw0.b
  bk0 <- model$bk0.b
  rvdn <- model$rvdn.b

  # reactions intended to change: successful
  fw0.yes <- fw0 & v>v0
  bk0.yes <- bk0 & v<v0
  rvdn.yes <- rvdn & abs(v)<abs(v0)
  yes <- which(fw0.yes|bk0.yes|rvdn.yes)

  # reactions intended to change: failed
  fw0.no <- fw0 & v<v0
  bk0.no <- bk0 & v>v0
  rvdn.no <- rvdn & abs(v)>abs(v0)
  no <- which(fw0.no|bk0.no|rvdn.no)

  # absolute differences between v and flux0 for the different sets of reactions
  ## reactions intended to change
  adv.yes <- ifelse(yes %in% which(rvdn.yes), abs(v0[yes])-abs(v[yes]), abs(v[yes]-v0[yes]))
  adv.no <- ifelse(no %in% which(rvdn.no), abs(v[no])-abs(v0[no]), abs(v[no]-v0[no]))
  ## reactions intended to remain steady
  adv.st <- abs(v[model$st]-v0[model$st])

  # score for reactions intended to change
  s.ch <- sum(adv.yes * model$w[yes]) + sum(v[model$fw.or.bk] * model$w[model$fw.or.bk]) - sum(adv.no * model$w[no]) # score.yes - score.no
  # score for reactions intended to remain steady
  s.st.uw <- sum(adv.st) # un-weighted
  s.st <- s.st.uw / model$n.st
  # raw score: just the (negated) optimal objective value
  #s.raw <- -lp.out$obj
  # adjusted score
  #s.adj <- s.raw + sum(v0[model$bk] * model$w[model$bk]) - sum(v0[model$fw] * model$w[model$fw]) + sum(v[model$fw.or.bk] * model$w[model$fw.or.bk])
  # recalculated score (should be the same as the adjusted score) ## yes, same
  #s.rec <- s.ch - s.st
  # ratio score like the original mta
  s.mta <- s.ch / s.st

  # return
  if (detail) {
    res <- data.table(solv.stat=lp.out$stat.str, v.opt=list(v), rxns.change.yes=list(yes), rxns.change.no=list(no), advs.change.yes=list(adv.yes), advs.change.no=list(adv.no), advs.steady=list(adv.st), score.change=s.ch, score.steady=s.st, score.mta=s.mta)
  } else {
    res <- data.table(solv.stat=lp.out$stat.str, score.change=s.ch, score.steady=s.st, score.mta=s.mta)
  }
  res
}

run.metal <- function(model, solv.pars, detail) {
  # solve the metal models for each rxn
  # return the summary of MTA scores for the rxns in a data.table
  solv.pars <- get.pars("lp", solv.pars)

  solv.out <- solve.model(model, pars=solv.pars)
  get.metal.score(model, solv.out, detail)
}


### --- additional variations of metal --- ###

mmetal <- function(model, flux0, dflux, rxns="all", nc=1L, n=1, seeds=1:n, detail=TRUE, solv.pars=list()) {
  # running multiple metal models including the original one, with the additional "control" models using: 1. dflux <- -dflux; 2. dflux <- a set of random orthogonal dflux vectors
  # n: the number of random orthogonal dflux vectors to generate (and thus the random metal models to run)
  # seeds: the seeds for generating the random orthogonal dflux vectors, its length is equal to n; or NULL meaning do not set seed
  # return like the usual metal model, but with additional columns in the summary table for "meta" mta scores:
  # score.meta1: the original mta score minus the mta score from using -dflux; score.meta2: the original score minus the median score across all additional control models

  if (n==0) {
    dfs <- list(-dflux)
  } else if (n>0) {
    if (is.null(seeds)) tmp <- lapply(1:n, function(i) get.ortho.vec(dflux)) else tmp <- lapply(seeds, function(i) get.ortho.vec(dflux, i))
    dfs <- c(list(-dflux), tmp)
  }

  # original model
  res <- metal(model, flux0, dflux, rxns, nc, detail, solv.pars)

  # control models
  ctrls <- sapply(dfs, function(x) {
    res <- metal(model, flux0, x, rxns, nc, detail=FALSE, solv.pars)
    res$result$score.mta
  })
  med <- apply(ctrls, 1, median)
  # add "meta" mta scores to res$scores
  res$result[, c("score.meta1","score.meta2"):=list(score.mta-ctrls[,1], score.mta-med)]
  res
}

get.ortho.vec <- function(x, seed=NULL) {
  # generate random orthogonal vectors to a given vector x that contains discretized values from {-1,0,1}

  # shuffle x
  set.seed(seed); y <- sample(x)
  # current inner product
  k <- sum(x*y)
  # if k happens to be 0, simply return y
  if (k==0) return(y)
  # if k is odd, we randomly swap a pair of y's where y==0 && x!=0 and y!=0 && x==0, so that we reduce one product-is-zero case; if it happens that no such case exists, we randomly swap a pair of y's where y==0 && x==0 and y!==0 && x!==0, so that we increase on product-is-zero case
  if (k%%2==1) {
    if (any(y==0 & x!=0)) {
      set.seed(seed+1); p1 <- sample(which(y==0 & x!=0), 1)
      set.seed(seed+2); p2 <- sample(which(y!=0 & x==0), 1)
    } else {
      set.seed(seed+3); p1 <- sample(which(y==0 & x==0), 1)
      set.seed(seed+4); p2 <- sample(which(y!=0 & x!=0), 1)
    }
    tmp <- y[p1]
    y[p1] <- y[p2]
    y[p2] <- tmp
  }
  k <- sum(x*y)
  # now the k should become even; randomly "flip" k/2 product-non-zero positions in the proper direction
  set.seed(seed+5); fp <- sample(which(x*y==sign(k)), abs(k)/2)
  y[fp] <- -y[fp]
  # then randomly "flip" equal number of positions of product +1 and -1
  set.seed(seed+6); nf <- sample(0:(sum(x*y!=0)/2), 1)
  set.seed(seed+7); fpp <- sample(which(x*y==1), nf)
  set.seed(seed+8); fpn <- sample(which(x*y==-1), nf)
  fp <- c(fpp, fpn)
  y[fp] <- -y[fp]
  y
}

rmetal <- function(model, flux0, dflux, rxns="all+ctrl", nc=1L, detail=TRUE, k=100, lp.pars=list(), qp.pars=list()) {
  # the function for running a "robust" version of metal similary to rMTA (basically, rMTA with metal rather than the original MIQP version of MTA); k is the rMTA-specific parameter
  # flux0 is the reference flux vector from sampling an iMAT output model
  # dflux is the flux change, i.e. output from de.dt2dflux(); I make it separate as usually we need to try different parameters in de.dt2dflux()
  # rxns are the indices or IDs or rxns to run MTA on, by default all rxns plus the control wild-type model; rxns="all" to run for all rxns w/o the ctrl; if using indices, 0 means ctrl; if using IDs, "ctrl", means ctrl
  # nc: number of cores to use for rxns
  # detail: whether to return more details or only the MTA scores
  # return both the mta.model, and the summary data.table of MTA scores for the rxns, in a list(mta.model, scores)

  # formulate metal model for either direction (dflux and -dflux)
  metal.model <- form.metal(model, flux0, dflux)
  metal.model0 <- form.metal(model, flux0, -dflux)
  # solve the MTA models across the rxns for both models
  res <- run.ko.screen(metal.model, rxns, run.metal, solv.pars=lp.pars, detail=detail, nc=nc)
  res0 <- run.ko.screen(metal.model0, rxns, run.metal, solv.pars=lp.pars, detail=FALSE, nc=nc)

  # MOMA for dflux
  tmpf <- function(x) get.metal.score(model=metal.model, x, detail=detail)
  res.moma <- moma(model, rxns, nc, flux0, obj=tmpf, solv.pars=qp.pars)

  # rMTA score
  res <- data.table(id=res$id, rxn=res$rxn, bTS=res$score.mta, wTS=res0$score.mta, mTS=res.moma$score.mta)
  res[, rTS:=ifelse(bTS>0 & mTS>0 & wTS<0, mta.pars$k*mTS*(bTS-wTS), mTS)]

  list(metal.model=metal.model, result.metal=res, result.moma=res.moma, result.rmetal=res)
}



