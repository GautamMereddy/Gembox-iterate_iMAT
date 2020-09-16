###### the metabolic transformation algorithm (MTA, the original MIQP version) ######


### --- prepare data for MTA --- ###

de.dt2dflux <- function(model, de.res, topn=100, padj.cutoff=0.05, na2zero=TRUE, old=TRUE) {
  # map a differential expression analysis result given in de.res to differential fluxes of reactions in a metabolic model
  # de.res: differential expresssion result as a data.table, with the columns named id (gene symbols), log.fc (log fold-change) and pval (raw p value)
  # topn and padj.cutoff: at most top topn DE genes in either direction with BH-adjusted p value < padj.cutoff are kept as the changed genes (i.e. with values -1/1), all other genes are set to 0 (unchanged)
  # the resulting discretized vector is passed to de2dflux to convert to discretized dflux values
  # if na2zero, NA's will be replaced by zeros in the final result

  de.res <- de.res[id %in% model$genes]
  de.res[, padj:=p.adjust(pval, method="BH")]
  gns.ch <- c(de.res[padj<padj.cutoff & log.fc<0][order(log.fc)][1:min(.N,topn), id], de.res[padj<padj.cutoff & log.fc>0][order(-log.fc)][1:min(.N,topn), id])
  de.res[, df:=ifelse(id %in% gns.ch, sign(log.fc), 0)]
  de.int <- de.res$df
  names(de.int) <- de.res$id

  # map to reactions
  df <- de2dflux(model, de.int, na.after=ifelse(na2zero,0,NA))
  if (old) df[model$rules==""] <- 0 # I later changed de2dflux such that rxns w/o genes are NA; previously such cases are 0; I added this to keep the old behaviour here for MTA
  npos <- sum(df>0, na.rm=TRUE)
  nneg <- sum(df<0, na.rm=TRUE)
  if (is.infinite(topn)) {
    message(sprintf("de.dt2dflux(): All the DE genes with p.adj<%g in the model selected, mapped to %d up-regulated reactions and %d down-regulated reactions.\n", padj.cutoff, npos, nneg))
  } else message(sprintf("de.dt2dflux(): At most top %g DE genes with p.adj<%g in the model selected, mapped to %d up-regulated reactions and %d down-regulated reactions.\n", topn, padj.cutoff, npos, nneg))
  df
}


### --- MTA --- ###

mta <- function(model, flux0, dflux, rxns="all+ctrl", ko=NULL, nc=1L, detail=TRUE, mta.pars=list(), solv.pars=list()) {
  # the main function for running MTA (the original MIQP version)
  # flux0 is the reference flux vector from sampling an iMAT output model
  # dflux is the flux change, i.e. output from de.dt2dflux(); I make it separate as usually we need to try different parameters in de.dt2dflux()
  # rxns are the indices or IDs or rxns to run MTA on, by default all rxns plus the control wild-type model; rxns="all" to run for all rxns w/o the ctrl; if using indices, 0 means ctrl; if using IDs, "ctrl", means ctrl
  # ko: NULL; or rxn indices or IDs to KO to combine with those in rxns -- these KO's will be added before screening for those in rxns, but after the metal.model has been formed using the original wildtype model (this is the reasonable way since an existent KO may change the formalization)
  # nc: number of cores to use for rxns
  # detail: whether to return more details or only the MTA scores
  # return both the mta.model, and the summary data.table of MTA scores for the rxns, in a list(mta.model, result)

  # formulate MTA model
  mta.model <- form.mta(model, flux0, dflux, mta.pars)

  # solve the MTA models across the rxns
  if (!is.null(ko)) {
    ko <- all2idx(mta.model, ko)
    mta.model$lb[ko] <- 0
    mta.model$ub[ko] <- 0
  }
  res <- run.ko.screen(mta.model, rxns, run.mta, solv.pars=solv.pars, detail=detail, nc=nc)

  list(mta.model=mta.model, result=res)
}


### --- helper functions for the individual internal steps of MTA --- ###

form.mta <- function(model, flux0, dflux, mta.pars) {
  # formulate MTA model (the original MIQP version)
  # flux0 is the reference flux vector from sampling an iMAT output model
  # dflux is the flux change, i.e. output from de.dt2dflux()
  # pars is mta.pars, the parameters for MTA
  pars <- get.pars("mta", mta.pars)

  n.mets <- nrow(model$S)
  n.rxns <- ncol(model$S)

  # constraint matrix
  ## for reactions to be changed in the forward direction
  rxns.fw <- which(flux0>=0 & dflux==1 | flux0<0 & dflux==-1)
  n.fw <- length(rxns.fw)
  m1 <- sparseMatrix(1:n.fw, rxns.fw, dims=c(2*n.fw, n.rxns))
  m2 <- rbind(cbind(Diagonal(n.fw, x=(-flux0[rxns.fw]-pars$epsil)), Diagonal(n.fw, x=(-pars$v.min))),
              cbind(Diagonal(n.fw), Diagonal(n.fw)))
  S <- rbind(cbind(model$S, sparseMatrix(NULL, NULL, dims=c(n.mets, 2*n.fw))), cbind(m1, m2))
  ## for reactions to be changed in the backward direction
  rxns.bk <- which((flux0<=0 & dflux==1 | flux0>0 & dflux==-1) & !(flux0==0 & model$lb==0))
  n.bk <- length(rxns.bk)
  m1 <- sparseMatrix(1:n.bk, rxns.bk, dims=c(2*n.bk, n.rxns+2*n.fw))
  m2 <- rbind(cbind(Diagonal(n.bk, x=(-flux0[rxns.bk]+pars$epsil)), Diagonal(n.bk, x=(-pars$v.max))),
              cbind(Diagonal(n.bk), Diagonal(n.bk)))
  S <- rbind(cbind(S, sparseMatrix(NULL, NULL, dims=c(n.mets+2*n.fw, 2*n.bk))), cbind(m1, m2))

  # constraints
  rowlb <- c(model$rowlb, rep(0:1, each=n.fw), rep(c(pars$v.min.c,1), each=n.bk))
  rowub <- c(model$rowub, rep(c(pars$v.max.c,1), each=n.fw), rep(0:1, each=n.bk))
  lb <- c(model$lb, rep(0, 2*n.fw+2*n.bk))
  ub <- c(model$ub, rep(1, 2*n.fw+2*n.bk))

  # objective function and others
  n <- ncol(S)
  rxns.st <- which(dflux==0 & model$c!=1 & !grepl("biomass", model$rxns, ignore.case=TRUE))
  tmp <- rep(0, n)
  tmp[rxns.st] <- 2*(1-pars$alpha)
  Q <- .sparseDiagonal(x=tmp) # need to use .sparseDiagonal() instead of Diagonal() since Rcplex requires Q to be of the class dsparseMatrix
  c <- rep(c(0, pars$alpha/2, 0, pars$alpha/2), c(n.rxns+n.fw, n.fw, n.bk, n.bk))
  c[rxns.st] <- -2*(1-pars$alpha)*flux0[rxns.st]
  vtype <- rep(c("C","I"), c(n.rxns, n-n.rxns))

  # return MTA model
  list(flux0=flux0, dflux=dflux,
       rxns.fw=rxns.fw, rxns.bk=rxns.bk, rxns.st=rxns.st,
       rxns=model$rxns, mets=model$mets, csense="min", c=c, Q=Q, S=S, rowlb=rowlb, rowub=rowub, lb=lb, ub=ub, vtype=vtype)
}

get.mta.score <- function(model, miqp.out, detail) {
  # given one MTA model and the output from solving the model, compute the MTA score and other associated info
  # return a 1-row data.table

  miqp.out <- miqp.out[[1]]
  if (length(miqp.out$xopt)==1 && is.na(miqp.out$xopt)) {
    if (detail) {
      return(data.table(solv.stat=miqp.out$stat, v.opt=NA, rxns.change.yes=NA, rxns.change.no=NA, rxns.change.overdo=NA, advs.change.yes=NA, advs.change.no=NA, advs.change.overdo=NA, advs.steady=NA, score.change=NA, score.steady=NA, score.mta=NA))
    } else {
      return(data.table(solv.stat=miqp.out$stat, score.change=NA, score.steady=NA, score.mta=NA))
    }
  }

  v0 <- model$flux0
  n <- length(v0)
  v <- miqp.out$xopt[1:n]
  fw <- (1:n) %in% setdiff(model$rxns.fw, model$rxns.bk)
  bk <- (1:n) %in% setdiff(model$rxns.bk, model$rxns.fw)
  fw.or.bk <- (1:n) %in% intersect(model$rxns.bk, model$rxns.fw) # these are those reactions intended to change in either direction with flux0=0

  # reactions intended to change: overdone (thus regared as failed)
  fw.overdo <- fw & v0<0 & v>(-v0)
  bk.overdo <- bk & v0>0 & v<(-v0)
  # reactions intended to change: successful
  fw.yes <- fw & v>v0 & !fw.overdo
  bk.yes <- bk & v<v0 & !bk.overdo
  fw.or.bk.yes <- fw.or.bk & v!=0
  # reactions intended to change: failed
  fw.no <- fw & v<v0
  bk.no <- bk & v>v0
  fw.or.bk.no <- fw.or.bk & v==0

  # calculate mta score
  # reactions intended to change
  yes <- which(fw.yes|bk.yes|fw.or.bk.yes)
  no <- which(fw.no|bk.no|fw.or.bk.no)
  overdo <- which(fw.overdo|bk.overdo)
  adv.yes <- abs(v[yes]-v0[yes])
  adv.no <- abs(v[no]-v0[no])
  adv.overdo <- abs(abs(v[overdo])-abs(v0[overdo]))
  s.ch <- sum(adv.yes)-sum(adv.no)-sum(adv.overdo)
  # reactions intended to stay steady
  adv.st <- abs(v[model$rxns.st]-v0[model$rxns.st])
  s.st <- sum(adv.st)
  s <- s.ch/s.st
  # return
  if (detail) {
    res <- data.table(solv.stat=miqp.out$stat, v.opt=list(v), rxns.change.yes=list(yes), rxns.change.no=list(no), rxns.change.overdo=list(overdo), advs.change.yes=list(adv.yes), advs.change.no=list(adv.no), advs.change.overdo=list(adv.overdo), advs.steady=list(adv.st), score.change=s.ch, score.steady=s.st, score.mta=s)
  } else {
    res <- data.table(solv.stat=miqp.out$stat, score.change=s.ch, score.steady=s.st, score.mta=s)
  }
  res
}

run.mta <- function(model, solv.pars, detail) {
  # solve the MTA MIQP models for each rxn
  # return the summary of MTA scores for the rxns in a data.table
  solv.pars <- get.pars("mip", solv.pars)
  solv.pars$trace <- 0

  solv.out <- solve.model(model, pars=solv.pars)
  get.mta.score(model, solv.out, detail)
}


### --- rMTA --- ###

rmta <- function(model, flux0, dflux, rxns="all+ctrl", ko=NULL, nc=1L, detail=TRUE, k=100, mta.pars=list(), mip.pars=list(), qp.pars=list()) {
  # the function for running rMTA (with MTA of the original MIQP version); k is the rMTA-specific parameter
  # flux0 is the reference flux vector from sampling an iMAT output model
  # dflux is the flux change, i.e. output from de.dt2dflux(); I make it separate as usually we need to try different parameters in de.dt2dflux()
  # rxns are the indices or IDs or rxns to run MTA on, by default all rxns plus the control wild-type model; rxns="all" to run for all rxns w/o the ctrl; if using indices, 0 means ctrl; if using IDs, "ctrl", means ctrl
  # ko: NULL; or rxn indices or IDs to KO to combine with those in rxns -- these KO's will be added before screening for those in rxns, but after the metal.model has been formed using the original wildtype model (this is the reasonable way since an existent KO may change the formalization)
  # nc: number of cores to use for rxns
  # detail: whether to return more details or only the MTA scores
  # return both the mta.model, and the summary data.table of MTA scores for the rxns, in a list(mta.model, result.mta, result.moma, result.rmta)

  # formulate MTA model for either direction (dflux and -dflux)
  mta.model <- form.mta(model, flux0, dflux, mta.pars)
  mta.model0 <- form.mta(model, flux0, -dflux, mta.pars)

  if (!is.null(ko)) {
    ko <- all2idx(model, ko)
    mta.model$lb[ko] <- 0
    mta.model$ub[ko] <- 0
    mta.model0$lb[ko] <- 0
    mta.model0$ub[ko] <- 0
    model$lb[ko] <- 0
    model$ub[ko] <- 0
  }

  # solve the MTA models across the rxns for both models
  message("rmta(): Running MTA for dflux.")
  res1 <- run.ko.screen(mta.model, rxns, run.mta, solv.pars=mip.pars, detail=detail, nc=nc)
  message("rmta(): Running MTA for -dflux.")
  res0 <- run.ko.screen(mta.model0, rxns, run.mta, solv.pars=mip.pars, detail=FALSE, nc=nc)

  # MOMA for dflux
  message("rmta(): Running MOMA.")
  tmpf <- function(x) get.mta.score(model=mta.model, x, detail=detail)
  res.moma <- moma(model, rxns, nc, flux0, obj=tmpf, solv.pars=qp.pars)

  # rMTA score
  res <- data.table(id=res1$id, rxn=res1$rxn, bTS=res1$score.mta, wTS=res0$score.mta, mTS=res.moma$score.mta)
  res[, rTS:=ifelse(bTS>0 & mTS>0 & wTS<0, k*mTS*(bTS-wTS), mTS)]

  list(mta.model=mta.model, result.mta=res1, result.moma=res.moma, result.rmta=res)
}




