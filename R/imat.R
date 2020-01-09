###### the iMAT algorithm (the original version) ######


### --- prepare data for iMAT --- ###

exprs2int <- function(model, x, q.lo=0.25, q.hi=0.75, na2zero=TRUE) {
  # transform a gene expression data from continuous values to discretized integer (-1/0/1) values for iMAT
  # x: a gene-by-sample expression matrix with gene names, or a vector of gene expression named by gene names
  # if x is a matrix, then will take the median across samples
  # then select only those genes in the model, then discretize based on q.lo and q.hi
  # if na2zero missing genes (i.e. model genes that are not in the expression data) will be set to 0, otherwise they remain as NA's

  if (is.matrix(x)) x <- apply(x, 1, median, na.rm=TRUE)
  if (!is.numeric(x) || is.null(names(x))) stop("Input data not provided in the correct format. Should be a named gene-by-sample matrix or a named gene expression vector.")
  # discretize gene expression values
  x <- x[model$genes]
  na.idx <- is.na(x)
  message(sprintf("%d model genes not in the expression data,", sum(na.idx)))
  qlo <- quantile(x, q.lo, na.rm=TRUE)
  qhi <- quantile(x, q.hi, na.rm=TRUE)
  x <- ifelse(x<qlo, -1L, ifelse(x>qhi, 1L, 0L))
  if (na2zero) {
    x[na.idx] <- 0L
    message("these gene values are set to 0.")
  } else message("these genes have NA values.")
  unname(x)
}


### --- iMAT --- ###

imat <- function(model, expr, imat.pars=list(), solv.pars=list(), samp.pars=list()) {
  # the main function for running iMAT (the original version)
  # expr is the output from exprs2int(); I make it separate as sometimes we need to examine and modify expr before running imat
  # set samp.pars to NULL to skip the sampling
  # return both the imat.model, and the result model (the updated metabolic model), in a list(imat.model, result.model)

  # formulate iMAT model
  imat.model <- form.imat(model, expr, imat.pars)

  # solve the iMAT model
  imat.model <- run.imat(imat.model, imat.pars, solv.pars)

  # update the original metabolic model based on iMAT result
  res.model <- update.model.imat(model, imat.model, imat.pars)

  if (!is.null(samp.pars)) {
    # sample the result model
    res.model <- sample.model(res.model, samp.pars)
  }

  list(imat.model=imat.model, result.model=res.model)
}


### --- helper functions for the individual internal steps of iMAT --- ###

form.imat <- function(model, expr, imat.pars) {
  # formulate iMAT model (the original version)
  # expr is the output from exprs2int()
  # imat.pars: the parameters for iMAT
  pars <- get.pars("imat", imat.pars)

  # get intended activity level of rxns (-1/0/1) from discretized expression data
  rxns.int <- exprs2fluxes(model, expr)
  rxns.int[model$lb==0 & model$ub==0] <- 0

  n.mets <- nrow(model$S)
  n.rxns <- ncol(model$S)
  S <- model$S

  # 1. Active reactions: specify the y+ indicator variables, representing activation in the forward direction (i.e. v>flux.act)
  rxns.act <- which(rxns.int==1)
  n.act <- length(rxns.act)
  if (n.act!=0) {
    m1 <- sparseMatrix(1:n.act, rxns.act, dims=c(n.act, n.rxns))
    m2 <- Diagonal(n.act, x=(-pars$flux.act-pars$flux.bound))
    S <- rbind(cbind(S, sparseMatrix(NULL, NULL, dims=c(n.mets, n.act))), cbind(m1, m2))
  }

  # 2. Reversible active reactions: for those reversible ones among the active reactions, specify the extra y- indicator variables, representing activation in the backward direction (i.e. v<-flux.act)
  # thus, an reversible active reaction has both the y+ and y- indicator variables, because it can be active in either direction (but never both, i.e. 1 XOR 2)
  rxns.act.rev <- which(rxns.int==1 & model$lb<0)
  n.act.rev <- length(rxns.act.rev)
  if (n.act.rev!=0) {
    m1 <- sparseMatrix(1:n.act.rev, rxns.act.rev, dims=c(n.act.rev, ncol(S)))
    m2 <- Diagonal(n.act.rev, x=pars$flux.act+pars$flux.bound)
    S <- rbind(cbind(S, sparseMatrix(NULL, NULL, dims=c(nrow(S), n.act.rev))), cbind(m1, m2))
  }

  # 3. Inactive reactions: specify the y0 indicator variables
  # 3a. specify inactivation in the forward direction (i.e. v<flux.inact)
  rxns.inact <- which(rxns.int==-1)
  n.inact <- length(rxns.inact)
  if (n.inact!=0) {
    m1 <- sparseMatrix(1:n.inact, rxns.inact, dims=c(n.inact, ncol(S)))
    m2 <- Diagonal(n.inact, x=pars$flux.bound-pars$flux.inact)
    S <- rbind(cbind(S, sparseMatrix(NULL, NULL, dims=c(nrow(S), n.inact))), cbind(m1, m2))
  }
  # 3b. for those reversible inactive reactions, need to further specify inactivation in the backward direction (i.e. v>-flux.inact)
  # note that a reversible inactive reaction has only one y0 indicator variable, because for these reactions we want -flux.inact<v<flux.inact (3a AND 3b)
  rxns.inact.rev <- which(rxns.int==-1 & model$lb<0)
  n.inact.rev <- length(rxns.inact.rev)
  if (n.inact.rev!=0) {
    m3 <- sparseMatrix(1:n.inact.rev, rxns.inact.rev, dims=c(n.inact.rev, ncol(S)-n.inact))
    m4 <- sparseMatrix(1:n.inact.rev, match(rxns.inact.rev, rxns.inact), x=pars$flux.inact-pars$flux.bound, dims=c(n.inact.rev, n.inact))
    S <- rbind(S, cbind(m3, m4))
  }

  # other parameters
  n <- n.act + n.act.rev + n.inact + n.inact.rev
  rowlb <- c(model$rowlb, rep(-pars$flux.bound, n))
  rowub <- c(model$rowub, rep(pars$flux.bound, n))
  n <- ncol(S) - n.rxns
  c <- rep(c(0,1), c(n.rxns, n))
  vtype <- ifelse(c==1, "I", "C")
  lb <- c(model$lb, rep(0, n))
  ub <- c(model$ub, rep(1, n))
  var.ind <- rep(c("v","y+","y-","y0"), c(n.rxns, n.act, n.act.rev, n.inact)) # iMAT variable type indicators (v: fluxex; y+/-/0: indicator variables)

  # iMAT model
  imat.model <- list(exprs.int=expr, fluxes.int=rxns.int,
                    rxns.act=rxns.act, rxns.act.rev=rxns.act.rev, rxns.inact=rxns.inact, rxns.inact.rev=rxns.inact.rev, var.ind=var.ind,
                    rxns=model$rxns, mets=model$mets, csense="max", c=c, S=S, rowlb=rowlb, rowub=rowub, lb=lb, ub=ub, vtype=vtype)
  if ("irxn.ids" %in% names(model)) imat.model$irxn.ids <- model$irxn.ids # for multicellular imat models, keep the irxns.ids variable (intracellular reactions indices)
  imat.model
}

run.imat <- function(imat.model, imat.pars, solv.pars) {
  # solve the iMAT MILP model
  # return the imat.model with the results from running iMAT added:
  # imat.model$fluxes.int.imat: a vector in the order of the model rxns, with values 0/9/1/-1 representing a rxn being inactive, activity level not enforced, active in the forward direction, and active in the backward direction as determined by iMAT
  # imat.model$solver.out: if imat.pars$mode==0, then this is the output of MILP solver; if mode==1, then this is a list(obj) where obj is a matrix of size #rxns-by-3 containing the MILP objective values for each rxn being forced inactive/foward-active/backward-active
  imat.pars <- get.pars("imat", imat.pars)
  solv.pars <- get.pars("mip", solv.pars)

  if (imat.pars$mode==0) {
    # mode 0: solving the model only once to determine (de)activated reactions
    milp.out <- solve.model(imat.model, pars=solv.pars)
    xopt <- get.xopt(imat.model, milp.out, imat.pars)
    imat.model$fluxes.int.imat <- get.imat.opt.flux.int(imat.model, xopt)
    imat.model$solver.out <- milp.out
  } else if (imat.pars$mode==1) {
    # mode 1: solve the iMAT MILP under the forced inactivaton/activation of each rxn, determine (de)activated rxns based on the resulting objective values
    tmp <- imat.mode1(imat.model, imat.pars, solv.pars)
    imat.model$fluxes.int.imat <- tmp$flux.int.imat
    imat.model$solver.out <- tmp$solver.out
  } else if (imat.pars$mode==2) {
    # mode 2: constrain the original MILP objective at the optimal value, then solve for min/max flux of each rxn
    tmp <- imat.mode2(imat.model, imat.pars, solv.pars)
    imat.model$fluxes.int.imat <- tmp$flux.int.imat
    imat.model$solver.out <- tmp$solver.out
  }
  imat.model
}

get.xopt <- function(imat.model, milp.out, pars) {
  # a helper function to extract a certain solution or a concensus solution by majority vote from milp.out
  # two parameters called sol and sol.major.cutoff should be in pars, sol specifies which solution to use, if 0 then pool together all solutions and get a consensus solution, that if >= sol.major.cutoff fraction of integer solution is 1, then set the consensus to 1, otherwise set the consensus to 0

  if (pars$sol==0) {
    xopt <- do.call(cbind, lapply(milp.out, function(x) {
      if (x$stat %in% c(101,102,128,129,130)) x$xopt else NULL
    }))
    if (is.null(xopt)) stop("Potential issues in all MILP solutions.")
    xopt <- rowMeans(xopt)
    xopt[imat.model$vtype=="I"] <- as.numeric(xopt[imat.model$vtype=="I"]>=pars$sol.major.cutoff)
  } else {
    if (!milp.out[[pars$sol]]$stat %in% c(101,102,128,129,130)) stop("Potential issues in the selected MILP solution.")
    xopt <- milp.out[[pars$sol]]$xopt
  }
  xopt
}

get.imat.opt.flux.int <- function(imat.model, xopt) {
  # given the formulated imat.model, and the xopt of solving that model, return a vector in the order of the model rxns, with values 0/9/1/-1 representing a rxn being inactive, activity level not enforced, active in the forward direction, and active in the backward direction as determined by iMAT

  yp <- xopt[imat.model$var.ind %in% c("y+","y+_1","y+_2")]
  ym <- xopt[imat.model$var.ind %in% c("y-","y-_1","y-_2")]
  y0 <- xopt[imat.model$var.ind %in% c("y0","y0_1","y0_2")]
  fw <- imat.model$rxns.act[yp==1]
  bk <- imat.model$rxns.act.rev[ym==1]
  inact <- imat.model$rxns.inact[y0==1]
  yres <- rep(9L, sum(imat.model$var.ind %in% c("v","v_1","v_2")))
  yres[fw] <- 1L
  yres[bk] <- -1L
  yres[inact] <- 0L
  yres
}

imat.mode1 <- function(imat.model, imat.pars, solv.pars) {
  # the mode 1 of imat (the fva-like approach): solve the iMAT MILP under the forced inactivaton/activation of each rxn, determine (de)activated rxns based on the resulting objective values
  # return a list(solver.out, flux.int.imat), solver.out is a data.table of the optimal iMAT objectives for all rxns, flux.int.imat is a vector in the order of the model rxns, with values 0/9/1/-1 representing a rxn being inactive, activity level not enforced, active in the forward direction, and active in the backward direction as determined by iMAT

  solv.pars$nsol <- 1
  solv.pars$tilim <- 120
  rxns <- which(imat.model$var.ind=="v")
  names(rxns) <- imat.model$rxns[rxns]
  obj <- rbindlist(parallel::mclapply(rxns, function(i) {
    lb <- imat.model$lb
    ub <- imat.model$ub
    lbi <- lb[i]
    ubi <- ub[i]
    # force inactivation
    if (lbi<0) lb[i] <- -imat.pars$flux.inact
    ub[i] <- imat.pars$flux.inact
    obj.inact <- solve.model(imat.model, lb=lb, ub=ub, pars=solv.pars)[[1]]$obj
    # force activation in the forward direction
    lb[i] <- imat.pars$flux.act
    ub[i] <- ubi
    obj.fw <- solve.model(imat.model, lb=lb, ub=ub, pars=solv.pars)[[1]]$obj
    # force activation in the backward direction
    if (lbi<0) {
      lb[i] <- lbi
      ub[i] <- -imat.pars$flux.act
      obj.bk <- solve.model(imat.model, lb=lb, ub=ub, pars=solv.pars)[[1]]$obj
    } else obj.bk <- -1
    data.table(inact=obj.inact, act.fw=obj.fw, act.bk=obj.bk)
  }, mc.cores=imat.pars$nc), idcol="rxn")
  obj <- cbind(data.table(id=rxns), obj)
  tmp <- apply(objs[,-1:-2], 1, function(x) {
    if (any(is.na(x))) {
      return(0)
    } else if (uniqueN(x)==3) {
      return(which.max(x))
    } else if (x[2]==x[3] && x[1]>x[2]) {
      return(1)
    } else if (x[1]==x[3] && x[2]>x[1]) {
      return(2)
    } else if (x[1]==x[2] && x[3]>x[1]) {
      return(3)
    } else return(0)
  })
  fint <- ifelse(tmp==1, 0L, ifelse(tmp==2, 1L, ifelse(tmp==3, -1L, 9L)))
  # result
  list(solver.out=obj, flux.int.imat=fint)
}

imat.mode2 <- function(imat.model, imat.pars, solv.pars) {
  # the mode 2 of imat: similar to mode 1, but after solving the original iMAT MILP, constrain the objective at the optimal value, then for each reaction solve a pair of MILP's to obtain its min/max fluxes
  # return a list(solver.out, flux.int.imat), solver.out is a data.table of the min/max fluxes for all rxns, flux.int.imat is a vector in the order of the model rxns, with values 0/9/1/-1 representing a rxn being inactive, activity level not enforced, active in the forward direction, and active in the backward direction as determined by iMAT

  solv.pars$nsol <- 1
  solv.pars$tilim <- 120
  obj.opt <- solve.model(imat.model, pars=solv.pars)[[1]]$obj
  imat.model.opt <- add.constraint(imat.model, 1:length(imat.model$c), imat.model$c, obj.opt, obj.opt)
  fva.res <- fva(imat.model.opt, rxns="all", nc=imat.pars$nc, solv.pars=solv.pars)
  fint <- ifelse(fva.res$vmin>=imat.pars$flux.act, 1L, ifelse(fva.res$vmax<= -imat.pars$flux.act, -1L, ifelse(fva.res$vmax<=imat.pars$flux.inact & fva.res$vmin>= -imat.pars$flux.inact, 0L, 9L)))
  # result
  list(solver.out=fva.res, flux.int.imat=fint)
}

update.model.imat <- function(model, imat.res, imat.pars) {
  # update metabolic model (model) based on the result of iMAT (imat.res), using parameters given in imat.pars
  # return the updated model
  pars <- get.pars("imat", imat.pars)

  if (pars$mode==2) {
    tmp <- !is.na(imat.res$solver.out$vmin)
    model$lb[tmp] <- imat.res$solver.out$vmin[tmp]
    tmp <- !is.na(imat.res$solver.out$vmax)
    model$ub[tmp] <- imat.res$solver.out$vmax[tmp]
  } else {
    fint <- imat.res$fluxes.int.imat
    model$lb[fint==1] <- pars$flux.act
    model$ub[fint==-1] <- -pars$flux.act
    model$ub[fint==0] <- pars$flux.inact
    model$lb[fint==0 & model$lb<0] <- -pars$flux.inact
  }
  model
}


