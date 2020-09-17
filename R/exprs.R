###### additional algorithms for incorporating expression data ######


### --- GIMME and E-Fmin --- ### 

gimme <- function(model, expr, rmfs="biomass", lbs=0.01, relative=FALSE, norm=TRUE, cutoff=1, mode=c(0,1), gap=NULL, agap=NULL, solv.pars=list(), samp.pars=NULL) {
  # function for running GIMME or E-Fmin; by default will run E-Fmin
  # expr is a vector of gene expression values (may be pre-normalized), either with names as gene symbols as in model$genes, or if unnamed the length and order should correspond to model$genes
  # rmfs: reaction indices or IDs (as in model$rxns) for "required metabolic functions"; or regex for biomass reaction; or NULL (if want to manually add constraints wrt rmf prior to calling this function)
  # lbs: minimal fluxes for rmfs in the corresponding order, or a single value; if relative=TRUE, lbs will be treated as the fraction of the maximal flux values; the default correspond to what's used in the E-Fmin paper
  # norm: whether to scale and shift the reaction expression so that they are within [0,1], default TRUE corresponds to E-Fmin
  # cutoff: gene expression cutoff above which the weight is 0 for GIMME; default 1 corresponds to E-Fmin
  # mode 0 means solve the model once and return list(gimme.model, fluxes) where fluxes is the vector of optimal fluxes;
  # mode 1 means first get optimal objective value, then constrain obj at optimal or near optimal (as defined by gap and agap), then run FVA, and return list(gimme.model, result.model, fluxes) where result.model is updated model with rxn bounds set to the FVA output
  # when mode==1 and samp.pars not NULL, will sample the resulting model
  
  mode <- match.arg(as.character(mode[1]), c("0","1"))
  
  # process expression values
  x <- exprs2fluxes(model, expr, na.after=NA)
  if (norm) {
  	xmin <- min(x, na.rm=TRUE)
  	xmax <- max(x, na.rm=TRUE)
  	if (xmin<0) x <- (x-xmin)/xmax else x <- x/xmax # original E-Fmin simply used x/xmax, seeming to assume that all expression values are positive and zero means no expression
  }
  w <- ifelse(x<cutoff, cutoff-x, 0)
  w[is.na(w)] <- 0

  # set the lower bounds of required reactions; for the reversible reactions, the constraints of |v|>lbs need to be handled with interger programming later
  if (!is.null(rmfs)) {
  	model <- set.required.rxns(model, rmfs, lbs, relative=relative)
  	rxns <- model$no.set.rxns
  	lbs <- model$no.set.lbs
  	if ("no.set.rxns" %in% names(model)) model <- model$model
  } else rxns <- lbs <- NULL

  # form model
  gimme.model <- form.gimme(model, w, rxns, lbs)

  # solve model and extract results
  if (any(isrev)) solv.pars <- get.pars("mip", solv.pars) else solv.pars <- get.pars("lp", solv.pars)
  if (mode=="0") {
  	solv.res <- solve.model(gimme.model, pars=solv.pars)[[1]]
  	gimme.model$solver.out <- solv.res
  	v <- solv.res$xopt[1:length(gimme.model$rxns)]
  	res <- list(gimme.model=gimme.model, fluxes=v)
  } else if (mode=="1") {
  	solv.res <- fva1(gimme.model, gap=gap, agap=agap, keep.solv.out=TRUE, solv.pars=solv.pars)
  	gimme.model$solver.out <- solv.res$solver.out
  	v <- solv.res$solver.out$xopt[1:length(gimme.model$rxns)]
  	res.model <- model
  	res.model$lb <- solv.res$fva.res$vmin
  	res.model$ub <- solv.res$fva.res$vmax
  	if (!is.null(samp.pars)) {
  	  res.model <- sample.model(res.model, samp.pars)
  	}
  	res <- list(gimme.model=gimme.model, result.model=res.model, fluxes=v)
  }
  
  res
}

form.gimme <- function(model, w, rmfs, rmf.lbs) {
  # formulate a GIMME-like model

  n.mets <- nrow(model$S)
  n.rxns <- ncol(model$S)

  # formulating the minimization of weighted sum of absolute values \Sum w|v|, by replacing v with two slack variables v1, v2 such that v=v1-v2, then |v| can be expressed as v1+v2
  # however, I still keep v in the formulation because it makes it easier to run FVA after initially solving the model, if needed
  rxns1 <- which(w!=0 & !is.na(w))
  w <- w[rxns1]
  n <- length(rxns1)
  S <- rbind(
    cbind( model$S,                          sparseMatrix(NULL, NULL, dims=c(n.mets, 2*n)) ),
    cbind( sparseMatrix(1:n, rxns1, x=1, dims=c(n, n.rxns)),  Diagonal(n,-1),  Diagonal(n) )
  )
  rowlb <- c(model$rowlb, rep(0, n))
  rowub <- c(model$rowub, rep(0, n))
  lb <- c(model$lb, rep(0, n))
  ub <- c(model$ub, rep(Inf, n))
  c <- c(rep(0, n.rxns), w, w)
  vtype <- rep("C", length(c))

  # handling |v|>lb constraints, if any, by adding an binary variable y and the constraint lb <= v + My <= M-lb, where M is a large constant, here I used M=1e4
  n.rmfs <- length(rmfs)
  if (n.rmfs>0) {
	S <- rbind(
	  cbind( S,                            sparseMatrix(NULL, NULL, dims=c(nrow(S), n.rmfs)) )
	  cbind( sparseMatrix(1:n.rmfs, rmfs, dims=c(n.rmfs, ncol(S))), Diagonal(n.rmfs, x=1e4) )
  	)
  	rowlb <- c(rowlb, rmf.lbs)
  	rowub <- c(rowub, 1e4-rmf.lbs)
  	lb <- c(lb, rep(0, n.rmfs))
  	ub <- c(ub, rep(1, n.rmfs))
  	c <- c(c, rep(0, n.rmfs))
  	vtype <- c(vtype, rep("I", n.rmfs))
  }

  # return model
  list(rxns=model$rxns, mets=model$mets, csense="min", c=c, S=S, rowlb=rowlb, rowub=rowub, lb=lb, ub=ub, vtype=vtype)
}


### --- Lee et al. BMC Syst Biol 2012 --- ### 

lee.smallbone <- function(model, expr, rmfs="biomass", lbs=0.01, relative=FALSE, norm=TRUE, cutoff=1, mode=c(0,1), gap=NULL, agap=NULL, solv.pars=list(), samp.pars=NULL) {
  # function for running the algorithm of Lee et al., the function name is from the surnames of the two co-first authors
  # expr is a vector of gene expression values (may be pre-normalized), either with names as gene symbols as in model$genes, or if unnamed the length and order should correspond to model$genes
  # rmfs: reaction indices or IDs (as in model$rxns) for "required metabolic functions"; or regex for biomass reaction; or NULL (if want to manually add constraints wrt rmf prior to calling this function)
  # lbs: minimal fluxes for rmfs in the corresponding order, or a single value; if relative=TRUE, lbs will be treated as the fraction of the maximal flux values; the default correspond to what's used in the E-Fmin paper
  # norm: whether to scale and shift the reaction expression so that they are within [0,1], default TRUE corresponds to E-Fmin
  # cutoff: gene expression cutoff above which the weight is 0 for GIMME; default 1 corresponds to E-Fmin
  # mode 0 means solve the model once and return list(gimme.model, fluxes) where fluxes is the vector of optimal fluxes;
  # mode 1 means first get optimal objective value, then constrain obj at optimal or near optimal (as defined by gap and agap), then run FVA, and return list(gimme.model, result.model, fluxes) where result.model is updated model with rxn bounds set to the FVA output
  # when mode==1 and samp.pars not NULL, will sample the resulting model
  
  mode <- match.arg(as.character(mode[1]), c("0","1"))
  if (length(rmfs)==1 && grepl("biomass", rmfs, ignore.case=TRUE)) {
    rxns <- get.biomass.idx(model, rgx=rmfs)
  } else rxns <- all2idx(model, rmfs)

  # process expression values
  x <- exprs2fluxes(model, expr, na.after=NA)
  if (norm) {
  	xmin <- min(x, na.rm=TRUE)
  	xmax <- max(x, na.rm=TRUE)
  	if (xmin<0) x <- (x-xmin)/xmax else x <- x/xmax # original E-Fmin simply used x/xmax, seeming to assume that all expression values are positive and zero means no expression
  }
  w <- ifelse(x<cutoff, cutoff-x, 0)
  w[is.na(w)] <- 0

  # set the lower bounds of required reactions; for the reversible reactions, the constraints of |v|>lbs need to be handled with interger programming later
  if (!is.null(rmfs)) {
  	if (length(lbs)==1) lbs <- rep(lbs, length(rxns))
    if (relative) lbs <- lbs * get.opt.fluxes(model, rxns, "max")
    isrev <- model$lb[rxns] < -lbs
    model$lb[rxns[!isrev]] <- lbs[!isrev]
    rxns <- rxns[isrev]
    lbs <- lbs[isrev]
  } else rxns <- lbs <- NULL

  # form model
  gimme.model <- form.gimme(model, w, rxns, lbs)

  # solve model and extract results
  if (any(isrev)) solv.pars <- get.pars("mip", solv.pars) else solv.pars <- get.pars("lp", solv.pars)
  if (mode=="0") {
  	solv.res <- solve.model(gimme.model, pars=solv.pars)[[1]]
  	gimme.model$solver.out <- solv.res
  	v <- solv.res$xopt[1:length(gimme.model$rxns)]
  	res <- list(gimme.model=gimme.model, fluxes=v)
  } else if (mode=="1") {
  	solv.res <- fva2(gimme.model, gap=gap, agap=agap, keep.solv.out=TRUE, solv.pars=solv.pars)
  	gimme.model$solver.out <- solv.res$solver.out
  	v <- solv.res$solver.out$xopt[1:length(gimme.model$rxns)]
  	res.model <- model
  	res.model$lb <- solv.res$fva.res$vmin
  	res.model$ub <- solv.res$fva.res$vmax
  	if (!is.null(samp.pars)) {
  	  res.model <- sample.model(res.model, samp.pars)
  	}
  	res <- list(gimme.model=gimme.model, result.model=res.model, fluxes=v)
  }
  
  res
}

form.gimme <- function(model, w, rmfs, rmf.lbs) {
  # formulate a GIMME-like model

  n.mets <- nrow(model$S)
  n.rxns <- ncol(model$S)

  # formulating the minimization of weighted sum of absolute values \Sum w|v|, by replacing v with two slack variables v1, v2 such that v=v1-v2, then |v| can be expressed as v1+v2
  # however, I still keep v in the formulation because it makes it easier to run FVA after initially solving the model, if needed
  rxns1 <- which(w!=0 & !is.na(w))
  w <- w[rxns1]
  n <- length(rxns1)
  S <- rbind(
    cbind( model$S,                          sparseMatrix(NULL, NULL, dims=c(n.mets, 2*n)) ),
    cbind( sparseMatrix(1:n, rxns1, x=1, dims=c(n, n.rxns)),  Diagonal(n,-1),  Diagonal(n) )
  )
  rowlb <- c(model$rowlb, rep(0, n))
  rowub <- c(model$rowub, rep(0, n))
  lb <- c(model$lb, rep(0, n))
  ub <- c(model$ub, rep(Inf, n))
  c <- c(rep(0, n.rxns), w, w)
  vtype <- rep("C", length(c))

  # handling |v|>lb constraints, if any, by adding an binary variable y and the constraint lb <= v + My <= M-lb, where M is a large constant, here I used M=1e4
  n.rmfs <- length(rmfs)
  if (n.rmfs>0) {
	S <- rbind(
	  cbind( S,                            sparseMatrix(NULL, NULL, dims=c(nrow(S), n.rmfs)) )
	  cbind( sparseMatrix(1:n.rmfs, rmfs, dims=c(n.rmfs, ncol(S))), Diagonal(n.rmfs, x=1e4) )
  	)
  	rowlb <- c(rowlb, rmf.lbs)
  	rowub <- c(rowub, 1e4-rmf.lbs)
  	lb <- c(lb, rep(0, n.rmfs))
  	ub <- c(ub, rep(1, n.rmfs))
  	c <- c(c, rep(0, n.rmfs))
  	vtype <- c(vtype, rep("I", n.rmfs))
  }

  # return model
  list(rxns=model$rxns, mets=model$mets, csense="min", c=c, S=S, rowlb=rowlb, rowub=rowub, lb=lb, ub=ub, vtype=vtype)
}

