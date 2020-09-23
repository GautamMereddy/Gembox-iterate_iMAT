###### additional algorithms for incorporating expression data ######


### --- GIMME and E-Fmin --- ### 

gimme <- function(model, expr, rmfs="biomass", lbs=0.01, relative=FALSE, norm=TRUE, cutoff=1, mode=c(0,1), nc=1L, gap=NULL, agap=NULL, solv.pars=list(), samp.pars=NULL) {
  # function for running GIMME or E-Fmin; by default will run E-Fmin
  # expr is a vector of gene expression values (may be pre-normalized), either with names as gene symbols as in model$genes, or if unnamed the length and order should correspond to model$genes
  # rmfs: reaction indices or IDs (as in model$rxns) for "required metabolic functions"; or regex for biomass reaction; or NULL (if want to manually add constraints wrt rmf prior to calling this function)
  # lbs: minimal fluxes for rmfs in the corresponding order, or a single value; if relative=TRUE, lbs will be treated as the fraction of the maximal flux values; the default correspond to what's used in the E-Fmin paper
  # norm: whether to scale and shift the reaction expression so that they are within [0,1], default TRUE corresponds to E-Fmin
  # cutoff: gene expression cutoff above which the weight is 0 for GIMME; default 1 corresponds to E-Fmin
  # mode 0 means solve the model once and return list(gimme.model, fluxes) where fluxes is the vector of optimal fluxes;
  # mode 1 means first get optimal objective value, then constrain obj at optimal or near optimal (as defined by gap and agap), then run FVA with n.cores=nc, and return list(gimme.model, result.model, fluxes) where result.model is updated model with rxn bounds set to the FVA output
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
  	if ("no.set.rxns" %in% names(model)) {
  	  model <- model$model
  	  solv.pars <- get.pars("mip", solv.pars)
  	} else solv.pars <- get.pars("lp", solv.pars)
  } else rxns <- lbs <- NULL

  # form model
  gimme.model <- form.gimme(model, w, rxns, lbs)

  # solve model and extract results
  if (mode=="0") {
  	message("Solving model...")
  	solv.res <- solve.model(gimme.model, pars=solv.pars)[[1]]
  	gimme.model$solver.out <- solv.res
  	v <- solv.res$xopt[1:length(gimme.model$rxns)]
  	res <- list(gimme.model=gimme.model, fluxes=v)
  } else if (mode=="1") {
  	message("Mode 1. Running FVA...")
  	solv.res <- fva1(gimme.model, nc=nc, gap=gap, agap=agap, keep.solv.out=TRUE, solv.pars=solv.pars)
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
  
  message("Done.")
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
  lb <- c(model$lb, rep(0, 2*n))
  ub <- c(model$ub, rep(Inf, 2*n))
  c <- c(rep(0, n.rxns), w, w)
  vtype <- rep("C", length(c))

  # handling |v|>lb constraints, if any, by adding an binary variable y and the constraint lb <= v + My <= M-lb, where M is a large constant, here I used M=1e4
  n.rmfs <- length(rmfs)
  if (n.rmfs>0) {
	S <- rbind(
	  cbind( S,                            sparseMatrix(NULL, NULL, dims=c(nrow(S), n.rmfs)) ),
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

exprs2fluxes.s <- function(model, x, s, x.na.before=NA, s.na.before=NA, x.na.after=NA, s.na.after=NA) {
  # map a numeric vector x of expression levels of genes to the flux levels of reactions in the model, together with a s vector of other values associated with each gene, e.g. sd of each gene, or some other score
  # x should contain either continuous expression values (then the output will also be continuous) or discrete values from {-1,0,1,NA} representing low/medium/high expression levels
  # x should be named by gene symbols as used in model$genes, or if it's unnamed and length being length(model$genes), assume it's already in the same order as model$genes
  # s should correspond to x in the same order
  # NA's in x will be kept and propagated during the conversion to flux values by default; or specify na.before and NA will be changed to this value
  # NA's in the result flux values will be changed to 0 by default; or specify na.after and NA will be changed to that value (or na.after=NA to keep NA)

  if (is.null(names(x))) {
    if (length(x)==length(model$genes)) {
      message("exprs2fluxes.s(): Assuming the input vector is in the same order as model$genes.")
    } else stop("Input vector and model$genes have different lengths!")
  } else {
  	tmp <- match(model$genes, names(x))
    x <- x[tmp]
    s <- s[tmp]
    if (all(is.na(x))) stop("Input doesn't contain any of the model genes!")
  }
  x[is.na(x)] <- x.na.before
  s[is.na(s)] <- s.na.before

  `&` <- function(a,b) min(a,b)
  `|` <- function(a,b) max(a,b)

  res <- sapply(model$rules, function(i) {
  	if (i=="") {
  	  res <- c(NA, NA)
  	} else {
  	  r1 <- eval(parse(text=i))
  	  if (!is.na(r1)) {
  	  	ids <- as.integer(stringr::str_extract_all(i, "[0-9]+")[[1]])
  	  	r2 <- s[ids][match(r1, x[ids])]
  	  } else r2 <- NA
  	  res <- c(r1,r2)
  	}
  	res
  })
  x.res <- res[1,]
  s.res <- res[2,]
  x.res[is.na(x.res)] <- x.na.after
  s.res[is.na(s.res)] <- s.na.after
  x.res[model$rules==""] <- NA # still NA for rxns w/o genes
  s.res[model$rules==""] <- NA
  list(x=unname(x.res), s=unname(s.res))
}

lee.smallbone <- function(model, expr, sd=1, mode=c(0,1), nc=1L, gap=NULL, agap=NULL, solv.pars=list(), samp.pars=NULL) {
  # function for running the algorithm of Lee et al., the function name is from the surnames of the two co-first authors
  # expr is a vector of non-negative gene expression values (may be pre-normalized), either with names as gene symbols as in model$genes, or if unnamed the length and order should correspond to model$genes
  # sd is the standard deviation of gene expression corresonding to expr in the same order, used to generate weights
  # mode 0 means solve the model once and return list(ls.model, fluxes) where fluxes is the vector of optimal fluxes;
  # mode 1 means first get optimal objective value, then constrain obj at optimal or near optimal (as defined by gap and agap), then run FVA (with n.cores=nc), and return list(ls.model, result.model, fluxes) where result.model is updated model with rxn bounds set to the FVA output
  # when mode==1 and samp.pars not NULL, will sample the resulting model
  
  mode <- match.arg(as.character(mode[1]), c("0","1"))
  solv.pars <- get.pars("lp", solv.pars)

  # process expression values
  if (length(sd)==1) sd <- rep(sd, length(expr)) else if (length(sd)!=length(expr)) stop("sd should correspond to expr in the same order.")
  tmp <- exprs2fluxes.s(model, expr, sd)
  x <- tmp$x
  w <- 1/tmp$s

  # form model and solve iteratively
  nr <- which(model$lb>=0)
  message("Solving model iteratively...")
  i <- 1L
  repeat {
  	message("Iteration # ", i)
  	ls.model <- form.lee.smallbone(model, nr, x[nr], w[nr])
  	solv.res <- fva1(ls.model, rxns=setdiff(1:length(model$rxns), nr), nc=nc, keep.solv.out=TRUE, solv.pars=solv.pars)
  	nr.new <- solv.res$fva.res[vmin>=0,id]
  	if (length(nr.new)==0) break else nr <- c(nr, nr.new)
  	i <- i+1
  }
  
  # solve model and extract results
  if (mode=="0") {
  	ls.model$solver.out <- solv.res$solver.out
  	v <- solv.res$solver.out$xopt[1:length(ls.model$rxns)]
  	res <- list(ls.model=ls.model, fluxes=v)
  } else if (mode=="1") {
  	message("Mode 1. Running FVA...")
  	solv.res <- fva1(ls.model, nc=nc, gap=gap, agap=agap, keep.solv.out=TRUE, solv.pars=solv.pars)
  	ls.model$solver.out <- solv.res$solver.out
  	v <- solv.res$solver.out$xopt[1:length(ls.model$rxns)]
  	res.model <- model
  	res.model$lb <- solv.res$fva.res$vmin
  	res.model$ub <- solv.res$fva.res$vmax
  	if (!is.null(samp.pars)) {
  	  res.model <- sample.model(res.model, samp.pars)
  	}
  	res <- list(ls.model=ls.model, result.model=res.model, fluxes=v)
  }

  message("Done.")
  res
}

form.lee.smallbone <- function(model, idx, x, w) {
  # formulate the model of Lee et al.

  n.mets <- nrow(model$S)
  n.rxns <- ncol(model$S)

  # formulating the minimization of weighted sum of absolute values \Sum w|v-x|, by adding two slack variables v1, v2 such that v-x=v1-v2, then |v-x| can be expressed as v1+v2
  tmp <- !is.na(x) & !is.na(w)
  idx <- idx[tmp]
  x <- x[tmp]
  w <- w[tmp]
  n <- length(idx)
  S <- rbind(
    cbind( model$S,                        sparseMatrix(NULL, NULL, dims=c(n.mets, 2*n)) ),
    cbind( sparseMatrix(1:n, idx, x=1, dims=c(n, n.rxns)),  Diagonal(n,-1),  Diagonal(n) )
  )
  rowlb <- c(model$rowlb, x)
  rowub <- c(model$rowub, x)
  lb <- c(model$lb, rep(0, 2*n))
  ub <- c(model$ub, rep(Inf, 2*n))
  c <- c(rep(0, n.rxns), w, w)

  # return model
  list(rxns=model$rxns, mets=model$mets, csense="min", c=c, S=S, rowlb=rowlb, rowub=rowub, lb=lb, ub=ub)
}

