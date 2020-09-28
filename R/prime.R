###### the PRIME algorithm (Yizhak et al. eLife 2014) ######


prime <- function(model, expr, gr, padj.cutoff=0.05, nc=1L, bm.rgx="biomass", default.bnd=1e3, bm.epsil=1e-4, min.step.size=0.1, bm.lb.rel=0.1, solv.pars=list()) {
  # function to run PRIME all steps at once
  # expr: a gene-by-sample matrix, need to has rownames of gene symbols as those in model$genes (but doesn't need to be of same length or in order)
  # gr: cell growth rate measures across samples, corresponding to the order or columns of expr
  # padj.cutoff: used to define significant growth-associated rxns
  # default.bnd, bm.epsil, min.step.size passed to shrink.model.bounds(); bm.lb.res passed to get.norm.range.min()
  # return a list of updated models, one for each sample in the corresponding order

  solv.pars <- get.pars("lp", solv.pars)
  message("Preparing base model...")
  model <- prepare.model.prime(model, default.bnd=default.bnd, bm.epsil=bm.epsil, min.step.size=min.step.size, bm.rgx=bm.rgx, solv.pars=solv.pars)
  message("Identifying growth-associated reactions...")
  prm.rxns <- get.prime.rxns(model, expr=expr, gr=gr, nc=nc, padj.cutoff=padj.cutoff)
  res <- run.prime(model, prm.rxns=prm.rxns, gr=gr, nc=nc, bm.rgx=bm.rgx, bm.lb.rel=bm.lb.rel, solv.pars=solv.pars)
}

prepare.model.prime <- function(model, default.bnd=1e3, bm.epsil=1e-4, min.step.size=0.1, bm.rgx="biomass", solv.pars=get.pars("lp", list())) {
  # step 1 of PRIME: prepare base model
  # default.bnd, bm.epsil, min.step.size, bm.rgx passed to shrink.model.bounds()

  # split each reversible rxn into two forward/backward rxns
  model <- convert.rev.rxns(model)
  # shrink rxn bounds
  model <- shrink.model.bounds(model, rxns="default", default.bnd=default.bnd, bm.epsil=bm.epsil, relative=FALSE, min.step.size=min.step.size, bm.rgx=bm.rgx, solv.pars=solv.pars)
}

get.prime.rxns <- function(model, expr, gr, nc=1L, padj.cutoff=0.05) {
  # step 2 of PRIME: select growth-associated rxns (i.e. rxns with significant correlation with the provided growth rate data across multiple samples)
  # expr: a gene-by-sample matrix, need to has rownames of gene symbols as those in model$genes (but doesn't need to be of same length or in order)
  # gr: cell growth rate measures across samples, corresponding to the order or columns of expr
  # return list(x, cor, i), x is the all-rxn-by-sample matrix mapped from expr, cor is a data.table containing results of correlation tests, i is the indices of the resulting growth-associated rxns with padj<padj.cutoff
  # this is a separate step and the results are returned as such to allow manual adjustment of selected reactions w/o recomputing the entire correlation
  # manually adjust result with something like: res$i <- res$cor[padj<new.cutoff, id], then pass result to the next step

  expr <- expr[match(model$genes, rownames(expr)), ]
  if (all(is.na(expr))) stop ("expr doesn't contain any of the model genes!")

  # map gene expr values to rxns values by taking the mean
  mat <- t(sapply(model$rules, function(x) {
  	i <- as.integer(stringr::str_extract_all(x, "[0-9]+")[[1]])
  	colMeans(expr[i,,drop=FALSE], na.rm=TRUE)
  }))
  mat[is.nan(mat)] <- NA
  if (!is.null(colnames(expr))) colnames(mat) <- colnames(expr)

  # correlation between rxn values and growth rates across samples for each rxn
  cor.res <- rbindlist(parallel::mclapply(1:nrow(mat), function(i) {
  	tryCatch({
  	  a <- cor.test(mat[i,], gr, method="spearman")
  	  data.table(rho=a$estimate, pval=a$p.value)
  	}, error=function(e) data.table(rho=NA, pval=NA))
  }, mc.cores=nc), idcol="id")
  cor.res[, padj:=p.adjust(pval,"BH")]
  tmp <- sum(cor.res$padj<padj.cutoff, na.rm=TRUE)
  if (tmp==0) stop("No significant growth-associated reactions with padj<", padj.cutoff, ".") else message("Found ", tmp, " growth-associated reactions with padj<", padj.cutoff, ".")

  list(x=mat, cor=cor.res, i=cor.res[padj<padj.cutoff, id])
}

run.prime <- function(model, prm.rxns, gr, nc=1L, bm.rgx="biomass", bm.lb.rel=0.1, solv.pars=get.pars("lp", list())) {
  # step 3 of PRIME (all the remaining steps)
  # prm.rxns: output from get.prime.rxns
  # gr: cell growth rate measures across samples, corresponding to the order or columns of mat
  # bm.lb.res passed to get.norm.range.min()
  # return a list of updated models, one for each sample in the corresponding order

  # compute the lb and ub of the normalization range
  rmin <- get.norm.range.min(model, bm.lb.rel=bm.lb.rel, nc=nc, bm.rgx=bm.rgx, solv.pars=solv.pars)
  message("Computing the upper bound of normalization range...")
  rmax <- get.norm.range.max(model, rxns=prm.rxns$i, range.min=rmin, nc=nc, bm.rgx=bm.rgx, solv.pars=solv.pars)
  # compute and set final rxn ub values for each sample
  message("Generating output models...")
  mat <- prm.rxns$x[prm.rxns$i,] * prm.rxns$cor[match(prm.rxns$i, id), sign(rho)]
  mat <- do.call(rbind, parallel::mclapply(1:nrow(mat), function(i) {
  	x <- mat[i,]
  	rng <- range(x, na.rm=TRUE)
  	(x-rng[1]) / diff(rng)
  }, mc.cores=nc))
  ubs <- mat * (rmax-rmin) + rmin
  models <- apply(ubs, 2, function(x) {
  	m <- model
  	m$ub[prm.rxns$i] <- x
  	m
  })
  # recombine splitted reversible rxns
  res <- lapply(models, revert.rev.rxns)
  message("Done.")
  res
}

get.norm.range.min <- function(model, bm.lb.rel=0.1, nc=1L, bm.rgx="biomass", solv.pars=get.pars("lp", list())) {
  # get the lower bound value for the PRIME normalization range
  # bm.lb.res: biomass lb cutoff relative to biomass.max

  # get essential rxns: those whose KO decrease biomass by >90% (default)
  message("Identifying essential reactions...")
  bm0 <- get.opt.flux(model, bm.rgx, solv.pars=solv.pars)
  bms <- unlist(parallel::mclapply(1:length(model$rxns), function(i) get.opt.flux(model, bm.rgx, ko=i, solv.pars=solv.pars), mc.cores=nc))
  ess.idx <- which(bms<bm.lb.rel*bm0)
  # vmins of essential rxns to support at least biomass.max*0.1 (default)
  message("Computing the lower bound of normalization range...")
  m <- set.biomass.bounds(model, rgx=bm.rgx, lb=bm.lb.rel, relative=TRUE, solv.pars=solv.pars)
  suppressMessages(ess.vmins <- get.opt.fluxes(m, rxns=ess.idx, dir="min", nc=nc, solv.pars=solv.pars))
  # max of vmins
  max(ess.vmins)
}

get.norm.range.max <- function(model, rxns, range.min, nc=1L, bm.rgx="biomass", solv.pars=get.pars("lp", list())) {
  # get the upper bound value for the PRIME normalization range
  # rxns: the indices of a set of growth-associated reactions, i.e. output $i of get.prime.rxns()
  # range.min: the lower bound value for the PRIME normalization range, i.e. output of get.norm.range.min()

  # biomass.max under series of ubs for gr.rxns
  bnd <- max(model$ub)
  xs <- seq(range.min, bnd, by=0.1)
  bms <- unlist(parallel::mclapply(xs, function(i) {
  	m <- set.rxn.bounds(model, rxns, ubs=i)
  	get.opt.flux(m, bm.rgx, solv.pars=solv.pars)
  }, mc.cores=nc))
  # this is how range.max was computed in the code provided by Yizhak et al., doesn't make sense to me...
  dd <- abs(diff(abs(diff(bms))))
  xs[which.max(dd)+1]
}



