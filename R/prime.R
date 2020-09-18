###### the PRIME algorithm (Yizhak et al. eLife 2014) ######


prime <- function(model, expr, gr, padj.cutoff=0.05, nc=1L, bm.rgx="biomass", default.bnd=1e3, bm.epsil=1e-4, min.step.size=0.1, bm.lb.rel=0.1, solv.pars=list()) {
  # function to run PRIME
  # expr: a gene-by-sample matrix, need to has rownames of gene symbols as those in model$genes (but doesn't need to be of same length or in order)
  # gr: cell growth rate measures across samples, corresponding to the order or columns of expr
  # padj.cutoff: used to define significant growth-associated rxns
  # default.bnd, bm.epsil, min.step.size passed to shrink.model.bounds(); bm.lb.res passed to get.norm.range.min()
  # return a list of updated models, one for each sample in the corresponding order

  solv.pars <- get.pars("lp", solv.pars)
  # split each reversible rxn into two forward/backward rxns
  model <- convert.rev.rxns(model)
  # shrink rxn bounds
  model <- shrink.model.bounds(model, rxns="default", default.bnd=default.bnd, bm.epsil=bm.epsil, relative=FALSE, min.step.size=min.step.size, bm.rgx=bm.rgx, solv.pars=solv.pars)
  # find growth-associated rxns and compute their rxn-level values for each sample
  prm.rxns <- get.prime.rxns(model, expr=expr, gr=gr, nc=nc, padj.cutoff=padj.cutoff)
  # compute the lb and ub of the normalization range
  rmin <- get.norm.range.min(model, bm.lb.rel=bm.lb.rel, nc=nc, bm.rgx=bm.rgx, solv.pars=solv.pars)
  rmax <- get.norm.range.max(model, rxns=prm.rxns$i, range.min=rmin, nc=nc, bm.rgx=bm.rgx, solv.pars=solv.pars)
  # compute and set final rxn ub values for each sample
  ubs <- prm.rxns$x * (rmax-rmin) + rmin
  models <- apply(ubs, 1, function(x) {
  	m <- model
  	m$ub[prm.rxns$i] <- x
  	m
  })
  # recombine splitted reversible rxns
  lapply(models, revert.rev.rxns)
}

get.prime.rxns <- function(model, expr, gr, nc=1L, padj.cutoff=0.05) {
  # select growth-associated rxns for PRIME (i.e. rxns with significant correlation with the provided growth rate data across multiple samples)
  # expr: a gene-by-sample matrix, need to has rownames of gene symbols as those in model$genes (but doesn't need to be of same length or in order)
  # gr: cell growth rate measures across samples, corresponding to the order or columns of expr
  # return list(i, x), i is the indices of the resulting growth-associated rxns, x is the rxn-by-sample matrix containing values that correspond to the (E-min(E))/(max(E)-min(E)) part in the Eq. (4) from Yizhak et al.

  expr <- expr[match(model$genes, rownames(expr)), ]
  if (all(is.na(expr))) stop ("expr doesn't contain any of the model genes!")

  # map gene expr values to rxns values by taking the mean
  mat <- sapply(model$rules, function(x) {
  	i <- as.integer(stringr::str_extract_all(x, "[0-9]+")[[1]])
  	colMeans(expr[i,,drop=FALSE], na.rm=TRUE)
  })
  mat[is.nan(mat)] <- NA

  # correlation between rxn values and growth rates across samples for each rxn
  cor.res <- rbindlist(mclapply(1:ncol(mat), function(i){
  	tryCatch({
  	  a <- cor.test(mat[,i], gr, method="spearman")
  	  data.table(rho=a$estimate, pval=a$p.value)
  	}, error=function(e) data.table(rho=NA, pval=NA))
  }, mc.cores=nc), idcol="id")
  cor.res[, padj:=p.adjust(pval,"BH")]
  # select rxns with signifcant correlations and compute the normalized rxn values
  cor.res <- cor.res[padj<padj.cutoff]
  if (nrow(cor.res)==0) stop("No significant growth-associated reactions with padj<", padj.cutoff, ".") else message("Found ", nrow(cor.res), " growth-associated reactions with padj<", padj.cutoff, ".")
  mat <- mat[, cor.res$id] * sign(cor.res$rho)
  mat <- do.call(rbind, mclpply(1:ncol(mat), function(i) {
  	x <- mat[,i]
  	rng <- range(x, na.rm=TRUE)
  	(x-rng[1]) / diff(rng)
  }, mc.cores=nc))
  if (!is.null(colnames(expr))) colnames(mat) <- colnames(expr)

  list(i=cor.res$id, x=mat)
}

get.norm.range.min <- function(model, bm.lb.rel=0.1, nc=1L, bm.rgx="biomass", solv.pars=get.pars("lp", list())) {
  # get the lower bound value for the PRIME normalization range
  # bm.lb.res: biomass lb cutoff relative to biomass.max

  # get essential rxns: those whose KO decrease biomass by >90% (default)
  bm0 <- get.opt.flux(model, bm.rgx, solv.pars=solv.pars)
  bms <- unlist(mclapply(1:length(model$rxns), function(i) get.opt.flux(model, bm.rgx, ko=i, solv.pars=solv.pars), mc.cores=nc))
  ess.idx <- which(bms<bm.lb.rel*bm0)
  # vmins of essential rxns to support at least biomass.max*0.1 (default)
  m <- set.biomass.bounds(model, rgx=bm.rgx, lb=bm.lb.rel, relative=TRUE, solv.pars=solv.pars)
  ess.vmins <- get.opt.fluxes(m, rxns=ess.idx, dir="min", nc=nc, solv.pars=solv.pars)
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
  bms <- unlist(mclapply(xs, function(i) {
  	m <- set.rxn.bounds(model, rxns, ubs=i)
  	get.opt.flux(m, bm.rgx, solv.pars=solv.pars)
  }, mc.cores=nc))
  # this is how range.max was computed in the code provided by Yizhak et al., doesn't make sense to me...
  dd <- abs(diff(abs(diff(bms))))
  xs[which.max(dd)+1]
}



