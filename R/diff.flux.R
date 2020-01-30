###### functions for differential flux analysis ######


df.wilcox <- function(mat0, mat1, by=c(1,2), nc, padj.cutoff, r.cutoff, df.cutoff) {
  # a helper function to perform differential flux analysis with wilcoxon's rank-sum tests
  # if by=1, then mat0, mat1 are matrices with rows being reactions in matched orders, columns are samples, i.e. each row vector represents a sampling distribution for the flux of a reaction
  # if by=2, then rxns in columns
  # the test will compare mat1 relative to mat0
  # nc: number of cores to use

  by <- match.arg(as.character(by[1]), c("1","2"))
  dflux.test <- function(s0, s1) {
    # an inner helper function to run wilcoxon's rank-sum test
    tryCatch({
      wilcox.res <- wilcox.test(s0, s1)
      # p value
      wilcox.p <- wilcox.res$p.value
      # effect size for wilcoxon's rank-sum test: rank biserial correlation
      wilcox.r <- unname(1 - 2 * wilcox.res$statistic / (length(s0)*length(s1)))
      # another effect size measure: difference of mean fluxes (I use mean instead of median to facilitate some specific downstream analyses)
      m0 <- mean(s0)
      m1 <- mean(s1)
      data.table(lb0=min(s0), ub0=max(s0), mean0=m0, lb1=min(s1), ub1=max(s1), mean1=m1, diff.mean=m1-m0, r=wilcox.r, pval=wilcox.p)
    }, error=function(e) {
      data.table(lb0=NA, ub0=NA, mean0=NA, lb1=NA, ub1=NA, mean1=NA, diff.mean=NA, r=NA, pval=NA)
    })
  }
  if (by=="1") {
    res <- rbindlist(parallel::mclapply(1:nrow(mat0), function(i) dflux.test(mat0[i,], mat1[i,]), mc.cores=nc))
  } else if (by=="2") {
    res <- rbindlist(parallel::mclapply(1:ncol(mat0), function(i) dflux.test(mat0[,i], mat1[,i]), mc.cores=nc))
  }
  res[, padj:=p.adjust(pval, method="BH")]
  res[, dir:=ifelse(!(padj<padj.cutoff & abs(r)>r.cutoff & abs(diff.mean)>df.cutoff), 0, ifelse(diff.mean>0, 1, -1))]
}

df.fva <- function(model0, model1, rxns, coefs, nc, df.cutoff) {
  # a helper function to perform differential flux analysis with FVA, comparing model1 to model0
  # rxns and coefs: either rxns being a vector of reactions and coefs being 1 (df of each of these single reactions), or both being lists in the matched order in the case of df of combined fluxes

  ub0 <- mcmapply(get.opt.flux, rxns, coefs, MoreArgs=list(model=model0, dir="max"), mc.cores=nc)
  lb0 <- mcmapply(get.opt.flux, rxns, coefs, MoreArgs=list(model=model0, dir="min"), mc.cores=nc)
  ub1 <- mcmapply(get.opt.flux, rxns, coefs, MoreArgs=list(model=model1, dir="max"), mc.cores=nc)
  lb1 <- mcmapply(get.opt.flux, rxns, coefs, MoreArgs=list(model=model1, dir="min"), mc.cores=nc)
  m0 <- (ub0+lb0)/2
  m1 <- (ub1+lb1)/2
  res <- data.table(lb0=lb0, ub0=ub0, mean0=m0, lb1=lb1, ub1=ub1, mean1=m1, diff.mean=m1-m0)
  # add summary of flux differences: positive value means flux value changes towards the positive side, vice versa; 0 means unchanged
  `%gt%` <- function(a,b) a-b > df.cutoff
  `%eq%` <- function(a,b) abs(a-b) <= df.cutoff
  res[, dir:=ifelse(ub1 %gt% ub0 & lb1 %gt% lb0, 3,
             ifelse(ub0 %gt% ub1 & lb0 %gt% lb1, -3,
             ifelse(ub1 %gt% ub0 & lb1 %eq% lb0 | ub1 %eq% ub0 & lb1 %gt% lb0, 2,
             ifelse(ub0 %gt% ub1 & lb1 %eq% lb0 | ub1 %eq% ub0 & lb0 %gt% lb1, -2,
             ifelse(mean1 %gt% mean0, 1,
             ifelse(mean0 %gt% mean1, -1, 0))))))]
}

df.fva.x2 <- function(model, rxns0, rxns1, coefs, nc, df.cutoff) {
  # a helper function to perform differential flux analysis with FVA for bicellular model
  # rxns0 and rxns1 are in matched order, comparing rxns1 to rxns0
  # rxns(0/1) and coefs: either rxns(0/1) being a vector of reactions and coefs being 1 (df of each of these single reactions), or both being lists in the matched order in the case of df of combined fluxes

  ub0 <- mcmapply(get.opt.flux, rxns0, coefs, MoreArgs=list(model=model, dir="max"), mc.cores=nc)
  lb0 <- mcmapply(get.opt.flux, rxns0, coefs, MoreArgs=list(model=model, dir="min"), mc.cores=nc)
  ub1 <- mcmapply(get.opt.flux, rxns1, coefs, MoreArgs=list(model=model, dir="max"), mc.cores=nc)
  lb1 <- mcmapply(get.opt.flux, rxns1, coefs, MoreArgs=list(model=model, dir="min"), mc.cores=nc)
  m0 <- (ub0+lb0)/2
  m1 <- (ub1+lb1)/2
  res <- data.table(lb0=lb0, ub0=ub0, mean0=m0, lb1=lb1, ub1=ub1, mean1=m1, diff.mean=m1-m0)
  # add summary of flux differences: positive value means flux value changes towards the positive side, vice versa; 0 means unchanged
  `%gt%` <- function(a,b) a-b > df.cutoff
  `%eq%` <- function(a,b) abs(a-b) <= df.cutoff
  res[, dir:=ifelse(ub1 %gt% ub0 & lb1 %gt% lb0, 3,
             ifelse(ub0 %gt% ub1 & lb0 %gt% lb1, -3,
             ifelse(ub1 %gt% ub0 & lb1 %eq% lb0 | ub1 %eq% ub0 & lb1 %gt% lb0, 2,
             ifelse(ub0 %gt% ub1 & lb1 %eq% lb0 | ub1 %eq% ub0 & lb0 %gt% lb1, -2,
             ifelse(mean1 %gt% mean0, 1,
             ifelse(mean0 %gt% mean1, -1, 0))))))]
}

get.diff.flux <- function(model0, model1, rxns="all", method=c("wilcox","fva","both"), nsamples=4000, nc=1L, padj.cutoff=10/nsamples, r.cutoff=0.2, df.cutoff=1e-6) {
  # do differential flux analysis: model1 compared to model0, for the reactions specified in rxns (can be either indices or IDs as in model$rxns)
  # method: df method to use
  # nsamples: a single number, meaning to use the last N samples
  # padj.cutoff, r.cutoff, df.cutoff are used to determine the significantly changed reactions; the default values are arbitrary

  if (length(rxns)==1 && rxns=="all") rxns <- 1:length(model0$rxns) else rxns <- all2idx(model0, rxns)
  method <- match.arg(method)
  if (method %in% c("wilcox","both")) {
    if (!"sample" %in% names(model0) || !"sample" %in% names(model1)) stop("No sampling result found in models, cannot use method 'wilcox' or 'both'.")
    ns <- ncol(model0$sample$pnts)
    if (ns-nsamples<1e3) stop("At least ", nsamples+1e3, " samples needed.")
    mat0 <- model0$sample$pnts[rxns, (ns-nsamples+1):ns]
    ns <- ncol(model1$sample$pnts)
    if (ns-nsamples<1e3) stop("At least ", nsamples+1e3, " samples needed.")
    mat1 <- model1$sample$pnts[rxns, (ns-nsamples+1):ns]
    res <- df.wilcox(mat0, mat1, 1, nc, padj.cutoff, r.cutoff, df.cutoff)
  } else {
    res <- df.fva(model0, model1, rxns, 1, nc, df.cutoff)
  }
  if (method=="both") {
    tmp <- df.fva(model0, model1, rxns, 1, nc, df.cutoff)
    res[, c("lb0","ub0","lb1","ub1","dir.fva"):=tmp[, .(lb0, ub0, lb1, ub1, dir)]]
    setnames(res, "dir", "dir.wilcox")
  }
  res <- cbind(data.table(id=rxns, rxn=model0$rxns[rxns]), res)
}

get.diff.flux.x2 <- function(model, c0="cell1", c1="cell2", rxns="all", method=c("wilcox","fva","both"), nsamples=4000, nc=1L, padj.cutoff=10/nsamples, r.cutoff=0.2, df.cutoff=1e-6) {
  # diff.flux for bicellular model, between the two cells, i.e. c1 (cell2) vs c0 (cell1)
  # by default, rxns=="all": do diff.flux for all reactions that are not exlusively in the extracellular space; or specify rxns by their indices or IDs as in unicellular model
  # method: df method to use
  # nsamples: a single number, meaning to use the last N samples
  # padj.cutoff, r.cutoff, df.cutoff are used to determine the significantly changed reactions; the default values are arbitrary

  method <- match.arg(method)

  r0 <- sort(grep(c0, model$rxns))
  r1 <- sort(grep(c1, model$rxns)) # sort just to be safe that r0 and r1 are in matched order
  if (length(rxns)==1 && rxns=="all") {
    rxns <- stringr::str_sub(model$rxns[r0],1,-7)
  } else {
    if (is.character(rxns)) rxns <- all2idx(model, paste0(rxns,"_cell1"))
    tmp <- rxns %in% r0
    if (!all(tmp)) stop("Not all given rxns matched to non-extracellular reactions, please check!")
    tmp <- r0 %in% rxns
    r0 <- r0[tmp]
    r1 <- r1[tmp]
    rxns <- stringr::str_sub(model$rxns[r0],1,-7)
  }

  if (method %in% c("wilcox","both")) {
    if (!"sample" %in% names(model)) stop("No sampling result found in models, cannot use method 'wilcox' or 'both'.")
    ns <- ncol(model$sample$pnts)
    if (ns-nsamples<1e3) stop("At least ", nsamples+1e3, " samples needed.")
    mat0 <- model$sample$pnts[r0, (ns-nsamples+1):ns]
    mat1 <- model$sample$pnts[r1, (ns-nsamples+1):ns]
    res <- df.wilcox(mat0, mat1, 1, nc, padj.cutoff, r.cutoff, df.cutoff)
  } else {
    res <- df.fva.x2(model, r0, r1, 1, nc, df.cutoff)
  }
  if (method=="both") {
    tmp <- df.fva.x2(model, r0, r1, 1, nc, df.cutoff)
    res[, c("lb0","ub0","lb1","ub1","dir.fva"):=tmp[, .(lb0, ub0, lb1, ub1, dir)]]
    setnames(res, "dir", "dir.wilcox")
  }
  res <- cbind(data.table(id=r0, rxn=rxns), res)
}

get.diff.comb.flux <- function(model0, model1, rxns, coefs, method=c("wilcox","fva","both"), nsamples=4000, nc=1L, padj.cutoff=10/nsamples, r.cutoff=0.2, df.cutoff=1e-6) {
  # do differential flux analysis (model1 compared to model0) for the linear combination of fluxes of rxns
  # rxns is a list, each element is a vector of reaction indices or IDs (as in model$rxns)
  # coefs is a list in the matched order with rxns, each element contains the coefficient for the linear combination
  # this function will do df for each case corresponding to each element of the rxns and coefs lists
  # method: df method to use
  # nsamples: a single number, meaning to use the last N samples
  # padj.cutoff, r.cutoff, df.cutoff are used to determine the significantly changed reactions; the default values are arbitrary

  if (!is.list(rxns) || !is.list(coefs)) stop("rxns and coefs should both be lists.")
  rxns <- lapply(rxns, all2idx, model=model0)
  method <- match.arg(method)
  if (method %in% c("wilcox","both")) {
    if (!"sample" %in% names(model0) || !"sample" %in% names(model1)) stop("No sampling result found in models, cannot use method 'wilcox' or 'both'.")
    ns <- ncol(model0$sample$pnts)
    if (ns-nsamples<1e3) stop("At least ", nsamples+1e3, " samples needed.")
    mat0 <- mapply(function(x,c) colSums(model0$sample$pnts[x, (ns-nsamples+1):ns, drop=FALSE]*c), rxns, coefs)
    ns <- ncol(model1$sample$pnts)
    if (ns-nsamples<1e3) stop("At least ", nsamples+1e3, " samples needed.")
    mat1 <- mapply(function(x,c) colSums(model1$sample$pnts[x, (ns-nsamples+1):ns, drop=FALSE]*c), rxns, coefs)
    res <- df.wilcox(mat0, mat1, 2, nc, padj.cutoff, r.cutoff, df.cutoff)
  } else {
    res <- df.fva(model0, model1, rxns, coefs, nc, df.cutoff)
  }
  if (method=="both") {
    tmp <- df.fva(model0, model1, rxns, coefs, nc, df.cutoff)
    res[, c("lb0","ub0","lb1","ub1","dir.fva"):=tmp[, .(lb0, ub0, lb1, ub1, dir)]]
    setnames(res, "dir", "dir.wilcox")
  }
  if (is.null(names(rxns))) tmp <- 1:length(rxns) else tmp <- names(rxns)
  res <- cbind(data.table(id=tmp, res))
}

get.diff.flux.by.met <- function(model0, model1, mets="all", nsamples=4000, nc=1L, padj.cutoff=10/nsamples, r.cutoff=0.2, df.cutoff=1e-6) {
  # do differential flux analysis (model1 compared to model0, with wilcoxon tests) for the flux through each metabolite (i.e. either the production or consumption, they should be the same in magnitude, S*v=0 is assumed)
  # nsamples: a single number, meaning to use the last N samples
  # padj.cutoff, r.cutoff, df.cutoff are used to determine the significantly changed reactions; the default values are arbitrary
  # note: here, unlike df by rxn, fva cannot be used because in the cases involving reversible reactions, we don't know which reactions to use to represent the flux through a metabolite,
  # e.g. (1) x->y, (2) y->z, (3) y<=>w, here to represent the flux through y, in general I don't know whether I should use v1 (equivalent to v2+v3; this is when v3>0) or v2 (equivalent to v1-v3; this is when v3<0), and fva on v1 and v2 can yield different results

  if (length(mets)==1 && mets=="all") mets <- 1:length(model0$mets) else mets <- all2idx(model0, mets)
  if (!"sample" %in% names(model0) || !"sample" %in% names(model1)) stop("No sampling result found in models, cannot proceed.")
  ns <- ncol(model0$sample$pnts)
  if (ns-nsamples<1e3) stop("At least ", nsamples+1e3, " samples needed.")
  mat0 <- abs(model0$S[mets,,drop=FALSE]) %*% abs(model0$sample$pnts[, (ns-nsamples+1):ns]*c) / 2
  ns <- ncol(model1$sample$pnts)
  if (ns-nsamples<1e3) stop("At least ", nsamples+1e3, " samples needed.")
  mat1 <- abs(model1$S[mets,,drop=FALSE]) %*% abs(model1$sample$pnts[, (ns-nsamples+1):ns]*c) / 2
  res <- df.wilcox(mat0, mat1, 1, nc, padj.cutoff, r.cutoff, df.cutoff)
  res <- cbind(data.table(id=mets, met=model0$mets[mets], res))
}

get.diff.flux.by.met.x2 <- function(model, c0="cell1", c1="cell2", mets="all", nsamples=4000, nc=1L, padj.cutoff=10/nsamples, r.cutoff=0.2, df.cutoff=1e-6) {
  # diff.flux.by.met for bicellular model, between the two cells, i.e. c1 (_cell2) vs c0 (_cell1)
  # will just re-use get.diff.flux.by.met; but perform analysis always for all intracellular metabolites
  # in the result, metabolite id correspond to those in the single cellular model

  m0 <- sort(grep(c0, model$mets))
  m1 <- sort(grep(c1, model$mets))
  if (length(mets)==1 && mets=="all") {
    mets <- stringr::str_sub(model$mets[m0],1,-7)
  } else {
    if (is.character(mets)) mets <- all2idx(model, paste0(mets,"_cell1"))
    tmp <- mets %in% m0
    if (!all(tmp)) stop("Not all given mets matched to intracellular metabolites, please check!")
    tmp <- m0 %in% mets
    m0 <- m0[tmp]
    m1 <- m1[tmp]
    mets <- stringr::str_sub(model$mets[m0],1,-7)
  }

  if (!"sample" %in% names(model)) stop("No sampling result found in models, cannot proceed.")
  ns <- ncol(model$sample$pnts)
  if (ns-nsamples<1e3) stop("At least ", nsamples+1e3, " samples needed.")
  mat0 <- abs(model$S[m0,,drop=FALSE]) %*% abs(model$sample$pnts[m0, (ns-nsamples+1):ns]*c) / 2
  mat1 <- abs(model$S[m1,,drop=FALSE]) %*% abs(model$sample$pnts[m1, (ns-nsamples+1):ns]*c) / 2
  res <- df.wilcox(mat0, mat1, 1, nc, padj.cutoff, r.cutoff, df.cutoff)
  res <- cbind(data.table(id=m0, met=mets, res))
}

check.diff.flux.of.met <- function(model, dflux.res, mets) {
  # given dflux.res from e.g. get.diff.flux() or get.diff.flux.x2() and a set of metabolite in mets (indices of IDs), for each met get the diff.flux result of the reactions associated with it
  # return a list, each element contains results for each met in mets

  mets <- all2idx(model, mets)
  names(mets) <- model$mets[mets]
  lapply(mets, function(i) {
    tmp <- model$S[i,]
    rxn.ids <- which(tmp!=0)
    met.coefs <- tmp[rxn.ids]
    res <- data.table(met=model$mets[i], rxn=model$rxns[rxn.ids], equation=get.rxn.equations(model, rxn.ids), met.coef=met.coefs, diff.mean=dflux.res[match(rxn.ids, id), diff.mean])
    res[, diff.mean.met:=diff.mean*met.coef]
    res[, subsystem:=model$subSystems[rxn.ids]]
    res[order(-diff.mean.met)]
  })
}

pathway.gsea <- function(dflux.res, pathways, value.name="r", id.name="id") {
  # metabolic pathway enrichment with gsea, from the result of differential flux analysis with get.diff.flux or get.diff.flux.by.met
  # value.name: the variable name in dflux.res for the measure of flux difference
  # id.name: the variable name in dflux.res for the reaction/metabolite id

  if (!requireNamespace("fgsea", quietly=TRUE)) {
    stop("Package fgsea needed for this function to work.")
  }

  vec <- dflux.res[[value.name]]
  names(vec) <- dflux.res[[id.name]]
  res <- fgsea::fgsea(pathways, vec, nperm=1e4)
  res <- res[order(padj, pval)]
}

