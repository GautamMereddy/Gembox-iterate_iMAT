###### basic metabolic model utils ######


### --- ID conversions and other conversion functions --- ###

all2idx <- function(model, x) {
  # convert rxns or mets to indices; if x is already numeric, return x as is; if x is logical, return(which(x))
  # for now, only convert rxns or mets; rxnNames, metNames and genes have duplications

  if (is.null(x)) return(NULL)
  if (is.numeric(x)) return(x)
  if (is.logical(x)) return(which(x))

  tmpf <- function(name, x) {
    if (any(x %in% model[[name]])) {
      message("all2idx(): Input as ", name, ".")
      res <- match(x, model[[name]])
      tmp <- is.na(res)
      if (any(tmp)) warning("NA's returned for these IDs not found: ", paste(x[tmp], collapse=", "), ".", call.=FALSE, immediate.=TRUE)
      res
    } else NULL
  }

  for (i in c("rxns","mets")) {
    res <- tmpf(i, x)
    if (!is.null(res)) return(res)
  }
  stop("None of the provided IDs is found in model.")
}

get.biomass.idx <- function(model, rgx="biomass") {
  # get the index of the biomass reaction from the model
  # rgx: regex used to find the biomass reaction from model$rxns

  res <- grep(rgx, model$rxns, ignore.case=TRUE)
  if (length(res)==0) stop("No reaction matched by biomass regex.")
  if (length(res)>1) stop("Biomass regex matched to more than one reactions: ", paste(model$rxns[res], collapse=", "))
  res
}

get.transport.info <- function(model, mets="all", c1="c", c2="e") {
  # get the info on metabolites that are transported across membrane, between two compartments
  # mets: metabolites of interest, give as ID as in model$mets but w/o compartment suffix; default to "all" (all possible metabolites across the two compartments as specificed with c1 and c2)
  # c1 and c2: two compartments
  # return a list by metabolite, named by metabolite IDs as in model$mets but w/o compartment suffix; each element (for each metabolite) is a data.table, with columns: id (rxn indices); rxn (rxn IDs as in model$rxns); coef (coefficient for the met in c1 in the rxn); equ (reaction equation); gene (transporter genes mapped to rxn)

  mets.in <- mets
  mets1.rgx <- paste0("(.*)(\\[",c1,"\\]$|_",c1,"$)")
  mets2.rgx <- paste0("(.*)(\\[",c2,"\\]$|_",c2,"$)")
  mets1 <- stringr::str_match(model$mets, mets1.rgx)
  mets2 <- stringr::str_match(model$mets, mets2.rgx)
  mets <- intersect(mets1[,2], mets2[,2])
  if (length(mets.in)==1 && mets.in!="all" || length(mets.in)>1) mets <- intersect(mets, mets.in)
  mets <- mets[!is.na(mets)]
  if (length(mets)==0) {
    warning("No transportation reaction is found, NULL returned.")
    return(NULL)
  }
  mets1 <- all2idx(model, mets1[match(mets, mets1[,2]), 1])
  mets2 <- all2idx(model, mets2[match(mets, mets2[,2]), 1])
  tmpf <- function(m1, m2) {
    rxns <- intersect(mets2rxns(model, m1)[[1]], mets2rxns(model, m2)[[1]])
    if (length(rxns)==0) return(NULL)
    coefs <- model$S[m1, rxns]
    data.table(id=rxns, rxn=model$rxns[rxns], coef=coefs, equ=get.rxn.equations(model, rxns), gene=rxns2genes(model, rxns))
  }
  res <- mapply(tmpf, mets1, mets2, SIMPLIFY=FALSE)
  names(res) <- mets
  res <- res[!sapply(res, is.null)]
  if (length(res)==0) {
    warning("No transportation reaction is found, NULL returned.")
    return(NULL)
  }
  res
}

rxns2mets <- function(model, x, mode=c(0,-1,1), rev.mode=c(1,0), out="idx") {
  # map reactions to metabolites, return a list
  # mode: 0: all; -1: reactants; 1: products
  # rev.mode: for a reversible reaction, 1: mode based on the direction it is written in; 0: always return all
  # out: output mode, default to indices; or "mets", "metNames" etc.

  # match.arg, below are workarounds since it cannot match numerical arguments
  mode <- match.arg(as.character(mode[1]), c("0","-1","1"))
  rev.mode <- match.arg(as.character(rev.mode[1]), c("1","0"))

  idx <- all2idx(model, x)
  names(idx) <- x
  lapply(idx, function(i) {
    if (model$lb[i]<0 && rev.mode=="0") typ <- "0" else typ <- mode
    res <- switch(typ,
                  `0`=which(model$S[,i]!=0),
                  `-1`=which(model$S[,i]<0),
                  `1`=which(model$S[,i]>0))
    if (out!="idx") res <- model[[out]][res]
    res
  })
}

mets2rxns <- function(model, x, mode=c(0,-1,1), rev.mode=c(1,0), out="idx") {
  # map metabolites to reactions, return a list
  # mode: 0: all; -1: reactants; 1: products
  # rev.mode: for a reversible reaction, 1: mode based on the direction it is written in; 0: always return all
  # out: output mode, default to indices; or "rxns", "rxnNames" etc.

  mode <- match.arg(as.character(mode[1]), c("0","-1","1"))
  rev.mode <- match.arg(as.character(rev.mode[1]), c("1","0"))

  idx <- all2idx(model, x)
  names(idx) <- x
  lapply(idx, function(i) {
    if (rev.mode=="1") {
      res <- switch(mode,
                    `0`=which(model$S[i,]!=0),
                    `-1`=which(model$S[i,]<0),
                    `1`=which(model$S[i,]>0))
    } else if (rev.mode=="0") {
      res <- switch(mode,
                    `0`=which(model$S[i,]!=0),
                    `-1`=which(model$S[i,]<0 | (model$S[i,]>0 & model$lb<0)),
                    `1`=which(model$S[i,]>0 | (model$S[i,]<0 & model$lb<0)))
    }
    if (out!="idx") res <- model[[out]][res]
    res
  })
}

rxns2genes <- function(model, x) {
  # map reactions to gene names, return as a list

  idx <- all2idx(model, x)
  idx[idx==0] <- NA # NA will be returned where x is 0
  lapply(stringr::str_extract_all(model$rules[idx], "[1-9][0-9]*"), function(x) unique(model$genes[as.numeric(x)]))
}

genes2rxns <- function(model, genes, mode=c(0,1), out="idx") {
  # map gene symbols to reaction indices, return as a list
  # mode==0 for any reactions involving the gene; mode==1 for reactions where the gene is essential (corresponding to that if the gene is removed, then the reaction cannot happen based on model$rules)
  # out: output mode, default to indices; or "rxns", "rxnNames" etc.

  mode <- match.arg(as.character(mode[1]), c("0","1"))

  if (mode=="0") {
    r2g <- rxns2genes(model, 1:length(model$rules))
    res <- lapply(genes, function(gi) which(sapply(r2g, function(gns) gi %in% gns)))
  } else if (mode=="1") {
    gind <- match(genes, model$genes)
    res <- lapply(gind, function(gi) {
      tmp <- rep(0, length(model$genes))
      tmp[gi] <- -1
      which(exprs2rxns(model, tmp, 0)==-1)
    })
  }

  if (out!="idx") res <- lapply(res, function(x) model[[out]][x])
  res
}

exprs2fluxes <- function(model, x, na2zero=TRUE) {
  # map a numeric vector x of expression levels of genes to the flux levels of reactions in the model
  # x should contain either continuous expression values (then the output will also be continuous) or discrete values from {-1,0,1,NA} representing low/medium/high expression levels
  # x should be named by gene symbols as used in model$genes, or if it's unnamed and length being length(model$genes), assume it's already in the same order as model$genes
  # NA's in x will be kept and propagated during the conversion to flux values; if na2zero==TRUE, NA in the final flux vector will be replaced with 0

  if (is.null(names(x))) {
    if (length(x)==length(model$genes)) {
      message("exprs2fluxes(): Assuming the input vector is in the same order as model$genes.")
    } else stop("Input vector and model$genes have different lengths!")
  } else {
    x <- x[model$genes]
    if (all(is.na(x))) stop("Input doesn't contain any of the model genes!")
  }

  if (all(x %in% c(-1,0,1,NA))) {
    `&` <- function(a,b) {
      if (isTRUE((is.na(a) && b==-1) || (is.na(b) && a==-1))) return(-1) # if one is NA and the other is -1, for sure the result is -1; all other NA cases are undetermined and NA will be returned
      min(a,b)
    }
    `|` <- function(a,b) {
      if (isTRUE((is.na(a) && b==1) || (is.na(b) && a==1))) return(1) # if one is NA and the other is 1, for sure the result is 1; all other NA cases are undetermined and NA will be returned
      max(a,b)
    }
  } else {
    `&` <- function(a,b) min(a,b)
    `|` <- function(a,b) max(a,b)
  }

  res <- sapply(model$rules, function(i) eval(parse(text=i)))
  if (na2zero) res[is.na(res)] <- 0
  unname(res)
}

de2dflux <- function(model, x, na2zero=TRUE) {
  # map a numeric vector x of differential gene expression values (from {-1,0,1,NA}, representing decreased/unchanged/increased expressions) to the flux changes of reactions in the model (also {-1,0,1})
  # x should be named by gene symbols as used in model$genes, or if it's unnamed and length being length(model$genes), assume it's already in the same order as model$genes
  # NA's in x will be kept and propagated during the conversion to flux values; if na2zero==TRUE, NA in the final flux vector will be replaced with 0
  # other than the propagation of NA, this function reproduces the behavior of the MATLAB evalExpRule function from the original MTA code, although I really think the behavior for de2dflux should be no different from that of exprs2fluxes

  if (!all(x %in% c(-1,0,1,NA))) stop("Expression values in x should be from {-1,0,1,NA}.")
  if (is.null(names(x))) {
    if (length(x)==length(model$genes)) {
      message("de2dflux(): Assuming the input vector is in the same order as model$genes.")
    } else stop("Input vector and model$genes have different lengths!")
  } else {
    x <- x[model$genes]
    if (all(is.na(x))) stop("Input doesn't contain any of the model genes!")
  }

  `&` <- function(a,b) {
    if (isTRUE(a==0 || b==0 || a!=b)) return(0) # if one is NA and the other is 0, for sure the result is 0
    min(a,b) # all other NA cases are undetermined and NA will be returned
  }
  `|` <- function(a,b) {
    if (isTRUE(a==b)) return(a)
    return(a+b) # all NA cases are undetermined and NA will be returned
  }

  res <- sapply(model$rules, function(i) eval(parse(text=i)))
  if (na2zero) res[is.na(res)] <- 0
  unname(res)
}

get.rxn.equations <- function(model, x, dir=1, use.names=FALSE) {
  # get equations of reactions
  # dir: if not default, should be a vector with the same length as x, representing their fluxes, then for reversible reations with negative fluxes the direction of the equations will be reversed
  # use.names: if TRUE, use model$metNames with attempts to add compartment label; otherwise ues model$mets

  idx <- all2idx(model, x)
  if (use.names) {
    mets <- model$metNames
    suffix <- stringr::str_match(model$mets, "[\\[_](.)\\]?$")[,2]
    suffix <- paste0("[",stringr::str_sub(suffix,1,1),"]")
    mets <- paste0(mets, suffix)
  } else mets <- model$mets
  mapply(function(i, coef) {
    x <- model$S[,i]*sign(ifelse(coef==0, 1, coef))
    rs <- paste(trimws(paste(ifelse(x[x<0]==-1,"",-x[x<0]), mets[x<0])), collapse=" + ")
    ps <- paste(trimws(paste(ifelse(x[x>0]==1,"",x[x>0]), mets[x>0])), collapse=" + ")
    if (model$lb[i]>=0) arrow <- "-->" else arrow <- "<==>"
    paste(rs, arrow, ps)
  }, idx, dir, SIMPLIFY=TRUE)
}

subsystems2gsets <- function(model, by=c("rxn","met"), exclude.mets=NULL, exclude.mets.rgx="default", exclude.mets.degree=ncol(model$S), name="subSystems") {
  # create a list of rxns or mets (specified in "by") sets from the "subSystems" field of a metabolic model
  # if by metabolite, exclude.mets contains indices or IDs or mets to exclude, exclude.mets.rgx are the regex of high-degree mets to exclude, also can set degree cutoff with exclude.mets.degree (metabolites with degree higher than this will be excluded)

  by <- match.arg(by)
  if (is.null(model[[name]])) stop("subSystems not in model.\n")
  tmp <- data.table(path=model$subSystems, rxn.id=model$rxns)
  tmp <- tmp[!is.na(path) & path!="", .(rxn.id=list(rxn.id)), by=path]
  gsets <- tmp$rxn.id
  names(gsets) <- tmp$path
  if (by=="rxn") return(gsets)

  # if by metabolite, exclude metabolites
  mets.rm <- get.exclude.mets(model, exclude.mets, exclude.mets.rgx, exclude.mets.degree)
  mets.rm <- model$mets[mets.rm]
  gsets <- lapply(gsets, function(x) {
    mets <- rxns2mets(model, x, out="mets")
    mets <- unique(unlist(mets))
    setdiff(mets, mets.rm)
  })
}


### --- model manipulation, e.g. subsetting models, adding/removing reactions, etc. --- ###

subset.model <- function(model, i, j, rm.extra=FALSE, update.genes=FALSE, row.vars=NULL, col.vars=NULL) {
  # subset model like a matrix, i for metabolites, j for reactions; i and j can be logical or indices or IDs as in model$mets and model$rxns, but negative indices cannot be used
  # rm.extra: after subsetting, whether to further remove "empty" reactions and metabolites involved in no reaction
  # update.genes: whether to remove genes no longer in the model, this will be done by setting those genes to NA to avoid the need to completely rewrite model$rules; but note that this will affect genes2fluxes()
  # the fields in the fields variable below are assumed to be present in the model, and only these fields of the model are kept and updated, unless:
  # names of extra row (metabolite-related) or column (reaction-related) to update and keep in the returned model are given in row.vars and col.vars
  fields <- c("rxns","rxnNames","lb","ub","c","rules","subSystems","genes","mets","metNames","metFormulas","rowlb","rowub","b","S")
  fields <- c(fields, row.vars, col.vars)
  fields <- intersect(fields, names(model))

  i <- all2idx(model, i)
  j <- all2idx(model, j)
  if (rm.extra) {
    tmp <- model$S
    tmp[-i,] <- 0
    tmp[,-j] <- 0 # note tmp[-i,-j] <- 0 will not give what is wanted
    tmp <- tmp!=0
    i <- rowSums(tmp)>0
    j <- colSums(tmp)>0
  }
  names(fields) <- fields
  res <- lapply(fields, function(field) {
    x <- model[[field]]
    if (field %in% c("rxns","rxnNames","lb","ub","c","rules","subSystems", col.vars)) x <- x[j]
    if (field %in% c("mets","metNames","metFormulas","rowlb","rowub","b", row.vars)) x <- x[i]
    if (field=="S") x <- x[i,j]
    x
  })
  if (update.genes) {
    gns.keep <- as.numeric(unlist(stringr::str_extract_all(res$rules, "[0-9]+")))
    gns.rm <- setdiff(1:length(res$genes), gns.keep)
    res$genes[gns.rm] <- NA
  }
  res
}

c.model <- function(model1, model2, c.vars=NULL, cx.vars=NULL) {
  # used to simply concatenate two models, treat all rxns, mets and genes as separate
  # the fields in the fields variable below are assumed to be present in the model, and only these fields of the model are kept, unless:
  # names of extra fields to combine are given in c.vars (for those to be simply c()) and cx.vars (for those to be suffixed with "_1" and "_2" then c())
  fields <- c("rxns","rxnNames","lb","ub","c","rules","subSystems","genes","mets","metNames","metFormulas","rowlb","rowub","b","S")
  fields <- c(fields, c.vars, cx.vars)
  fields <- intersect(fields, intersect(names(model1), names(model2)))
  names(fields) <- fields
  lapply(fields, function(field) {
    x1 <- model1[[field]]
    x2 <- model2[[field]]
    if (field=="S") x <- rbind(cbind(x1, sparseMatrix(NULL, NULL, dims=c(nrow(x1), ncol(x2)))),
                               cbind(sparseMatrix(NULL, NULL, dims=c(nrow(x2), ncol(x1))), x2))
    if (field %in% c("lb","ub","c","rowlb","rowub","b","metFormulas", c.vars)) x <- c(x1, x2)
    if (field %in% c("rxns","rxnNames","subSystems","genes","mets","metNames", cx.vars)) x <- c(paste0(x1,"_1"), paste0(x2,"_2"))
    if (field=="rules") {
      x2[x2!="0"] <- stringr::str_replace_all(x2[x2!="0"], "[0-9]+", function(x) as.integer(x)+length(model1$genes))
      x <- c(x1, x2)
    }
    x
  })
}

add.constraint <- function(model, rxns, coefs, rowlb, rowub) {
  # add one new constraint to model (i.e. a row of model$S)
  # rxns: indices or IDs (as in model$rxns) of rxns in the contraint; coefs: the corresponding coefficients of rxns; rowlb, rowub: the bounds of the contraint

  rxns <- all2idx(model, rxns)
  tmp <- rep(0, ncol(model$S))
  tmp[rxns] <- coefs
  model$S <- rbind(model$S, tmp)
  model$rowlb <- c(model$rowlb, rowlb)
  model$rowub <- c(model$rowub, rowub)
  model
}

rm.rxns <- function(model, x, rm.extra=FALSE, update.genes=FALSE, row.vars=NULL, col.vars=NULL) {
  # remove reactions from model, given as logical or indices or IDs as in model$rxns
  # rm.extra: after subsetting, whether to further remove "empty" reactions and metabolites involved in no reaction
  # update.genes: whether to remove genes no longer in the model, this will be done by setting those genes to NA to avoid the need to completely rewrite model$rules; but note that this will affect genes2fluxes()
  # the default fields as in subset.model() are assumed to be present in the model, and only these fields of the model are kept and updated, unless:
  # names of extra row (metabolite-related) or column (reaction-related) to update and keep in the returned model are given in row.vars and col.vars

  x <- all2idx(model, x)
  j <- setdiff(1:ncol(model$S), x)
  subset.model(model, 1:nrow(model$S), j, rm.extra, update.genes, row.vars, col.vars)
}

add.rxn <- function(model, rxn, mets, coefs, lb, ub, rule="0", rxnName=rxn, subSystem=NA, metNames=mets, metFormulas=NA) {
  # add one new reaction to model
  # rxn is the ID of the new reaction to be added to model$rxns; rxnName will be added to model$rxnNames; and subSystem to model$subSystems
  # mets are the IDs of metabolites (as in model$mets) within the new reaction, new mets will be added automatically to model$mets; metNames are the names corresponding to mets, similar for metFormulas; coefs are their corresponding coefficients
  # for the new mets, their rowlb, rowub and b will be set to zeros
  # rule: a character vector of length one, specifying the mapping of genes to the rxn, in the format of model$rules but using gene symbols, e.g. "GENE1 & (GENE2 | GENE3)"; default to "0" meaning no gene mapping
  # These fields below are assumed to be present in the model, and only these fields of the model are kept and updated (i.e. the returned model will have other fields discarded):
  fields <- c("rxns","rxnNames","lb","ub","c","rules","subSystems","genes","mets","metNames","metFormulas","rowlb","rowub","b","S")

  model <- model[fields[fields %in% names(model)]]
  model$rxns <- c(model$rxns, rxn)
  model$rxnNames <- c(model$rxnNames, rxnName)
  model$subSystems <- c(model$subSystems, subSystem)
  model$lb <- c(model$lb, lb)
  model$ub <- c(model$ub, ub)
  model$c <- c(model$c, 0)
  tmp <- !mets %in% model$mets
  mets.add <- mets[tmp]
  nmets.add <- length(mets.add)
  model$mets <- c(model$mets, mets.add)
  model$metNames <- c(model$metNames, metNames[tmp])
  model$metFormulas <- c(model$metFormulas, metFormulas[tmp])
  model$S <- rbind(model$S, matrix(0, nrow=nmets.add, ncol=ncol(model$S)))
  model$rowlb <- c(model$rowlb, rep(0, nmets.add))
  model$rowub <- c(model$rowub, rep(0, nmets.add))
  model$b <- c(model$b, rep(0, nmets.add))
  x <- rep(0, nrow(model$S))
  x[match(mets, model$mets)] <- coefs
  model$S <- cbind(model$S, x)
  if (rules!="0") {
    gns <- stringr::str_extract_all(rules, "[^()&| ]+")[[1]]
    model$genes <- c(model$genes, gns[!gns %in% model$genes])
    rules <- stringr::str_replace_all(rules, "[^()&| ]+", function(x) paste0("x[",match(x,model$genes),"]"))
  }
  model$rules <- c(model$rules, rules)
  model
}

set.rxn.bounds <- function(model, rxns, lbs=NULL, ubs=NULL, relative=FALSE, nc=1L, solv.pars=get.pars("lp", list())) {
  # set the lb's and ub's of reactions, given as indices or IDs as in model$rxns
  # to set only lb's or only ub's, pass the other argument as NULL
  # relative: if TRUE, lb's and ub's will be set to v_max*lbs and v_max*ubs, respectively; in the case of reversible reactions (actual v_min<0), lb's will be set to v_min*lbs (v_min and v_max determined by FVA)
  # nc: number of cores passed to get.opt.fluxes(); solv.pars: solver parameters passed to get.opt.fluxes(), by default LP

  x <- all2idx(model, rxns)
  if (relative) {
    vmaxs <- get.opt.fluxes(model, x, "max", nc, solv.pars)
    if (!is.null(lbs)) {
      vmins <- get.opt.fluxes(model, x, "min", nc, solv.pars)
      model$lb[x] <- lbs * ifelse(vmins<0, vmins, vmaxs)
    }
    if (!is.null(ubs)) model$ub[x] <- ubs * vmaxs
  } else {
    if (!is.null(lbs)) model$lb[x] <- lbs
    if (!is.null(ubs)) model$ub[x] <- ubs
  }
  model
}

set.medium <- function(model, medium, cells.per.ml=2e5, cell.gdw=4e-10, dbl.hr=24) {
  # set model constraint based on medium composition
  # medium: a data.table with the first column being names of metabolites as in model$mets (but w/o compartment suffix), second column being concentration in mM in the medium
  # cells.per.ml: cell count per mL; cell.gdw: dry weight per cell in grams; dbl.hr: cell doubling time in hours (default to approximate values for HeLa cells in a reasonably normal culture)
  # for metabolites transported by a single reaction, the corresponding lb or ub will be adjusted; for metabolites transported by multiple reactions, will add rows to model$S constraining the summed transport fluxes
  # metabolites whose info is not given in medium is not touched (regarded as unknown rather than definitely not present in the medium)

  m <- copy(as.data.table(medium))
  setnames(m, c("met", "conc"))
  tx <- get.transport.info(model, m$met, "e", "c")
  conc <- m[match(names(tx), met), -conc] / (1000*cells.per.ml*cell.gdw*dbl.hr)

  for (i in 1:length(tx)) {
    x <- tx[[i]]
    if (nrow(x)==1) {
      bnd <- conc[i] / x$coef
      if (x$coef>0 && bnd>model$lb[x$id]) model <- set.rxn.bounds(model, x$id, lbs=bnd)
      if (x$coef<0 && bnd<model$ub[x$id]) model <- set.rxn.bounds(model, x$id, ubs=bnd)
    } else {
      model <- add.constraint(model, x$id, x$coef, conc[i], 1e10)
    }
  }
  model 
}

rm.blocked.rxns <- function(model, nc=1L, rm.extra=FALSE, update.genes=FALSE, solv.pars=get.pars("lp", list())) {
  # remove blocked/deadend reactions, i.e. reactions with fixed flux of 0
  # These fields below are assumed to be present in the model, and only these fields of the model are kept and updated using subset.model() (i.e. the returned model will have other fields discarded):
  # "rxns","rxnNames","lb","ub","c","rules","subSystems","genes","mets","metNames","metFormulas","rowlb","rowub","b","S"
  # solv.pars: solver parameters passed to fva(), by default LP

  fva.res <- fva(model, nc=nc, solv.pars=solv.pars)
  x <- fva.res$vmin==0 & fva.res$vmax==0
  rm.rxns(model, x, rm.extra, update.genes)
}

convert.rev.rxns <- function(model) {
  # convert each reversible rxn in model into two irreversible rxns representing the positive and negative halves, both with lb==0
  # the positive half will keep its original ID (rxns) and name (rxnNames), the ID and name of the negative half will be the original name suffixed with "_REV"
  # These fields below are assumed to be present in the model
  # "rxns","rxnNames","lb","ub","c","rules","subSystems","S"

  rev.idx <- which(model$lb<0)
  nrev <- length(rev.idx)
  model$S <- cbind(model$S, -model$S[, rev.idx])
  rev.lb <- model$lb[rev.idx]
  model$lb[rev.idx] <- 0
  model$lb <- c(model$lb, rep(0, nrev))
  model$ub <- c(model$ub, -rev.lb)
  model$rxns <- c(model$rxns, paste0(model$rxns[rev.idx],"_REV"))
  model$rxnNames <- c(model$rxnNames, paste0(model$rxnNames[rev.idx],"_REV"))
  model$c <- c(model$c,  model$c[rev.idx])
  model$rules <- c(model$rules, model$rules[rev.idx])
  model$subSystems <- c(model$subSystems, model$subSystems[rev.idx])
  model
}


### --- functions involving metabolic network as a graph --- ###

get.exclude.mets <- function(model, mets=NULL, rgx="default", degree=ncol(model$S)) {
  # a helper function used by others for getting a set of metabolites that will be excluded from metabolic network; return a vector of metabolite indices
  # mets: mets to be excluded (either indices or IDs or names); default to nothing
  # rgx: regex of mets to be excluded; default to some high degree metabolites (see below, works for mets format like "h_c" or "h[c]")
  # degree: exclude metabolites with degree (number of edges to reactions) greater than this; default effect is that no further mets are excluded

  exclude.mets <- all2idx(model, mets)
  exclude.mets <- c(exclude.mets, which(Matrix::rowSums(model$S!=0)>degree))
  if (!(is.null(rgx) || is.na(rgx) || rgx=="")) {
    if (rgx=="default") rgx <- "^h[\\[_].\\]?$|^oh1[\\[_].\\]?$|^h2o[\\[_].\\]?$|^atp[\\[_].\\]?$|^adp[\\[_].\\]?$|^pi[\\[_].\\]?$|^ppi[\\[_].\\]?$|^coa[\\[_].\\]?$|^o2[\\[_].\\]?$|^co2[\\[_].\\]?$|^nadp[\\[_].\\]?$|^nadph[\\[_].\\]?$|^nad[\\[_].\\]?$|^nadh[\\[_].\\]?$|^fad[\\[_].\\]?$|^fadh2[\\[_].\\]?$|^na1[\\[_].\\]?$|^so4[\\[_].\\]?$|^nh4[\\[_].\\]?$|^cl[\\[_].\\]?$"
    emd <- grep(rgx, model$mets)
    exclude.mets <- unique(c(exclude.mets, emd))
  }
  # print a message about the removed metabolites
  tmp <- model$mets[exclude.mets]
  tmp <- unique(stringr::str_replace(tmp, "[\\[_].\\]?$", ""))
  message("get.exclude.mets(): The following metabolites are excluded:")
  message(paste(tmp, collapse=", "))
  exclude.mets
}

s2igraph <- function(model, exclude.mets=NULL, exclude.mets.rgx="default", exclude.mets.degree=ncol(model$S)) {
  # create an igraph bipartite graph from the model S matrix (directed graph, unweighted)
  # exclude.mets are mets to be excluded (either indices or IDs or names)
  # exclude.mets.rgx: regex of some high degree metabolites to be excluded by default; if "default", will use the default as in get.exclude.mets()
  # exclude.mets.degree: exclude metabolites with degree greater than this

  s <- as.matrix(model$S) # convert to "dense" matrix, since the igraph::graph_from_incidence_matrix function contains bugs working with sparse matrix
  # exclude metabolites
  exclude.mets <- get.exclude.mets(model, exclude.mets, exclude.mets.rgx, exclude.mets.degree)
  s[exclude.mets, ] <- 0
  rownames(s) <- model$mets
  colnames(s) <- model$rxns

  # create directed bipartite graph between mets and rxns
  # reversible reactions: bi-directional edges
  tmp <- s
  tmp[, model$lb>=0] <- 0
  tmp <- abs(tmp)
  g0 <- igraph::graph_from_incidence_matrix(tmp, directed=TRUE, mode="all")
  # non-reversible reactions: one-way edges
  ## edges from rxns to the product mets
  tmp <- s
  tmp[, model$lb<0] <- 0
  tmp[tmp<0] <- 0
  g1 <- igraph::graph_from_incidence_matrix(tmp, directed=TRUE, mode="in")
  ## edges from the reactant mets to rxns
  tmp <- s
  tmp[, model$lb<0] <- 0
  tmp[tmp>0] <- 0
  tmp <- -tmp
  g2 <- igraph::graph_from_incidence_matrix(tmp, directed=TRUE, mode="out")
  # combined graph (will be bipartite)
  g <- igraph::union(g0, g1, g2)
  # fix the node type attributes
  `%nna%` <- function(a, b) ifelse(is.na(a), b, a)
  igraph::V(g)$type <- igraph::V(g)$type_1 %nna% igraph::V(g)$type_2 %nna% igraph::V(g)$type_3
  g <- igraph::delete_vertex_attr(g, "type_1")
  g <- igraph::delete_vertex_attr(g, "type_2")
  g <- igraph::delete_vertex_attr(g, "type_3")
}

get.neighborhood <- function(model, ids, order=1, exclude.mets=NULL, exclude.mets.rgx="default", exclude.mets.degree=ncol(model$S)) {
  # given the IDs (as in model$rxns or model$mets) of either a set of rxns or a set of mets, return the IDs of the rxns or mets with distance<=order from each of the given ones (as a list). by default order=1 means the rxns sharing a met or the mets within the same rxn. Whether rxns or mets are provided will be decided automatically.
  # exclude.mets are mets to be excluded (either indices or IDs or names)
  # exclude.mets.rgx: regex of some high degree metabolites to be excluded by default; if "default", will use the default as in get.exclude.mets()
  # exclude.mets.degree: exclude metabolites with degree greater than this

  gb <- s2igraph(model, exclude.mets, exclude.mets.rgx, exclude.mets.degree)

  if (any(ids %in% model$rxns)) {
    # project bipartite graph into the graph of rxns
    gp <- igraph::bipartite_projection(gb, which="true")
  } else if (any(ids %in% model$mets)) {
    # project bipartite graph into the graph of mets
    gp <- igraph::bipartite_projection(gb, which="false")
  }
  # get the neighborhood
  res <- lapply(igraph::ego(gp, order=order, nodes=ids), names)
  names(res) <- ids
  res
}

get.path <- function(model, from, to, shortest=TRUE, exclude.mets=NULL, exclude.mets.rgx="default", exclude.mets.degree=ncol(model$S)) {
  # if shortest=TRUE, print the shortest (directed) path(s) between a pair of nodes (mets or rxns, can be mixed) in the metabolic network, given as IDs as in model$mets or model$rxns. The path(s) will contain both reactions and metabolites along the way. Whether rxns or mets are provided will be decided automatically.
  # will also return a list of shortest paths, each element per path being also a list with $path containing a vector of the entire path, $mets containing the path containing only mets, and $rxns containing the path containing only rxns
  # if shortest=FALSE, will get all simple paths, with the current implementation with igraph::all_simple_paths() this is *extremely* slow that it's virually impractical.
  # exclude.mets are mets to be excluded (either indices or IDs or names)
  # exclude.mets.rgx: regex of some high degree metabolites to be excluded by default; if "default", will use the default as in get.exclude.mets()
  # exclude.mets.degree: exclude metabolites with degree greater than this

  gb <- s2igraph(model, exclude.mets, exclude.mets.rgx, exclude.mets.degree)

  if (shortest) {
    tmp <- igraph::all_shortest_paths(gb, from, to, mode="out", weights=NA)
    tmp <- lapply(tmp$res, function(x) {
      x <- names(x)
      xx <- x
      # format and print results
      if (x[1] %in% model$mets) {
        mets <- x[seq(1, length(x), 2)]
        rxns <- x[seq(2, length(x), 2)]
        xx[seq(2, length(xx), 2)] <- paste0("--(", xx[seq(2, length(xx), 2)], ")->")
      } else if (x[1] %in% model$rxns) {
        rxns <- x[seq(1, length(x), 2)]
        mets <- x[seq(2, length(x), 2)]
        xx[seq(1, length(xx), 2)] <- paste0("--(", xx[seq(1, length(xx), 2)], ")->")
      }
      xx <- paste(xx, collapse=" ")
      list(prt=xx, path=x, mets=mets, rxns=rxns)
    })
    # print result
    print(sapply(tmp, function(x) x$prt))
    res <- lapply(tmp, function(x) x[c("path", "mets", "rxns")])
  } else {
    tmp <- igraph::all_simple_paths(gb, from, to, mode="out")
    # extract result: not implemented yet
    res <- tmp # a place-holder
  }
  invisible(res)
}

