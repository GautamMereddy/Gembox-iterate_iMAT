###### reading/preparing metabolic models from files ######


### --- reading models from files --- ###

read.matlab.model <- function(fn, gene.rgx="[A-Za-z0-9_.-]+") {
  # read a model in the MATLAB (.mat) format
  # the resulting model is a list, and will contain a minimal set of necessary elements including:
  fields <- c("rxns", "mets", "S", "lb", "ub", "rowlb", "rowub", "genes", "rules", "c")
  # this function is a beta version, may not always work since different matlab models can have different formats; also, the `genes` field of the resulting model may need to be manually converted to gene symbols; will give an error if not all of the fields above are in the resulting model.

  if (!requireNamespace(c("R.matlab","rlist"), quietly=TRUE)) {
    stop("Packages \"R.matlab\" and \"rlist\" needed for this function to work.")
  }

  tmp <- R.matlab::readMat(fn)
  # tmp[[1]] contains the model but as an array, each element of the array being a field of the MATLAB struct... so need to use apply
  # for some weird reason each element of the array is a deeply nested list, so use rlist::list.flatten to make each of them a simple list of only one level
  # cannot use unlist here because it can flattern things directly to vectors and discarding NULLs/empty values.
  model <- apply(tmp[[1]], 1, function(x) rlist::list.flatten(x))
  # then manually unlist each of the simple lists, so that NULLs/empty values are not lost
  model <- lapply(model, function(x) {
    # if of length 1, just extract x[[1]], simplify to vector when appropriate
    if (length(x)==1) {
      res <- x[[1]]
      if (is.matrix(res) && (nrow(res)==1 || ncol(res)==1)) res <- as.vector(res)
    } else { # if length >1, carefully unlist into a vector so that NULLs/empty values are not lost
      res <- unname(sapply(x, function(y) {
        y <- as.vector(y) # y can be a (single-element or empty) matrix/array, so as.vector
        if (length(y)==0 || is.null(y) || y=="") y <- NA
        y
      }))
    }
    return(res)
  })

  # the genes are usually given as some IDs rather than gene symbols, I save these instead in the $gene.ids field
  model$gene.ids <- model$genes
  # need to manually add gene symbols to the $genes field, e.g. something like this (the example below is for converting human entrez IDs into gene symbols):
  #gid <- str_split(model$gene.ids, "_", simplify=TRUE)[,1]
  #library("org.Hs.eg.db")
  #mapp <- select(org.Hs.eg.db, keys=gid, columns=c("ENTREZID","SYMBOL"), keytype="ENTREZID") # many:1 mapping (with NA's)
  #all(gid==mapp$ENTREZID) # TRUE
  #model$genes <- ifelse(is.na(mapp$SYMBOL), model$gene.ids, mapp$SYMBOL)

  # rules mapping genes to reactions
  if ("rules" %in% names(model)) {
    model$rules <- stringr::str_replace_all(model$rules, "x\\([0-9]+\\)", function(x) paste0("x[",stringr::str_sub(x,3,-2),"]"))
    model$rules[is.na(model$rules)] <- ""
  }

  if ("grRules" %in% names(model) && !"rules" %in% names(model)) {
    #gene.rgx <- paste(stringr::str_replace_all(model$gene.ids, "(\\W)", "\\\\\\1"), collapse="|") # this has a potential issue, e.g. "A|A1" will match only A and not A1
    rules <- sapply(model$grRules, function(x) {
      if (is.na(x)) {
        x <- ""
      } else {
        x <- stringr::str_replace_all(x, " and ", " & ")
        x <- stringr::str_replace_all(x, " or ", " \\| ")
        x <- stringr::str_replace_all(x, gene.rgx, function(x) paste0("x[", match(x, model$gene.ids), "]"))
      }
      x
    })
    model$rules <- unname(rules)
  }

  if (any(grepl("x\\[NA\\]", model$rules))) warning("model$rules contain x[NA], please check manually!")

  if (!"rowlb" %in% names(model) && !"rowub" %in% names(model)) {
    if ("b" %in% names(model)) {
      model$rowlb <- model$b
      model$rowub <- model$b
    } else {
      model$rowlb <- rep(0, nrow(model$S))
      model$rowub <- rep(0, nrow(model$S))
    }
  }

  model$S <- Matrix(model$S, sparse=TRUE)

  # check for necessary fields
  tmp <- !fields %in% names(model)
  if (any(tmp)) warning("Failed to create these fields in the model, please fix manually: ", paste0(fields[tmp], collapse=", "), ".")
  message("Please double-check whether the essential fields in the model have the correct format.")

  model
}

read.sbml.model <- function(fn) {
  # read a model in the SBML (.xml) format
  # the resulting model is a list, and will contain a minimal set of necessary elements including:
  fields <- c("rxns", "mets", "S", "lb", "ub", "rowlb", "rowub", "genes", "rules", "c")
  # this function is a beta version, may not always work since different models can have different formats; also, the `genes` field of the resulting model may need to be manually converted to gene symbols; will give an error if not all of the fields above are in the resulting model.

  if (!requireNamespace("sybilSBML", quietly=TRUE)) {
    stop("Package \"sybilSBML\" needed for this function to work.")
  }

  #validateSBMLdocument(fn)
  tmp <- sybilSBML::readSBMLmod(fn)

  model <- list()
  model$description <- tmp@mod_desc
  model$modelName <- tmp@mod_name
  model$modelVersion <- tmp@version
  model$comps <- tmp@mod_compart
  model$mets <- tmp@met_id
  model$metNames <- tmp@met_name
  model$metCharges <- tmp@met_attr$charge
  model$metFormulas <- tmp@met_attr$chemicalFormula
  model$rxns <- tmp@react_id
  model$rxnNames <- tmp@react_name
  model$S <- tmp@S
  model$lb <- tmp@lowbnd
  model$ub <- tmp@uppbnd
  model$c <- tmp@obj_coef
  model$rowlb <- rep(0, length(tmp@met_id))
  model$rowub <- rep(0, length(tmp@met_id))
  model$b <- rep(0, length(tmp@met_id))
  model$rxnGeneMat <- tmp@rxnGeneMat
  model$subSystems <- tmp@subSys
  model$rules <- tmp@gprRules

  # the genes are usually given as some IDs rather than gene symbols, I save these instead in the $gene.ids field
  model$gene.ids <- tmp@allGenes
  # need to manually add gene symbols to the $genes field, e.g. something like this (the example below is for converting human HGNC IDs into gene symbols):
  #library(biomaRt)
  #mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  #tmp1 <- getBM(attributes=c("hgnc_id", "hgnc_symbol"), filters="hgnc_id", values=tmp@allGenes, mart=mart)
  #anyDuplicated(tmp1$hgnc_id) # 0
  #setdiff(tmp@allGenes, tmp1$hgnc_id) # some bad IDs, and others somehow not matched by biomart, manually map them to gene symbols based on the HGNC database
  #tmp1 <- rbind(tmp1, data.frame(hgnc_id=c("HGNC:HGNC:987","HGNC:HGNC:2898","HGNC:8780","HGNC:8789","HGNC:8790","HGNC:8903","HGNC:14423","HGNC:11240"), hgnc_symbol=c("BCKDHB","DLD","PDE4A","PDE6G","PDE6H","PGLS","RDH8","SPHK1")))
  #all(tmp@allGenes %in% tmp1$hgnc_id) # now TRUE
  #tmp1 <- as.data.table(tmp1)
  #model$genes <- tmp1[match(tmp@allGenes, hgnc_id), hgnc_symbol]
  model$genes <- model$gene.ids # for now as a placeholder

  # check for necessary fields
  tmp <- !fields %in% names(model)
  if (any(tmp)) warning("Failed to create these fields in the model, please fix manually: ", paste0(fields[tmp], collapse=", "), ".")
  message("Please double-check whether the essential fields in the model have the correct format.")

  model
}


### --- creating a conjunct (multicellular) model from a base model --- ###

make.conjunct.model <- function(model, n=2) {
  # join two copies of model into a conjunct multicellular model, where the cells share the extracellular space (metabolites and reactions in there)
  # n: the number of cells
  # model should contain at least the fields below in proper formats:
  fields <- c("rxns", "mets", "S", "lb", "ub", "rowlb", "rowub", "genes", "rules", "c")

  # check for necessary fields
  tmp <- !fields %in% names(model)
  if (any(tmp)) stop("These necessary items are missing in the model: ", paste0(fields[tmp], collapse=", "), ".")

  res <- model

  # mets in the extracellular space (emet.ids)
  tmp <- grepl("\\[e\\]$|_e$", model$mets)
  emet.ids <- which(tmp)
  # and the rest (i.e. intracellular)
  imet.ids <- which(!tmp)
  # cell compartment info of mets
  res$met.cell.ids <- ifelse(tmp, "e", "cell1") # base model as cell 1
  # save the original rowlb, rowub and b for emet.ids
  em.rowlb <- model$rowlb[emet.ids]
  em.rowub <- model$rowub[emet.ids]
  if ("b" %in% names(model)) em.b <- model$b[emet.ids]

  # exclusively extracellular reactions, these include the extracellular boundary reactions (those with names EX_.*, and equations like "xxx[e]-->") and those involving only extracellular metabolites
  tmp <- apply(model$S, 2, function(x) all(which(x!=0) %in% emet.ids))
  erxn.ids <- which(tmp)
  # and the rest
  not.erxn.ids <- which(!tmp)
  # cell compartment info of mets
  res$rxn.cell.ids <- ifelse(tmp, "e", "cell1") # base model as cell 1

  # add these info too
  res$imet.ids <- imet.ids
  res$emet.ids <- emet.ids
  res$irxn.ids <- not.erxn.ids
  res$erxn.ids <- erxn.ids

  # label the base model as cell1
  tmpf <- function(a, b="not.erxn.ids") {
    if (a %in% names(model)) {
      if (is.null(b)) {
        res[[a]] <<- paste0(model[[a]], "_cell1")
      } else {
        b <- get(b)
        res[[a]][b] <<- paste0(model[[a]][b], "_cell1")
      }
    }
  }
  for (i in c("rxns","rxnNames","subSystems")) tmpf(i)
  for (i in c("mets","metNames")) tmpf(i, "imet.ids")
  for (i in c("genes","gene.ids")) tmpf(i, NULL)

  # function for adding one cell at a time
  add.cell <- function(i) {

    # 1. simple concatenation or concatenation with added suffix indicating cell index
    tmpf <- function(a, b="not.erxn.ids", sfx=TRUE) {
      if (a %in% names(model)) {
        if (sfx) {
          sfx <- paste0("_cell",i)
          if (is.null(b)) {
            res[[a]] <<- c(res[[a]], paste0(model[[a]], sfx))
          } else {
            b <- get(b)
            res[[a]] <<- c(res[[a]], paste0(model[[a]][b], sfx))
          }
        } else {
          if (is.null(b)) {
            res[[a]] <<- c(res[[a]], model[[a]])
          } else {
            b <- get(b)
            res[[a]] <<- c(res[[a]], model[[a]][b])
          }
        }
      }
      NULL
    }
    for (x in c("lb","ub","c")) tmpf(x, sfx=FALSE)
    for (x in c("metFormulas","rowlb","rowub","b")) tmpf(x, "imet.ids", FALSE)
    for (x in c("rxns","rxnNames","subSystems")) tmpf(x)
    for (x in c("mets","metNames")) tmpf(x, "imet.ids")
    for (x in c("genes","gene.ids")) tmpf(x, NULL)
    # for rowlb, rowub and b, adjust for the fact that multiple cells are now producing/consuming metabolites
    res$rowlb[emet.ids] <<- res$rowlb[emet.ids] + em.rowlb
    res$rowub[emet.ids] <<- res$rowub[emet.ids] + em.rowub
    if ("b" %in% names(model)) res$b[emet.ids] <<- res$b[emet.ids] + em.b
    # cell compartment info
    res$rxn.cell.ids <<- c(res$rxn.cell.ids, rep(paste0("cell",i), length(not.erxn.ids)))
    res$met.cell.ids <<- c(res$met.cell.ids, rep(paste0("cell",i), length(imet.ids)))

    # 2. S matrix
    # the order of reactions is as above; for the (a) part, the matrix is just the original S unchanged but added more rows of zeros for the intracellular metabolites of the second cell; for the (b) part, the matrix is based on the not.erxn.ids columns of the original S matrix just with the rows corresponding to imet.ids "cut and pasted" to bottom
    # this is the (b) part:
    tmp <- res$S[, not.erxn.ids]
    tmp[imet.ids, ] <- 0
    tmp <- rbind(tmp, model$S[imet.ids, not.erxn.ids])
    # cbind (a) and (b)
    res$S <<- cbind(rbind(res$S, sparseMatrix(NULL, NULL, dims=c(length(imet.ids),ncol(res$S)))), tmp)

    # 3. rules (not using "grRules" or "rxnGeneMat", and these are also not included in the resulting model)
    # just need to correct the gene indices for the added cell
    tmp <- model$rules[not.erxn.ids]
    tmp <- stringr::str_replace_all(tmp, "[0-9]+", function(x) {
      x <- as.integer(x)
      x+length(res$genes)-length(model$genes) # need to -length(model$genes) because at this point the genes from the last cell have already been added to res$genes
    })
    res$rules <<- c(res$rules, tmp)
    # for erxn.ids that are mapped to genes, then need to recreate rules for these reactions as sth like "(rules.cell1) | (rules.cell2)"
    tmp <- rxns2genes(model, erxn.ids)
    dupg.erxn.ids <- erxn.ids[sapply(tmp, length)!=0]
    tmp <- stringr::str_replace_all(model$rules[dupg.erxn.ids], "[0-9]+", function(x) {
      x <- as.integer(x)
      x+length(res$genes)-length(model$genes) # need to -length(model$genes) because at this point the genes from the last cell have already been added to res$genes
    })
    res$rules[dupg.erxn.ids] <<- paste0("(",res$rules[dupg.erxn.ids],") | (",tmp,")")
  }

  # add the remaining (n-1) cells
  for (i in 2:n) add.cell(i)

  res
}

set.cell.fractions <- function(model.sc, model.mc, cell.fracs, nc=1L) {
  # adjust the fractions of different cells in a multicellular model
  # this is achieved by scaling the reaction bounds of the different cells, computed from the single-cell base model with FVA
  # nc: number of cores used by fva()
  # model.sc: the single-cell base model from which the multi-cellular model was constructed
  # model.mc: the multi-cellular model constructed using make.conjunct.model()
  # cell.fracs: a vector of cell fractions corresponding to the cells in model.mc (in the order of cell1, cell2, etc.)

  # bounds of reactions that are not exclusively extracellular
  bnds <- fva(model.sc, which(model.mc$rxn.cell.ids=="cell1"), nc=nc)
  res <- model.mc
  cells <- model.mc$rxn.cell.ids
  ncells <- uniqueN(cells[cells!="e"])
  if (length(cell.fracs)!=ncells) stop("length of cell.fracs not equal to the number of cells in the model.")
  cell.fracs <- cell.fracs/sum(cell.fracs)
  for (i in 1:length(cell.fracs)) {
    x <- cell.fracs[i]
    res <- set.rxn.bounds(res, paste0(bnds$rxn,"_cell",i), bnds$vmin*x, bnds$vmax*x)
  }
  res
}

join.models <- function(..., cell.fracs=1, nc=1L) {
  # join several models provided in ... into a conjunct multicellular model, where the cells share the extracellular space (metabolites and reactions in there)
  # all models in ... need to be based on the same base model, i.e. these models differ only by the reaction bound values
  # cell.fracs: a vector of cell fractions of the models, in the corresponding order, will rescale such that the total is 1; default 1 means each model will have fraction 1/n; reaction bounds will be adjusted based on the rescaled cell.fracs
  # nc: number of cores used by fva()
  # model should contain at least the fields below in proper formats:
  fields <- c("rxns", "mets", "S", "lb", "ub", "rowlb", "rowub", "genes", "rules", "c")

  models <- list(...)
  ncells <- length(models)
  for (m in models) {
    # check for necessary fields
    tmp <- !fields %in% names(m)
    if (any(tmp)) stop("These necessary items are missing in the model: ", paste0(fields[tmp], collapse=", "), ".")
  }

  # rescale model bounds based on cell.fracs
  if (length(cell.fracs)==1) cell.fracs <- rep(1, ncells) else if (length(cell.fracs)!=ncells) stop("length of cell.fracs not equal to the number of provided models.")
  cell.fracs <- cell.fracs/sum(cell.fracs)
  models <- lapply(1:ncells, function(i) {
    m <- models[[i]]
    bnds <- fva(m, nc=nc)
    m$lb <- bnds$vmin*cell.fracs[i]
    m$ub <- bnds$vmax*cell.fracs[i]
    m
  })
  
  res <- models[[1]]

  # mets in the extracellular space (emet.ids)
  tmp <- grepl("\\[e\\]$|_e$", res$mets)
  emet.ids <- which(tmp)
  # and the rest (i.e. intracellular)
  imet.ids <- which(!tmp)
  # cell compartment info of mets
  res$met.cell.ids <- ifelse(tmp, "e", "cell1") # first model as cell 1

  # exclusively extracellular reactions, these include the extracellular boundary reactions (those with names EX_.*, and equations like "xxx[e]-->") and those involving only extracellular metabolites
  tmp <- apply(res$S, 2, function(x) all(which(x!=0) %in% emet.ids))
  erxn.ids <- which(tmp)
  # and the rest
  not.erxn.ids <- which(!tmp)
  # cell compartment info of mets
  res$rxn.cell.ids <- ifelse(tmp, "e", "cell1") # first model as cell 1

  # add these info too
  res$imet.ids <- imet.ids
  res$emet.ids <- emet.ids
  res$irxn.ids <- not.erxn.ids
  res$erxn.ids <- erxn.ids

  # label the first model as cell1
  tmpf <- function(a, b="not.erxn.ids") {
    if (a %in% names(res)) {
      if (is.null(b)) {
        res[[a]] <<- paste0(res[[a]], "_cell1")
      } else {
        b <- get(b)
        res[[a]][b] <<- paste0(res[[a]][b], "_cell1")
      }
    }
  }
  for (i in c("rxns","rxnNames","subSystems")) tmpf(i)
  for (i in c("mets","metNames")) tmpf(i, "imet.ids")
  for (i in c("genes","gene.ids")) tmpf(i, NULL)

  # function for adding one cell at a time
  add.cell <- function(i) {
    model <- models[[i]]

    # 1. simple concatenation or concatenation with added suffix indicating cell index
    tmpf <- function(a, b="not.erxn.ids", sfx=TRUE) {
      if (a %in% names(model)) {
        if (sfx) {
          sfx <- paste0("_cell",i)
          if (is.null(b)) {
            res[[a]] <<- c(res[[a]], paste0(model[[a]], sfx))
          } else {
            b <- get(b)
            res[[a]] <<- c(res[[a]], paste0(model[[a]][b], sfx))
          }
        } else {
          if (is.null(b)) {
            res[[a]] <<- c(res[[a]], model[[a]])
          } else {
            b <- get(b)
            res[[a]] <<- c(res[[a]], model[[a]][b])
          }
        }
      }
      NULL
    }
    for (x in c("lb","ub","c")) tmpf(x, sfx=FALSE)
    for (x in c("metFormulas","rowlb","rowub","b")) tmpf(x, "imet.ids", FALSE)
    for (x in c("rxns","rxnNames","subSystems")) tmpf(x)
    for (x in c("mets","metNames")) tmpf(x, "imet.ids")
    for (x in c("genes","gene.ids")) tmpf(x, NULL)
    # for lb and ub of erxn.ids, and rowlb, rowub and b of emet.ids, treat as additive
    # the additive lb and ub is due to the consideration that all rxn bounds for each cell have been adjusted based on cell.fracs
    res$lb[erxn.ids] <<- res$lb[erxn.ids] + model$lb[erxn.ids]
    res$ub[erxn.ids] <<- res$ub[erxn.ids] + model$ub[erxn.ids]
    res$rowlb[emet.ids] <<- res$rowlb[emet.ids] + model$rowlb[emet.ids]
    res$rowub[emet.ids] <<- res$rowub[emet.ids] + model$rowub[emet.ids]
    if ("b" %in% names(model)) res$b[emet.ids] <<- res$b[emet.ids] + model$b[emet.ids]
    # cell compartment info
    res$rxn.cell.ids <<- c(res$rxn.cell.ids, rep(paste0("cell",i), length(not.erxn.ids)))
    res$met.cell.ids <<- c(res$met.cell.ids, rep(paste0("cell",i), length(imet.ids)))

    # 2. S matrix
    # the order of reactions is as above; for the (a) part, the matrix is just the original S unchanged but added more rows of zeros for the intracellular metabolites of the second cell; for the (b) part, the matrix is based on the not.erxn.ids columns of the original S matrix just with the rows corresponding to imet.ids "cut and pasted" to bottom
    # this is the (b) part:
    tmp <- res$S[, not.erxn.ids]
    tmp[imet.ids, ] <- 0
    tmp <- rbind(tmp, model$S[imet.ids, not.erxn.ids])
    # cbind (a) and (b)
    res$S <<- cbind(rbind(res$S, sparseMatrix(NULL, NULL, dims=c(length(imet.ids),ncol(res$S)))), tmp)

    # 3. rules (not using "grRules" or "rxnGeneMat", and these are also not included in the resulting model)
    # just need to correct the gene indices for the added cell
    tmp <- model$rules[not.erxn.ids]
    tmp <- stringr::str_replace_all(tmp, "[0-9]+", function(x) {
      x <- as.integer(x)
      x+length(res$genes)-length(model$genes) # need to -length(model$genes) because at this point the genes from the last cell have already been added to res$genes
    })
    res$rules <<- c(res$rules, tmp)
    # for erxn.ids that are mapped to genes, then need to recreate rules for these reactions as sth like "(rules.cell1) | (rules.cell2)"
    tmp <- rxns2genes(model, erxn.ids)
    dupg.erxn.ids <- erxn.ids[sapply(tmp, length)!=0]
    tmp <- stringr::str_replace_all(model$rules[dupg.erxn.ids], "[0-9]+", function(x) {
      x <- as.integer(x)
      x+length(res$genes)-length(model$genes) # need to -length(model$genes) because at this point the genes from the last cell have already been added to res$genes
    })
    res$rules[dupg.erxn.ids] <<- paste0("(",res$rules[dupg.erxn.ids],") | (",tmp,")")
  }

  # add the remaining (n-1) cells
  for (i in 2:ncells) add.cell(i)

  res
}
