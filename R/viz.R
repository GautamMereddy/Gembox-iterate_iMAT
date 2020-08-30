###### functions for visualization of metabolic models ######


map.colors <- function(x, cols=c("blue2","grey70","red2"), trim=FALSE, lims=NULL, mid=NULL) {
  # map a numeric vector x to colors (i.e. color-code for x)
  # cols: a vector of colors (i.e. a spectrum of arbitrary length), corresponding to low to high x values
  # trim: whether to "remove" more extreme values in x; if not FALSE, can provide a numerical value, meaning that values outside of the range mean(x) +/- trim*sd(x) will be set to the corresponding boundary values; if TRUE corresponding to trim=3
  # lims: vector of two, low and high values correspond to the two extreme colors, default to range(x), if smaller than range(x), x values outside of the range will be mapped to the two extreme colors
  # mid: value corresponding to the mid-point color in cols (if cols have an even number colors the mid-point color will be extrapolated); default (NULL) means not forcing the mid-point color to correspond to any particular value
  # return a vector of colors corresponding to x

  if (isTRUE(trim)) trim <- 3
  if (is.numeric(trim)) {
    meanx <- mean(x)
    sdx <- sd(x)
    lb <- meanx - trim*sdx
    ub <- meanx + trim*sdx
    x[x<lb] <- lb
    x[x>ub] <- ub
  }
  if (is.null(lims)) {
    lims <- range(x)
  } else {
    if (lims[1]>min(x)) x[x<lims[1]] <- lims[1]
    if (lims[2]<max(x)) x[x>lims[2]] <- lims[2]
  }

  mapc <- function(x, cols, lims) {
    # helper function, map x to colors
    n <- uniqueN(x)
    if (n==1) {
      if (lims[1]==lims[2]) {
        cols <- colorRampPalette(cols)(2*length(cols)+1)
        res <- rep(cols[(length(cols)+1)/2], length(x))
      } else {
        cols <- colorRampPalette(cols)(20)
        bins <- as.numeric(cut(c(x[1],lims), 20))[1]
        res <- rep(cols[bins], length(x))
      }
    } else {
      nco <- round(n*diff(lims)/diff(range(x)))
      if (nco<50) nco <- 50
      if (nco>500) nco <- 500
      cols <- colorRampPalette(cols)(nco)
      bins <- as.numeric(cut(c(lims, x), nco))[-1:-2]
      res <- cols[bins]
    }
    res
  }

  if (!is.null(mid)) {
    if (length(cols) %% 2 == 0) {
      tmp <- colorRampPalette(cols[c(length(cols)/2, length(cols)/2+1)])(3)[2]
      cols <- c(cols[1:(length(cols)/2)], tmp, cols[(length(cols)/2+1):length(cols)])
    }
    if (mid<=lims[1]) {
      res <- mapc(x, cols[((length(cols)+1)/2):length(cols)], c(mid, lims[2]))
    } else if (mid>=lims[2]) {
      res <- mapc(x, cols[1:((length(cols)+1)/2)], c(lims[1], mid))
    } else {
      if (mid<mean(lims)) lims[1] <- mid-(lims[2]-mid) else lims[2] <- mid+mid-lims[1]
      res <- rep("x", length(x))
      # lower half
      l <- x<=mid
      res[l] <- mapc(x[l], cols[1:((length(cols)+1)/2)], c(lims[1], mid))
      # upper half
      u <- x>=mid
      res[u] <- mapc(x[u], cols[((length(cols)+1)/2):length(cols)], c(mid, lims[2]))
    }
  } else res <- mapc(x, cols, lims)

  res
}


plot.model <- function(model, rxns=NULL, fluxes=NULL, dfluxes=NULL, mets=NULL, exclude.mets.rgx="default", dup.mets.rgx="default", use.aes=c("both","color","width"), abs.dflux=FALSE, cols=c("blue2","grey70","red2"), lwds=c(5,20), label.value=c("none","flux","dflux"), layout=c("","layout_with_fr","layout_nicely","layout_randomly","layout_as_star","layout_as_tree","layout_as_bipartite","layout_in_circle","layout_on_sphere","layout_on_grid","layout_with_dh","layout_with_gem","layout_with_graphopt","layout_with_kk","layout_with_lgl","layout_with_mds","layout_with_sugiyama"), seed=NULL, width=NULL, height=NULL) {
  # generate an interactive network plot for a metabolic model, can also incorporate fluxes and dfluxes data
  # model: the base metabolic model
  # rxns: reactions to plot; mets: metabolites to include in the plot; if provide rxns but not mets, default to mets in the rxns, vice versa
  # fluxes: the flux values corresponding to the reactions in rxns; if plotting this, the arrows of the reversible reactions will reflect the direction of fluxes, otherwise draw double arrows for reversible reactions
  # dfluxes: the values of flux changes corresponding to the reactions in rxns
  # exclude.mets.rgx: regex for names of metabolites (as in model$mets) to be excluded from the plots; the default regex works for some of the high-degree mets in recon1 and iMM1415
  # dup.mets.rgx: after keeping the mets in mets.ids and excluding those in exclude.mets.rgx, for the remaining mets, use dup.mets to specify regex of mets to be plot as separate nodes for each reaction, when they are recurrent in multiple reactions; the default is the same as exclude.mets, so these will be excluded; to duplicate these instead of removing them, set exclude.mets to NULL
  # use.aes: use line color, or line width, or both; for dfluxes, width can only represent magnitude; if both fluxes and dfluxes are given, then will always use line width for flux and color for dflux
  # abs.dflux: whether to treat dfluxes as change in absolute fluxes (i.e. magnitude of fluxes) or "raw" changes (i.e. dependent on the direction of reactions); if FALSE and also plotting fluxes, for reversible reactions the color will be adjusted accordingly, if FALSE and not plotting fluxes, the dfluxes of reversible reactions will always be plotted using the "positive" colors and the arrows of the reactions will correspond to the direction of change
  # cols: a vector of colors for negative -> positive values if plotting dfluxes; if plotting fluxes, the mid-point color and right-half of the colors will be used for low -> high fluxes
  # lwds: a vector of length 2, range of line widths
  # label.value: whether to print the flux or dflux values in the visualization
  # layout: graph layout for visNetwork::visIgraphLayout
  # seed: random seed for layout; default is not set (random)
  # width and height: width and height of visualization
  
  if (!requireNamespace("visNetwork", quietly=TRUE)) {
    stop("Package \"visNetwork\" needed for this function to work.")
  }

  use.aes <- match.arg(use.aes)
  label.value <- match.arg(label.value)
  layout <- match.arg(layout)
  if (is.null(rxns) && is.null(mets)) stop("Need to provide rxns and/or mets.")

  # mets
  if (is.null(mets)) met.ids <- unique(unlist(rxns2mets(model, rxn.ids))) else met.ids <- all2idx(model, mets)
  rm.mets <- get.exclude.mets(model, mets=NULL, rgx=exclude.mets.rgx, degree=ncol(model$S))
  met.ids <- setdiff(met.ids, rm.mets)
  mets <- model$mets[met.ids]
  tmp <- stringr::str_match(mets, "[\\[_](.)\\]?$")[,2]
  tmp <- paste0("[",stringr::str_sub(tmp,1,1),"]")
  met.ns <- paste0(model$metNames[met.ids], tmp)
  dup.mets <- get.exclude.mets(model, mets=NULL, rgx=dup.mets.rgx, degree=ncol(model$S))
  md.ids <- intersect(met.ids, dup.mets)

  # rxns
  if (is.null(rxns)) {
    rxn.ids <- unique(unlist(mets2rxns(model, met.ids)))
    tmp <- rep(TRUE, length(rxn.ids))
  } else {
    rxn.ids <- all2idx(model, rxns)
    tmp <- rxn.ids %in% unique(unlist(mets2rxns(model, met.ids)))
    rxn.ids <- rxn.ids[tmp]
  }
  # directions of reactions for visualization
  rv.rxns <- model$lb[rxn.ids]<0
  dirs <- as.numeric(!rv.rxns)
  # fluxes
  if (!is.null(fluxes)) {
    fluxes <- fluxes[tmp]
    fluxes[is.na(fluxes)] <- 0
    dirs[fluxes<0] <- -1
    v <- abs(fluxes)
    rxn.equs <- get.rxn.equations(model, rxn.ids, dir=fluxes, use.names=TRUE)
  } else {
    v <- NULL
    rxn.equs <- get.rxn.equations(model, rxn.ids, use.names=TRUE)
  }
  if (!is.null(dfluxes)) {
    dfluxes <- dfluxes[tmp]
    dfluxes[is.na(dfluxes)] <- 0
  }
  rxns <- model$rxns[rxn.ids]
  rxn.ns <- model$rxnNames[rxn.ids]
  
  # helper function of mapping line widths
  maplw <- function(x, lwds) {
    medx <- median(x)
    madx <- mad(x)
    lb <- medx - 2*madx
    ub <- medx + 2*madx
    x[x<lb] <- lb
    x[x>ub] <- ub
    if (uniqueN(x)==1) {
      res <- rep(1, length(x))
    } else {
      res <- x / diff(range(x)) * diff(lwds)
      res <- res - min(res) + lwds[1]
    }
    res
  }

  # colors and line widths
  if (!is.null(dfluxes)) {
    if (!abs.dflux) {
      if (is.null(fluxes)) {
        dirs[rv.rxns] <- sign(dfluxes[rv.rxns])
        dfluxes[rv.rxns] <- abs(dfluxes[rv.rxns])
      } else {
        dfluxes[rv.rxns] <- dfluxes[rv.rxns] * sign(fluxes[rv.rxns])
      }
    }
    if (!is.null(fluxes)) {
      cols <- map.colors(dfluxes, cols=cols, trim=TRUE, mid=0)
      lwds <- maplw(v, lwds)
    } else {
      if (use.aes %in% c("color","both")) cols <- map.colors(dfluxes, cols=cols, trim=TRUE, mid=0) else cols <- rep("#666666", length(rxns)) # grey40
      if (use.aes %in% c("width","both")) lwds <- maplw(abs(dfluxes), lwds) else lwds <- rep(8, length(rxns))
    }
  } else {
    if (!is.null(fluxes)) {
      if (use.aes %in% c("color","both")) cols <- map.colors(v, cols=cols, trim=TRUE, mid=0) else cols <- rep("#666666", length(rxns)) # grey40
      if (use.aes %in% c("width","both")) lwds <- maplw(v, lwds) else lwds <- rep(8, length(rxns))
    } else {
      cols <- rep("#1A1A1A", length(rxns)) # grey10
      lwds <- rep(2, length(rxns))
    }
  }

  # network data for visualization
  viz.dat <- lapply(1:length(rxn.ids), function(i) {
    # get mets in the reaction
    x <- model$S[met.ids, rxn.ids[i]]
    mi <- x!=0
    mets.i <- ifelse(met.ids[mi] %in% md.ids, paste0(mets[mi],i), mets[mi])
    rs <- mets.i[x[mi]<0] # reactants
    ps <- mets.i[x[mi]>0] # products
    # node (both mets and rxn) info
    nd <- data.table(id=c(mets.i, rxns[i]), label=c(mets[mi], rxns[i]),
                     title=c(met.ns[mi], sprintf("<p><b>%s</b><br>%s<br>v=%.4g<br>dv=%.4g</p>",rxn.ns[i],rxn.equs[i],ifelse(is.null(v[i]),NA,v[i]),ifelse(is.null(dfluxes[i]),NA,dfluxes[i]))),
                     group=c(rep("met", sum(mi)), "rxn"))
    # edge from reactants to reaction
    if (length(rs)!=0) {
      if (dirs[i]==1) ed <- data.table(from=rs, to=model$rxns[rxn.ids[i]], arrows="middle", smooth=TRUE, color=cols[i], width=lwds[i])
      if (dirs[i]==-1) ed <- data.table(from=model$rxns[rxn.ids[i]], to=rs, arrows="to", smooth=TRUE, color=cols[i], width=lwds[i])
      if (dirs[i]==0) ed <- data.table(from=rs, to=model$rxns[rxn.ids[i]], arrows="from;to", smooth=TRUE, color=cols[i], width=lwds[i])
    } else ed <- NULL
    # edge from reaction to products
    if (length(ps)!=0) {
      if (dirs[i]==-1) ed <- rbind(ed, data.table(from=ps, to=model$rxns[rxn.ids[i]], arrows="middle", smooth=TRUE, color=cols[i], width=lwds[i]))
      if (dirs[i]==1) ed <- rbind(ed, data.table(from=model$rxns[rxn.ids[i]], to=ps, arrows="to", smooth=TRUE, color=cols[i], width=lwds[i]))
      if (dirs[i]==0) ed <- rbind(ed, data.table(from=model$rxns[rxn.ids[i]], to=ps, arrows="from;to", smooth=TRUE, color=cols[i], width=lwds[i]))
    }
    # flux/dflux value label
    if (label.value=="flux" && !is.null(fluxes) || label.value=="dflux" && !is.null(dfluxes)) {
      lab <- switch(label.value, flux=sprintf("%.2g",v[i]), dflux=sprintf("%.2g",dfluxes[i]))
      ed[1, c("label","font.color"):=list(lab, cols[i])]
    }
    list(ed=ed, nd=nd)
  })

  # collect all edges
  eds <- rbindlist(lapply(viz.dat, function(x) x$ed))
  # collect all nodes
  nds <- unique(rbindlist(lapply(viz.dat, function(x) x$nd)))
  # draw network
  `%>%` <- visNetwork::`%>%`
  vis <- visNetwork::visNetwork(nds, eds, width=width, height=height) %>%
    visNetwork::visEdges(font=list(size=30)) %>%
    visNetwork::visGroups(groupname="met", shape="dot", size=15, color=list(border="#1A1A1A", background="#4169E1"), borderWidth=1.5, font=list(size=30, color="#000000")) %>% # grey10; background royalblue
    visNetwork::visGroups(groupname="rxn", shape="text", font=list(size=36, color="#000000")) %>%
    visNetwork::visOptions(highlightNearest=TRUE, nodesIdSelection=TRUE, collapse=TRUE)
  if (layout=="") vis <- vis %>% visNetwork::visLayout(randomSeed=seed) else vis <- vis %>% visNetwork::visIgraphLayout(layout=layout, randomSeed=seed)
  vis
}


plot.model1 <- function(model, rxns, fluxes=rep(1, length(rxns)), dfluxes=rep(0, length(rxns)), mets, exclude.mets.rgx="default", dup.mets.rgx=exclude.mets.rgx, use.flux=c("dflux","flux"), use=c("both","color","width"), cols=c("green4","grey","red3"), sizes=c(0.5,5), layout=c("neato","fdp","dot","circo","twopi"), margins=c(150,150,150,150)) {
  # plotting a metabolic model, can also incorporate fluxes and dfluxes data; this function doesn't work well in practice
  # model: the base metabolic model
  # rxns: reactions to plot
  # fluxes: the flux values corresponding to the reactions in rxns
  # dfluxes: the values of flux changes corresponding to the reactions in rxns
  # mets: metabolites to include in the plot
  # exclude.mets.rgx: regex for names of metabolites (as in model$mets) to be excluded from the plots; the default regex works for some of the high-degree mets in recon1 and iMM1415
  # dup.mets.rgx: after keeping the mets in mets.ids and excluding those in exclude.mets.rgx, for the remaining mets, use dup.mets to specify regex of mets to be plot as separate nodes for each reaction, when they are recurrent in multiple reactions; the default is the same as exclude.mets, so these will be excluded; to duplicate these instead of removing them, set exclude.mets to NULL
  # use.flux: choose to plot fluxes or dfluxes
  # use: use line color, or line width, or both to represent the flux or dflux values
  # cols: a vector of length 3, colors for decreased, unchanged, and increased reactions respectively if plotting dfluxes; if plotting fluxes, the 2nd and 3rd colors will be used for low and high fluxes respectively; if not using colors, the 2nd color will be used for all reactions
  # sizes: a vector of length 2, range of line widths
  # layout: graph layout for hyperdraw::graphLayout
  # margins: plot margins on the up, bottom, left, and right
  
  if (!requireNamespace(c("hypergraph","hyperdraw","RColorBrewer"), quietly=TRUE)) {
    stop("Packages \"hypergraph\", \"hyperdraw\" and \"RColorBrewer\" needed for this function to work.")
  }

  use.flux <- match.arg(use.flux)
  use <- match.arg(use)
  layout <- match.arg(layout)

  # build hyperedges from rxns
  rxn.ids <- all2idx(model, rxns)
  rxns <- model$rxns[rxn.ids]
  met.ids <- all2idx(model, mets)
  rm.mets <- get.exclude.mets(model, mets=NULL, rgx=exclude.mets.rgx, degree=ncol(model$S))
  met.ids <- setdiff(met.ids, rm.mets)
  mets <- model$mets[met.ids]
  dup.mets <- get.exclude.mets(model, mets=NULL, rgx=dup.mets.rgx, degree=ncol(model$S))
  md.ids <- intersect(met.ids, dup.mets)
  hypeds <- lapply(1:length(rxn.ids), function(i) {
    x <- model$S[met.ids, rxn.ids[i]]
    mi <- x!=0
    mets.i <- ifelse(met.ids[mi] %in% md.ids, paste0(mets[mi],i), mets[mi])
    rs <- mets.i[x[mi]<0] # reactants
    ps <- mets.i[x[mi]>0] # products
    if (length(rs)==0) rs <- paste0("EX_",ps)
    if (length(ps)==0) ps <- paste0("EX_",rs)
    # by default, create a hyperedge for a reaction, with arrows pointing from reactants to products:
    res <- hypergraph::DirectedHyperedge(rs, ps, label=rxns[i])
    # when plotting flux, if a reversible reaction is going backwards (i.e. flux<0), draw the arrow in the corresponding direction:
    if (use.flux=="flux" && fluxes[i]<0) res <- hypergraph::DirectedHyperedge(ps, rs, label=rxns[i])
    # when plotting dflux, (only) for the reversible reactions, draw the direction of the arrow according to the direction of dflux (also note: when using color, these will always be plotted in red or the color representing increase):
    if (use.flux=="dflux" && model$lb[rxn.ids][i]<0 && dfluxes[i]<0) res <- hypergraph::DirectedHyperedge(ps, rs, label=rxns[i])
    res
  })
  hypeds <- hypeds[!sapply(hypeds, is.null)]
  # now arrows for reversible reactions have been sorted out, we update flux and dflux values for plotting
  fluxes <- abs(fluxes)
  dfluxes[model$lb[rxn.ids]<0] <- abs(dfluxes[model$lb[rxn.ids]<0])

  # if plotting flux:
  if (use.flux=="flux") {
    v <- fluxes
    # trimming large flux values to facilitate plotting
    v0 <- median(v) + 2*mad(v)
    v[v>v0] <- v0
    # line colors:
    if (use %in% c("color","both")) {
      nc <- uniqueN(c(0,v))
      if (nc==1) {
        cols <- rep(cols[2],length(v))
      } else {
        cols <- colorRampPalette(cols[c(2,3)])(nc)
        bid <- as.numeric(cut(c(0,v), nc))
        cols <- cols[bid[-1]]
      }
    } else cols <- rep(cols[2],length(v))
    # line widths:
    if (use %in% c("width","both")) lwds <- v / (max(v)-min(v)) * diff(sizes) - min(v) + sizes[1] else lwds <- rep(2,length(v))
  }
  # if plotting dflux:
  if (use.flux=="dflux") {
    v <- dfluxes
    # if dfluxes values do not range from -1 and 1, trimm large positive and negative dflux values to facilitate plotting
    if (any(v>1) || any(v< -1)) {
      v0 <- median(v[v>0]) + 2*mad(v[v>0])
      v[v>v0] <- v0
      v0 <- median(v[v<0]) - 2*mad(v[v<0])
      v[v<v0] <- v0
    }
    # line colors:
    if (use %in% c("color","both")) {
      unqv <- unique(c(0,v))
      nc1 <- sum(unqv>=0)
      nc2 <- sum(unqv<=0)
      idx1 <- v>=0
      idx2 <- v<0
      if (nc1==1) {
        c1 <- rep(cols[2], sum(idx1))
      } else {
        c1 <- colorRampPalette(cols[2:3])(nc1)
        bid <- as.numeric(cut(c(0,v[idx1]), nc1))
        c1 <- c1[bid[-1]]
      }
      if (nc2==1) {
        c2 <- rep(cols[2], sum(idx2))
      } else {
        c2 <- colorRampPalette(cols[1:2])(nc2)
        bid <- as.numeric(cut(c(0,v[idx2]), nc2))
        c2 <- c2[bid[-1]]
      }
      cols <- c()
      cols[idx1] <- c1
      cols[idx2] <- c2
    } else cols <- rep(cols[2],length(v))
    # line widths:
    if (use %in% c("width","both")) {
      v <- abs(v)
      lwds <- v / (max(v)-min(v)) * diff(sizes) - min(v) + sizes[1]
    } else lwds <- rep(2,length(v))
  }

  # build graph object
  node.names <- unique(unlist(lapply(hypeds, function(x) c(x@head, x@tail))))
  hg <- hypergraph::Hypergraph(node.names, hypeds)
  testbph <- hyperdraw::graphBPH(hg)
  my.graph <- hyperdraw::graphLayout(testbph, layoutType=layout)
  # various plot parameters
  hyperdraw::nodeDataDefaults(my.graph, "shape") <- "box"
  hyperdraw::nodeDataDefaults(my.graph, "margin") <- 'unit(3, "mm")'  
  hyperdraw::edgeDataDefaults(my.graph, "lwd") <- 2
  hyperdraw::graphDataDefaults(my.graph, "arrowLoc") <- "end"
  # set line widths and colors
  for (i in 1:length(rxn.ids)) {
    rxn <- rxns[i]
    lwd <- as.character(lwds[i])
    col <- cols[i]
    lapply(my.graph@edgeNodeIO$outgoing[[rxn]], function(x) hyperdraw::edgeData(my.graph, rxn, x, "lwd") <- lwd)
    lapply(my.graph@edgeNodeIO$incoming[[rxn]], function(x) hyperdraw::edgeData(my.graph, x, rxn, "lwd") <- lwd)
    lapply(my.graph@edgeNodeIO$outgoing[[rxn]], function(x) hyperdraw::edgeData(my.graph, rxn, x, "color") <- col)
    lapply(my.graph@edgeNodeIO$incoming[[rxn]], function(x) hyperdraw::edgeData(my.graph, x, rxn, "color") <- col)
  }
  # plot margins
  my.graph@graph@boundBox@upRight@y <- my.graph@graph@boundBox@upRight@y + margins[1] # top
  my.graph@graph@boundBox@botLeft@y <- my.graph@graph@boundBox@botLeft@y - margins[2] # bottom
  my.graph@graph@boundBox@botLeft@x <- my.graph@graph@boundBox@botLeft@x - margins[3] # left 
  my.graph@graph@boundBox@upRight@x <- my.graph@graph@boundBox@upRight@x + margins[4] # right  
  
  # plot
  plot(my.graph)
  #return(my.graph)
}
