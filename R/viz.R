###### functions for visualization of metabolic models ######


plot.model <- function(model, rxns, fluxes=rep(1, length(rxns)), dfluxes=rep(0, length(rxns)), mets, exclude.mets.rgx="default", dup.mets.rgx=exclude.mets.rgx, use.flux=c("dflux","flux"), use=c("both","color","width"), cols=c("green4","grey","red3"), sizes=c(0.5,5), layout=c("neato","fdp","dot","circo","twopi"), margins=c(150,150,150,150)) {
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
        cols <- RColorBrewer::colorRampPalette(cols[c(2,3)])(nc)
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
        c1 <- RColorBrewer::colorRampPalette(cols[2:3])(nc1)
        bid <- as.numeric(cut(c(0,v[idx1]), nc1))
        c1 <- c1[bid[-1]]
      }
      if (nc2==1) {
        c2 <- rep(cols[2], sum(idx2))
      } else {
        c2 <- RColorBrewer::colorRampPalette(cols[1:2])(nc2)
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

check.mets.plot <- function(model, fluxes=rep(1, length(rxn.ids)), dfluxes=rep(0, length(rxn.ids)), mets, exclude.mets.rgx="default") {
  # generate an interactive network plot involving a set of metabolites given in mets
  # the plot will show the interconversion among the mets, and the interconversion rates can be plotted based on data in fluxes or dfluxes
  # this function is not fully implemented yet...
  if (!requireNamespace("visNetwork", quietly=TRUE)) {
    stop("Package \"visNetwork\" needed for this function to work.")
  }

  s <- model$S
  rm.mets <- get.exclude.mets(model, mets=NULL, rgx=exclude.mets.rgx, degree=ncol(model$S))
  s[rm.mets, ] <- 0
  mets <- all2idx(model, mets)
  rxn.ids <- unique(unlist(mets2rxns(model, mets)))
  # network data
  tmp <- lapply(rxn.ids, function(i) {
    from <- which(s[,i]<0)
    to <- which(s[,i]>0)
    # edge info
    ed <- as.data.table(expand.grid(from=from, to=to))[, c("title","arrows","smooth"):=list(model$rxns[i], "to", TRUE)]
    # if only 2 or 3 reactants/products, also pull them together by edges
    if (length(from)==2) ed <- rbind(ed, data.table(from=from[1], to=from[2], title=model$rxns[i], arrows=NA, smooth=TRUE))
    if (length(from)==3) ed <- rbind(ed, as.data.table(expand.grid(from=from, to=from))[c(2,3,6)][, c("title","arrows","smooth"):=list(model$rxns[i], NA, TRUE)])
    if (length(to)==2) ed <- rbind(ed, data.table(from=to[1], to=to[2], title=model$rxns[i], arrows=NA, smooth=TRUE))
    if (length(to)==3) ed <- rbind(ed, as.data.table(expand.grid(from=to, to=to))[c(2,3,6)][, c("title","arrows","smooth"):=list(model$rxns[i], NA, TRUE)])
    # node info
    nd <- data.table(id=c(from,to), label=model$mets[c(from,to)])
    list(ed=ed, nd=nd)
  })
  # collect all edges
  eds <- rbindlist(lapply(tmp, function(x) x$ed))
  # collapse
  # collect all nodes
  nds <- rbindlist(lapply(tmp, function(x) x$nd))
  nds <- unique(nds)
  nds[, group:=ifelse(id %in% met.ids, "main", "adjacent")]
  visNetwork::visNetwork(nds, eds) %>% visNetwork::visIgraphLayout(layout="layout_with_fr")
}
