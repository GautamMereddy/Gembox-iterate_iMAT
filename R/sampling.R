###### functions for sampling from a model ######


sample.model <- function(model, pars=list()) {
  # sample a model
  # pars: parameters for sampling in a list
  pars <- get.pars("samp", pars)

  if (pars$method=="achr") {
    if ("sample" %in% names(model)) {
      message("Will use the warmup points and status stored in the model.")
      message("Will sample ", pars$n.sample, " points.")
      res <- achr(model, model$sample$stat, model$sample$warmup.pnts, pars$n.sample, pars$steps.per.pnt)
      model$sample$stat <- res$stat
      model$sample$pnts <- cbind(model$sample$pnts, res$pnts)
    } else {
      warmup.pnts <- sample.warmup.pnts(model, pars$n.warmup, pars$nc)
      center.pnt <- rowMeans(warmup.pnts)
      init.stat <- list(center.pnt=center.pnt, prev.pnt=center.pnt, n.tot.steps=0)
      message("Will sample ", pars$n.sample, " points.")
      res <- achr(model, init.stat, warmup.pnts, pars$n.sample, pars$steps.per.pnt)
      model$sample <- list()
      model$sample$warmup.pnts <- warmup.pnts
      model$sample$pnts <- res$pnts
      model$sample$stat <- res$stat
    }
  }
  model
}

sample.warmup.pnts <- function(model, n, nc) {
  # sample warmup points for ACHR
  # n: the number of warmup points
  # nc: number of cores to use

  n.rxns <- ncol(model$S)
  if (n<2*n.rxns) {
    n <- 2*n.rxns
    warning(sprintf("#{warmup points} should be at least 2*#{reactions}=%d.", 2*n.rxns), immediate.=TRUE)
  }
  message("Begin generating ", n, " warmup points...")
  orth.pnts <- get.orth.pnts(model, n, nc)
  rand.pnts <- get.rand.pnts(model, n, nc)
  r <- rep(runif(n), each=n.rxns)
  dim(r) <- c(n.rxns, n)
  res <- orth.pnts*r + rand.pnts*(1-r)
  message("Done generating warmup points.")
  res
}

get.orth.pnts <- function(model, n, nc) {
  # sample orthogonal points for ACHR
  # n: the number of orthogonal points
  # nc: number of cores to use

  n.rxns <- ncol(model$S)
  mat <- cbind(Diagonal(n.rxns), Diagonal(n.rxns, x=-1))
  if (n<=2*n.rxns) {
    mat <- mat[, sample(2*n.rxns, n)]
  } else {
    mat <- cbind(mat[, sample(2*n.rxns)], mat[, sample(2*n.rxns, n-2*n.rxns, replace=TRUE)])
  }
  message("1. Generating orthogonal points, progress:")
  get.opt.pnts(model, mat, nc)
}

get.rand.pnts <- function(model, n, nc) {
  # sample random points for ACHR
  # n: the number of random points
  # nc: number of cores to use

  n.rxns <- ncol(model$S)
  cs <- runif(n.rxns*n) - 0.5
  dim(cs) <- c(n.rxns, n)
  message("2. Generating random points, progress:")
  get.opt.pnts(model, cs, nc)
}

get.opt.pnts <- function(model, mat, nc) {
  # a helper function to get optimal points corresponding to running LPs with objective function coefficients being the columns of mat, for ACHR

  do.call(cbind, pbmcapply::pbmclapply(1:ncol(mat), function(i) solve.model(model, csense="max", c=mat[,i]/norm(mat[,i],"2"), pars=.pkg.const$lp[[.pkg.var$solver]])[[1]]$xopt, mc.cores=nc))
}
