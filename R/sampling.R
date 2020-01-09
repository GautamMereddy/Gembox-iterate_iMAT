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
  message("Will generate ", n, " warmup points.")
  message("Begin generating warmup points...")
  orth.pnts <- get.orth.pnts(model, n, nc)
  rand.pnts <- get.rand.pnts(model, n, nc)
  r <- rep(runif(n), each=n.rxns)
  dim(r) <- c(n.rxns, n)
  res <- orth.pnts*r + rand.pnts*(1-r)
  message("Finished generating warmup points.")
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
  cl <- parallel::makeCluster(nc, type="FORK")
  res <- parallel::parApply(cl, mat, 2, get.opt.pnt, model=model)
  parallel::stopCluster(cl)
  res
}

get.rand.pnts <- function(model, n, nc) {
  # sample random points for ACHR
  # n: the number of random points
  # nc: number of cores to use

  n.rxns <- ncol(model$S)
  cs <- runif(n.rxns*n) - 0.5
  dim(cs) <- c(n.rxns, n)
  cl <- parallel::makeCluster(nc, type="FORK")
  res <- parallel::parApply(cl, cs, 2, get.opt.pnt, model=model)
  parallel::stopCluster(cl)
  res
}

get.opt.pnt <- function(model, c) {
  # a helper function to get an optimal point corresponding to running LP with objective function coefficients being c (normalized), for ACHR
  c <- c / norm(c,"2")
  res <- solve.model(model, csense="max", c=c, pars=.pkg.const$lp[[.pkg.var$solver]])[[1]]$xopt
}

