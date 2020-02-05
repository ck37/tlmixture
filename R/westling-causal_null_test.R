# All code in this file via Ted Westling. Minor edits were made for package checking compatability.

# Y: n x 1 numeric outcome vector
# A: n x 1 numeric exposure vector
# W: n x p data.frame of confounders
# g.hats: n x 1 numeric vector of estimated  G(da | w) / F(da) values (can be obtained using con.mixed.dens.SL)
# mu.hat: function taking arguments a, w, and returning estimated outcome regression value.
# p: numeric vector of powers p to use for the test, i.e. p=c(1,2,Inf)
# alpha: testing level
# n.sim: number of simulations to use for limit Gaussian process estimation
# return.Omega: whether to return the estimated primitive function and confidence interval/band
# cv.folds: if k-fold cross-fitting was used for g.hat and mu.hat, then cv.folds should be a length-k list of indices in each fold
#' @import mvtnorm
#' @importFrom sets interval_contains_element interval interval_complement
causal.null.test <- function(Y, A, W, g.hats, mu.hat, p=2, alpha = .05, n.sim = 1e4, return.Omega = FALSE, cv.folds = NULL) {
#  library(sets)
#  library(mvtnorm)
  n <- length(Y)
  a.vals <- sort(unique(A))

  if(is.null(cv.folds)) {
    ord <- order(A)
    A <- A[ord]
    Y <- Y[ord]
    W <- W[ord,]
    g.hats <- g.hats[ord]
    a.ecdf <- stats::ecdf(A)
    a.weights <- sapply(a.vals, function(a0) mean(A == a0))
    A.a.val <- sapply(A, function(a0) which(a.vals == a0))
    u.vals <- a.ecdf(a.vals)
    mu.hats.a.vals <- sapply(a.vals, function(a0) mu.hat(a0, W)) #rows index W, columns index a.vals
    mu.hats <- mu.hats.a.vals[,A.a.val]
    #g.hats <- c(g.hat(A, W))
    theta.a.vals <- colMeans(mu.hats.a.vals)
    theta.A <- theta.a.vals[A.a.val]

    mu.hats.data <- diag(mu.hats)

    partial.mu.means <- t(apply(mu.hats, 1, cumsum)) / n

    gamma.hat <- mean(mu.hats)
    #gamma.hat <- mean(Y)#mean(diag(mu.hats))#

    Omega.a.vals <- sapply(a.vals, function(a0) mean(as.numeric(A <= a0) * theta.A)) - gamma.hat * u.vals

    IF.vals <- sapply(a.vals, function(a0) {
      #mumean.vals <- apply(mu.hats, 1, function(row) mean(as.numeric(A <= a0) * row))
      mumean.vals <- partial.mu.means[,max(which(A <= a0))]
      #as.numeric(A <= a0) * ((Y - mu.hats.data) / g.hats + theta.A - gamma.hat) + mumean.vals - Y * a.ecdf(a0) - (mean(as.numeric(A <= a0) * theta.A) - gamma.hat * a.ecdf(a0))
      (as.numeric(A <= a0) - a.ecdf(a0)) * ((Y - mu.hats.data) / g.hats + theta.A - gamma.hat) + mumean.vals - partial.mu.means[,n] * a.ecdf(a0) - 2 * Omega.a.vals[which(a.vals == a0)]
    })

    Omega.hat <- colMeans(IF.vals) + Omega.a.vals

    # plot(a.vals, Omega.hat, type='l', ylim=c(-.2, .5))
    # lines(a.vals, Omega.a.vals, col='blue')
    # lines(a0.vals, Omega0.vals, col='red')

    Sigma.hat <- sapply(1:length(a.vals), function(s) sapply(1:length(a.vals), function(t) {
      mean(IF.vals[,s] * IF.vals[,t])
    }))
    #paths <- matrix(rnorm(n * n.sim), ncol=n) %*% IF.vals

    # sds <- sqrt(colMeans(IF.vals.cent^2))
    # lines(a.vals, Omega.hat + quantile(apply(abs(paths), 1, max), .975) / sqrt(n) , col='red')
    # lines(a.vals, Omega.hat - quantile(apply(abs(paths), 1, max), .975) / sqrt(n) , col='red')

  } else {
    n.folds <- length(cv.folds)
    fold.Omega.hats <- matrix(NA, nrow = n.folds, ncol = length(a.vals))
    IF.vals <- vector(length=n.folds, mode='list')
    for(j in 1:n.folds) {
      Nv <- length(cv.folds[[j]])
      A.test <- A[cv.folds[[j]]]
      Y.test <- Y[cv.folds[[j]]]
      W.test <- W[cv.folds[[j]],]
      ord <- order(A.test)
      A.test <- A.test[ord]
      Y.test <- Y.test[ord]
      W.test <- W.test[ord,]
      g.hats.test <- g.hats[[j]][ord]
      a.ecdf <- stats::ecdf(A.test)
      a.weights <- sapply(a.vals, function(a0) mean(A.test == a0))
      A.a.val <- sapply(A.test, function(a0) which(a.vals == a0))
      u.vals <- a.ecdf(a.vals)
      mu.hats.a.vals <- sapply(a.vals, function(a0) mu.hat[[j]](a0, W.test)) #rows index W, columns index a.vals
      mu.hats <- mu.hats.a.vals[,A.a.val]
      #g.hats <- c(g.hat(A, W))
      theta.a.vals <- colMeans(mu.hats.a.vals)
      theta.A <- theta.a.vals[A.a.val]

      mu.hats.data <- diag(mu.hats)

      partial.mu.means <- t(apply(mu.hats, 1, cumsum)) / Nv

      gamma.hat <- mean(mu.hats)
      #gamma.hat <- mean(Y)#mean(diag(mu.hats))#

      Omega.a.vals <- sapply(a.vals, function(a0) mean(as.numeric(A.test <= a0) * theta.A)) - gamma.hat * u.vals

      IF.vals[[j]] <- sapply(a.vals, function(a0) {
        #mumean.vals <- apply(mu.hats, 1, function(row) mean(as.numeric(A <= a0) * row))
        mumean.vals <- partial.mu.means[,max(which(A.test <= a0))]
        #as.numeric(A <= a0) * ((Y - mu.hats.data) / g.hats + theta.A - gamma.hat) + mumean.vals - Y * a.ecdf(a0) - (mean(as.numeric(A <= a0) * theta.A) - gamma.hat * a.ecdf(a0))
        (as.numeric(A.test <= a0) - a.ecdf(a0)) * ((Y.test - mu.hats.data) / g.hats.test + theta.A - gamma.hat) + mumean.vals - partial.mu.means[,ncol(partial.mu.means)] * a.ecdf(a0) - 2 * Omega.a.vals[which(a.vals == a0)]
      })

      fold.Omega.hats[j,] <- colMeans(IF.vals[[j]]) + Omega.a.vals
    }

    Omega.hat <- colMeans(fold.Omega.hats)
    Sigma.hat <- sapply(1:length(a.vals), function(s) sapply(1:length(a.vals), function(t) {
      mean(unlist(lapply(IF.vals, function(IF) mean(IF[,s] * IF[,t]))))
    }))
    #paths <- matrix(rnorm(n * n.sim), ncol=n) %*% IF.vals
  }

  paths <- mvtnorm::rmvnorm(n.sim, sigma=Sigma.hat)

  a.weights <- sapply(a.vals, function(a) mean(A == a))
  ret <- t(sapply(p, function(pp) {
    stat <- ifelse(pp < Inf, (sum(abs(Omega.hat )^pp * a.weights))^{1/pp}, max(abs(Omega.hat)))

    if(pp < Inf) {
      stats <- (apply(abs(paths)^pp, 1, function(row) sum(row * a.weights)))^{1/pp}
    } else {
      stats <- apply(abs(paths), 1, max)
    }

    p.val <- mean(stats / sqrt(n) > stat)

    q <- quantile(stats, (1 - alpha / 2))
    ci.ll <- max(stat - q / sqrt(n), 0)
    ci.ul <- stat + q / sqrt(n)

    res <- c(stat, p.val, ci.ll, ci.ul)
    res
  }))
  ret.df <- data.frame(p = p, obs.stat = ret[,1], p.val = ret[,2])
  if(!return.Omega) return(ret.df)
  if(return.Omega) {
    ret.df$ci.ll <- ret[,3]
    ret.df$ci.ul <- ret[,4]
    ret.list <- list(test = ret.df, Omega.hat = Omega.hat, IF.vals = IF.vals, paths = paths)
    return(ret.list)
  }
}


find.bin <- function(x, bins) {
  mat <- t(sapply(x-1e-10, function(x0) {
    unlist(lapply(bins, function(bin) {
      sets::interval_contains_element(bin, x0)
    }))
  }))
  if(any(rowSums(mat) > 1)) stop("Overlapping bins")
  if(any(rowSums(mat) == 0)) stop("Element outside all bins")

  apply(mat, 1, function(row) which(row))
  # ret <- rep(NA, length(x))
  # if(any(bins$mass.pt)) {
  #   for(row in which(bins$mass.pt)) {
  #     ret[ x == bins$lower[row] ] <- bins$bin[row]
  #   }
  # }
  # breaks <- sort(c(bins$lower[!bins$mass.pt], max(bins$upper[!bins$mass.pt])))
  # ret[is.na(ret)] <- findInterval(x[is.na(ret)], vec = breaks, all.inside = TRUE, left.open = TRUE)
  # return(ret)
}


# A: n x 1 numeric exposure vector
# W: n x p data.frame of confounders
# n.bins: numeric vector of number of bins to use, i.e. 2:5 or c(2,4,6,8,10)
# SL.library: super learner library to use for each bin-specific estimate
# verbose: whether to report progress
# n.folds: number of folds for cross-validation across bins
#' @import Rsolnp
con.mixed.dens.SL <- function(A, W, n.bins, SL.library, verbose=FALSE, n.folds = 10) {
  n <- nrow(W)

  sorted  <- rep(1:n.folds, length.out = n)
  folds <- sample(sorted, n, replace=FALSE)
  valid.rows <- lapply(1:n.folds, function(v) which(folds == v))

  #library(Rsolnp)

  fits <- NULL
  for(b in n.bins) {
    if(verbose) cat("\nEstimating models with", b, "bins... ")
    fits[[paste0('dens.fit.', b, 'bins')]] <- con.mixed.dens.one.bin(A, W, n.bins = b, SL.library = SL.library, verbose = verbose, valid.rows = valid.rows, n.folds = n.folds)
  }

  algs.per.bin <- ncol(fits[[1]]$cv.library.densities)
  n.algs <- length(n.bins) * algs.per.bin
  cv.library.densities <- library.densities <- matrix(NA, nrow=n, ncol=n.algs)
  library.names <- NULL
  start.col <- 1
  for(b in n.bins) {
    end.col <- start.col + algs.per.bin - 1
    cv.library.densities[,start.col:end.col] <- fits[[paste0('dens.fit.', b, 'bins')]]$cv.library.densities
    library.densities[,start.col:end.col] <- fits[[paste0('dens.fit.', b, 'bins')]]$library.densities
    library.names <- c(library.names, fits[[paste0('dens.fit.', b, 'bins')]]$alg.names)
    start.col <- end.col + 1
  }

  if(verbose) cat("\nOptimizing model weights...\n")

  # Remove algs with errors in cv predictions
  errors.in.library <- apply(cv.library.densities, 2, function(col) any(is.na(col)))
  if(any(errors.in.library)) warning(paste0("Errors in the following candidate algorithms: ", library.names[which(errors.in.library)]))
  n.include <- sum(!errors.in.library)

  # Do SL log-likelihood optimization
  cv_risk <- function(beta) -mean(log(cv.library.densities[,!errors.in.library] %*% beta))
  utils::capture.output(solnp_solution <- solnp(rep(1/n.include, n.include), cv_risk, eqfun=sum, eqB=1, ineqfun=function(beta) beta, ineqLB=rep(0,n.include), ineqUB=rep(1, n.include)))
  coef <- rep(0, n.algs)
  coef[!errors.in.library] <- solnp_solution$pars
  if(verbose) {
    cat("Top five learners by weight: \n")
    for(j in 1:5) {
      cat(library.names[order(coef, decreasing = TRUE)[j]], " (weight ", sort(coef, decreasing = TRUE)[j], ")\n", sep='')
    }
  }
  SL.density <- c(library.densities[,!errors.in.library] %*% solnp_solution$pars)

  return(list(n.bins = n.bins, fits = fits, cv.library.densities = cv.library.densities, library.densities = library.densities, SL.densities = SL.density, coef = coef, library.names = library.names, a.ecdf = stats::ecdf(A)))


}

#' @import Rsolnp
#' @import SuperLearner
# @import sets
con.mixed.dens.one.bin <- function(A, W, n.bins, SL.library, verbose=FALSE, valid.rows=NULL, n.folds=10) {
  a.ecdf <- stats::ecdf(A)
  U <- a.ecdf(A)
  n <- nrow(W)
  W <- as.data.frame(W)
  U <- as.numeric(U)
  #library(Rsolnp)
  #library(SuperLearner)
  #library(sets)
  tab <- table(U)
  un.U <- as.numeric(names(tab))
  un.U.frac <- as.numeric(tab) / length(U)
  if(n.bins <= 1) stop("Number of bins must be > 1")
  if(length(un.U) < n.bins) stop("Number of bins must not be larger than number of unique values of U.")
  if(length(un.U) == n.bins) {
    mass.pts <- un.U
    bins <- data.frame(bin = 1:n.bins, lower = un.U, upper = un.U, bin.length = 0, mass.pt = TRUE)
  }
  if(length(un.U) > n.bins) {
    if(any(un.U.frac >= 1/n.bins)) {
      mass.pts <- un.U[un.U.frac >= 1/n.bins]
      n.mass.pts <- length(mass.pts)
      mass.pt.lowers <- sapply(mass.pts, function(x) max(c(U[x - U > 1/(10*n)], 0)))
      mass.intervals <- lapply(1:n.mass.pts, function(j) {
        sets::interval(mass.pt.lowers[j], mass.pts[j], bounds="(]")
      })
      cont.intervals <- data.frame(lower=c(0,mass.pts), upper=c(mass.pt.lowers, 1))
      cont.intervals$length <- cont.intervals$upper - cont.intervals$lower
      cont.intervals <- subset(cont.intervals, length > 0)
    } else {
      mass.pts <- NULL
      n.mass.pts <- 0
      mass.intervals <- NULL
      cont.intervals <- data.frame(lower=0, upper=1, length=1)
    }

    n.cont.bins <- n.bins - n.mass.pts
    delta <- sum(cont.intervals$length) / n.cont.bins
    delta <- round(delta, digits=ceiling(log10(n)) + 2)
    cont.bin.endpts <- matrix(NA, nrow=n.cont.bins, ncol=2)

    for(j in 1:n.cont.bins) {
      if(j == 1) start <- cont.intervals$lower[1]
      else start <- end
      start.interval <- max(which(cont.intervals$lower <= start & start <= cont.intervals$upper))
      if(start == cont.intervals$upper[start.interval]) {
        start.interval <- start.interval + 1
        start <- cont.intervals$lower[start.interval]
      }
      end <- start + delta
      end.interval <- start.interval
      if(!(all.equal(end, cont.intervals$upper[end.interval]) == TRUE) && end > cont.intervals$upper[end.interval]) {
        length.used <- cont.intervals$upper[end.interval] - start
        length.left <- delta - length.used
        end.interval <- end.interval + 1
        end <- cont.intervals$lower[end.interval] + length.left
      }
      while(!(all.equal(end, cont.intervals$upper[end.interval]) == TRUE) && end > cont.intervals$upper[end.interval]) {
        length.used <- length.used + cont.intervals$upper[end.interval] - cont.intervals$lower[end.interval]
        length.left <- delta - length.used
        end.interval <- end.interval + 1
        end <- cont.intervals$lower[end.interval] + length.left
      }
      end <- round(end, digits=ceiling(log10(n)) + 3)
      if(j == n.cont.bins) end <- cont.intervals$upper[nrow(cont.intervals)]
      cont.bin.endpts[j,] <- c(start, end)
    }

    cont.intervals <- lapply(1:n.cont.bins, function(j) {
      if(j == 1) int <- sets::interval(cont.bin.endpts[j, 1], cont.bin.endpts[j, 2], bounds="[]")
      else int <- sets::interval(cont.bin.endpts[j, 1], cont.bin.endpts[j, 2], bounds="(]")
      if(n.mass.pts > 0) {
        for(k in 1:n.mass.pts) {
          int <- sets::interval_complement(mass.intervals[[k]], int)
        }
      }
      return(int)
    })

    bins <- c(mass.intervals, cont.intervals)
    bin.sizes <- unlist(lapply(bins, interval_measure))

    # ycuts <- seq(0,1,by=1/n.cont.bins)
    # bin.cuts <- as.numeric(quantile(U[!(U %in% mass.pts)], ycuts, type=1))
    # bins <- data.frame(bin = 1:n.cont.bins , lower = bin.cuts[1:n.cont.bins], upper = bin.cuts[-1])
    # bins$bin.length <- bins$upper - bins$lower
    # bins$mass.pt <- FALSE
    # if(any(bins$bin.length == 0)) {
    #   bins$mass.pt[bins$bin.length == 0] <- TRUE
    # }
    # if(n.mass.pts > 0) {
    #   mass.pt.bins <- data.frame(bin=(n.cont.bins + 1):n.bins, lower = mass.pts, upper = mass.pts, bin.length = 0, mass.pt = TRUE)
    #   bins <- rbind(bins, mass.pt.bins)
    # }
  }

  # bins$weight <- ifelse(bins$mass.pt, 1 / sapply(bins$lower, function(l) mean(U == l)), 1 / bins$bin.length)
  #
  disc.U <- find.bin(U, bins)
  #disc.U.num <- as.numeric(disc.U)
  #labs <- levels(disc.U)

  bin.fits <- NULL
  bin.probs <- matrix(NA, nrow=n, ncol=n.bins)
  for(bin in 1:n.bins) {
    if(verbose) cat("bin", bin, "... ")
    utils::capture.output(bin.fit <- try(SuperLearner(Y=as.numeric(disc.U==bin), X=W, family='binomial', SL.library = SL.library, method='method.NNloglik', cvControl = list(V=n.folds, validRows=valid.rows)), silent=TRUE))
    if(class(bin.fit) == "try-error") {
      utils::capture.output(bin.fit <- try(SuperLearner(Y=as.numeric(disc.U==bin), X=W, family='binomial', SL.library = SL.library, method='method.NNLS', cvControl = list(V=n.folds, validRows=valid.rows)), silent=TRUE))
    }
    if(class(bin.fit) == "try-error") {
      utils::capture.output(bin.fit <- try(SuperLearner(Y=as.numeric(disc.U==bin), X=W, family='binomial', SL.library = SL.library, method='method.NNLS2', cvControl = list(V=n.folds, validRows=valid.rows)), silent=TRUE))
    }
    if(class(bin.fit) != "try-error") {
      bin.fits[[paste0("bin", bin, ".SL")]] <- bin.fit
      bin.probs[,bin] <- bin.fit$SL.predict
    } else {
      bin.mean <- mean(as.numeric(disc.U==bin))
      if(class(SL.library) == "character") n.algs <- length(SL.library)
      else n.algs <- sum(unlist(lapply(SL.library, function(sl) length(sl) - 1)))
      bin.fits[[paste0("bin", bin, ".SL")]] <- list(Z = matrix(bin.mean, nrow=n,ncol=n.algs), library.predict = matrix(bin.mean, nrow=n,ncol=n.algs))
      bin.probs[,bin] <- bin.mean
    }
  }

  SL.bin.probs <- make.doubly.stochastic(bin.probs, row.sums = rep(1, n), col.sums = bin.sizes * n)#bin.probs / rowSums(bin.probs)

  SL.densities <- SL.bin.probs[cbind(1:n, disc.U)] / bin.sizes[disc.U]

  n.alg <- ncol(bin.fits[["bin1.SL"]]$Z)
  cv.library.densities <- library.densities <- matrix(NA, nrow=n, ncol=n.alg)
  for (j in 1:n.alg) {
    cv.bin.probs <- library.bin.probs <- matrix(NA, nrow=n, ncol=n.bins)
    for (bin in 1:n.bins) {
      cv.bin.probs[, bin] <- bin.fits[[paste0("bin", bin, ".SL")]]$Z[,j]
      library.bin.probs[, bin] <- bin.fits[[paste0("bin", bin, ".SL")]]$library.predict[,j]
    }
    if(any(is.na(cv.bin.probs)) | any(colSums(cv.bin.probs) == 0) | any(rowSums(cv.bin.probs) == 0)) {
      cv.library.densities[,j] <- rep(NA, n)
    } else {
      cv.bin.probs <- make.doubly.stochastic(cv.bin.probs, row.sums = rep(1, n), col.sums = bin.sizes * n)
      cv.library.densities[,j] <- cv.bin.probs[cbind(1:n, disc.U)] / bin.sizes[disc.U]
    }

    if(any(is.na(library.bin.probs)) | any(colSums(library.bin.probs) == 0) | any(rowSums(library.bin.probs) == 0)) {
      library.densities[,j] <- rep(NA, n)
    } else {
      library.bin.probs <- make.doubly.stochastic(library.bin.probs, row.sums = rep(1, n), col.sums = bin.sizes * n)
      library.densities[,j] <- library.bin.probs[cbind(1:n, disc.U)] / bin.sizes[disc.U]
    }

  }

  alg.names <- paste0(bin.fits[["bin1.SL"]]$libraryNames, "_", n.bins, "bins")

  return(list(bins = bins, bin.fits = bin.fits, a.ecdf = a.ecdf, SL.bin.probs = SL.bin.probs, SL.densities = SL.densities, cv.library.densities = cv.library.densities, library.densities = library.densities, alg.names = alg.names))

}



make.doubly.stochastic <- function(mat, row.sums, col.sums, tol = .001) {
  ret <- mat
  while(sum(abs(rowSums(ret) - row.sums)) > tol | sum(abs(colSums(ret) - col.sums)) > tol) {
    ret <- ret / (rowSums(ret) / row.sums)
    ret <- t( t(ret) / (colSums(ret) / col.sums))
  }
  return(ret)
}


predict.con.mixed.dens.SL <- function(fit, new.A, new.W, threshold = .001) {
  new.W <- as.data.frame(new.W)
  new.U <- fit$a.ecdf(new.A)
  trunc.coef <- fit$coef
  trunc.coef[trunc.coef < threshold] <- 0
  trunc.coef <- trunc.coef / sum(trunc.coef)
  nonzero <- which(trunc.coef > 0)
  lib.name.splits <- strsplit(fit$library.names, "_")
  lib.name.nbins <- unlist(lapply(lib.name.splits, function(l) as.numeric(strsplit(l[3], "bins")[[1]])))
  lib.name.alg <- unlist(lapply(lib.name.splits, function(l) paste0(l[1:2], collapse="_")))
  bins.to.fit <- unique(lib.name.nbins[nonzero])
  pred.densities <- matrix(NA, nrow=length(new.A), ncol=length(fit$library.names))
  for(bin in bins.to.fit) {
    ind <- which(fit$n.bins == bin)
    new.bins <- find.bin(new.U, bins = fit$fits[[ind]]$bins)
    bin.sizes <- unlist(lapply(fit$fits[[ind]]$bins, interval_measure))
    pred.probs <- matrix(NA, nrow = length(new.U), ncol = length(unique(lib.name.alg)))
    for(k in 1:length(fit$fits[[ind]]$bin.fits)) {
      if(any(new.bins == k)) {
        pred.probs[new.bins == k,] <- predict.SuperLearner(fit$fits[[ind]]$bin.fits[[k]], newdata = new.W[new.bins == k,])$library.predict
      }
    }
    pred.densities[,which(lib.name.nbins == bin)] <- pred.probs / bin.sizes[new.bins]
  }
  c(pred.densities[,nonzero] %*% trunc.coef[nonzero])
}


