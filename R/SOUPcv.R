#' Cross validation for SOUP
#' 
#' @param expr A cell-by-gene expression matrix, either the raw counts or log-transformed expressions. 
#' @param type "log" if \code{expr} has been normalized and log-transformed (default),
#'     or "count" (default) \code{expr} contains the raw counts.
#' @param Ks A sequence of cluster numbers
#' @param nfold (optional) Number of folds, default is 10
#' @param nCV (optional) Number of repetitions, default is 10
#' @param mc.cores (optional) Number of cores for parallelization, default is 1 without parallelization
#' @param seeds (optional) A list of seeds to be used, with length nCV, default is NULL
#' @param verbose (optional) Whether to print progress, default is TRUE
#' 
#' @export
#' 
cvSOUP <- function(expr, type="log", Ks=c(2:10), 
                   nfold=10, nCV=10, mc.cores=1,
                   ext.prop=NULL, pure.prop=0.5,
                   seeds=NULL, verbose=TRUE) {
  
  cv.errors = matrix(NA, nrow=nCV, ncol=length(Ks))
  cv.sds = matrix(NA, nrow=nCV, ncol=length(Ks))

  for (i.cv in c(1:nCV)) {
    if (verbose) {
      cat("CV", i.cv, "...")
    }
    if (!is.null(seeds) && length(seeds)==nCV) {
      seed = seeds[i.cv]
    } else {
      seed = NULL
    }
    cv.out = cv.error.SOUP(expr=expr, type=type, 
                           seed=seed,
                           nfold=nfold, Ks=Ks,
                           ext.prop=ext.prop, pure.prop=pure.prop)
    cv.errors[i.cv, ] = cv.out$cvm
    cv.sds[i.cv, ] = cv.out$cvsd
  }
  
  if (verbose) {
    cat("done.\n")
  }
  
  cvm = colMeans(cv.errors)
  K.cv = Ks[which.min(cvm)]

  return(list(cv.errors = cv.errors,
              cv.sds = cv.sds,
              cvm = cvm,
              K.cv = K.cv))
}

#' Compute cross validation errors
#' 
#' @param expr A cell-by-gene expression matrix, either the raw counts or log-transformed expressions. 
#' @param type "log" if \code{expr} has been normalized and log-transformed (default),
#'     or "count" (default) \code{expr} contains the raw counts.
#' @param Ks A sequence of cluster numbers
#' @param nfold (optional) Number of folds, default is 10
#' @param seed (optional) random seed, default is NULL
#' @param mc.cores (optional) Number of cores for parallelization, default is 1 without parallelization
#' 
#' @export
cv.error.SOUP <- function(expr, type="log", Ks=c(2:10), 
                          nfold=10, seed=NULL, mc.cores=1,
                          ext.prop=NULL, pure.prop=0.5) {
  
  ## cross validation
  doCV <- function(fold, nfold, i.permute.ind, 
                   expr, Ks, type, ext.prop, pure.prop) {
    n.sc = nrow(expr)
    n.gene = ncol(expr)
    n.K = length(Ks)
    fold.size = floor(n.sc / nfold)
    
    ## test set
    i.start = (fold.size*(fold-1)+1)
    i.end = min(fold.size*fold, n.sc)
    i.test = i.permute.ind[c(i.start:i.end)]
    n.test = length(i.test)
    
    ## SOUP on training set
    soup.out = SOUP(expr=expr[-i.test, ], Ks=Ks, type=type,
                    ext.prop=ext.prop, pure.prop=pure.prop,
                    nstart=50, verbose=FALSE)
    
    ## Predict on test set
    cv.error = rep(NA, n.K)
    for (i.K in c(1:n.K)) {
      train.centers = soup.out$centers[[i.K]]
      test.theta = predictTheta(new.expr=expr[i.test, ], 
                                t.centers=t(train.centers))
      test.pred = test.theta %*% train.centers
      cv.error[i.K] = sum((expr[i.test, ] - test.pred)^2) / (n.test*n.gene)
    }
    
    return(cv.error)
  }
  
  ## permute indices
  n.sc = nrow(expr)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  i.permute.ind = base::sample(c(1:n.sc), size=n.sc, replace=FALSE)
  
  ## cross validation in parallel
  cv.errors = parallel::mclapply(c(1:nfold), FUN=doCV, 
                                   nfold=nfold, i.permute.ind=i.permute.ind, 
                                   expr=expr, Ks=Ks, type=type,
                                   ext.prop=ext.prop, pure.prop=pure.prop,
                                   mc.cores=mc.cores, mc.preschedule=TRUE)
  cv.errors = matrix(unlist(cv.errors), nrow=nfold, byrow=TRUE)
  
  return(list(cvm = colMeans(cv.errors),
              cvsd = apply(cv.errors, 2, sd)
  ))
}


#' Predict the membership for new data points
#' 
#' @param new.expr cell-by-gene expression matrix
#' @param t.centers transposed center matrix, n.gene-by-K
#' 
#' @return The predicted membership matrix.
#' @export
#' 
predictTheta <- function(new.expr, t.centers) {
  n.cell = nrow(new.expr)
  K = ncol(t.centers)
  
  ## solve for t(new.expr) = t(centers) %*% t(theta)
  theta = matrix(0, nrow=n.cell, ncol=K)
  for (i in c(1:n.cell)) {
    theta[i, ] = limSolve::lsei(A = t.centers, 
                                B = new.expr[i, ], 
                                type=2)$X
  }
  row.names(theta) <- row.names(new.expr)
  
  ## clean up
  theta = projMembership(theta)
  
  return(theta)
}



