#' Cross Validation for SOUP
#' 
#' @export
#' 
cvSOUP <- function(expr, type="count", 
                   nfold=10, nCV=10, Ks=c(2:10)) {
  
  cv.errors = matrix(NA, nrow=nCV, ncol=length(Ks))
  cv.sds = matrix(NA, nrow=nCV, ncol=length(Ks))
  cv.aics = matrix(NA, nrow=nCV, ncol=length(Ks))
  cv.bics = matrix(NA, nrow=nCV, ncol=length(Ks))
  
  for (i.cv in c(1:nCV)) {
    cv.out = cv.error.SOUP(expr=expr, type=type, 
                           nfold=nfold, Ks=Ks)
    cv.errors[i.cv, ] = cv.out$cvm
    cv.sds[i.cv, ] = cv.out$cvsd
    cv.aics[i.cv, ] = cv.out$cvaic
    cv.bics[i.cv, ] = cv.out$cvbic
  }
  
  cvm = colMeans(cv.errors)
  cvaic = colMeans(cv.aics)
  cvbic = colMeans(cv.bics)
  
  K.aic = Ks[which.min(cvaic)]
  K.bic = Ks[which.min(cvbic)]
  K.cv = Ks[which.min(cvm)]
  K.cv.each = Ks[apply(cv.errors, 1, which.min)]
  K.majority = names(sort(table(K.cv.each), decreasing = TRUE)[1])
  
  return(list(cv.errors = cv.errors,
              cvm = cvm,
              cvsd = colMeans(cv.sds),
              cvaic = cvaic,
              cvbic = cvbic,
              K.cv = Ks[which.min(cvm)],
              K.aic = K.aic,
              K.bic = K.bic,
              K.majority=K.majority))
}

#' Compute Cross Validation Errors
#' 
#' @param expr A cell-by-gene expression matrix, either the raw counts or log-transformed expressions. 
#' @param type "log" if \code{expr} has been normalized and log-transformed (default),
#'     or "count" (default) \code{expr} contains the raw counts.
#' @param nfold Number of folds
#' @param Ks A sequence of cluster numbers
#' 
#' @export
cv.error.SOUP <- function(expr, type="log", 
                          nfold=10, Ks=c(2:10), seed=NULL) {
  n.sc = nrow(expr)
  cv.error = matrix(0, nrow=nfold, ncol=length(Ks))
  cv.aic = matrix(0, nrow=nfold, ncol=length(Ks))
  cv.bic = matrix(0, nrow=nfold, ncol=length(Ks))
  
  ## permute indices
  if (!is.null(seed)) {
    set.seed(seed)
  }
  i.permute.ind = base::sample(c(1:n.sc), size=n.sc, replace=FALSE)
  fold.size = floor(n.sc / nfold)
  
  for (fold in c(1:nfold)) {
    if (fold < nfold) {
      i.test = c( (fold.size*(fold-1)+1) : (fold.size*fold))
    } else {
      i.test = c((fold.size*(fold-1)+1) : n.sc)
    }
    i.train = c(1:n.sc)[-i.test]
    n.test = length(i.test)
    
    ## SOUP on training set
    soup.out = SOUP(expr[i.train, ], Ks=Ks, type=type)
    
    ## Predict remaining
    for (i.K in c(1:length(Ks))) {
      train.centers = soup.out$centers[[i.K]]
      test.theta = predictTheta(expr[i.test, ], 
                                t.centers=t(train.centers))
      test.pred = test.theta %*% train.centers
      cv.error[fold, i.K] = sum((expr[i.test, ] - test.pred)^2) / n.test
    }
    
    ## AIC = 2K + n*ln(2 RSS)
    cv.aic[fold, ] = 2 * Ks + n.test * log(2 * cv.error[fold, ])
    ## BIC
    cv.bic[fold, ] = log(n.test) * Ks + n.test * log(2 * cv.error[fold, ])
  }
  
  return(list(cvm = colMeans(cv.error),
              cvsd = apply(cv.error, 2, sd),
              cvaic = colMeans(cv.aic),
              cvbic = colMeans(cv.bic)
  ))
}


#' Predict the membership for new data points
#' 
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



