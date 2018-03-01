
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
  theta[theta < 0] = 0
  theta = scaleRowSums(theta)
  theta[is.na(theta)] = 0
  
  return(theta)
  
}

#' Cross Validation Errors
#' 
#' @export
cv.error.SOUP <- function(expr, type="count", 
                    nfold=10, Ks=c(2:10)) {
  n.sc = nrow(expr)
  cv.error = matrix(0, nrow=nfold, ncol=length(Ks))
  
  ## permute
  i.permute.ind = base::sample(c(1:n.sc), size=n.sc, replace=FALSE)
  fold.size = floor(n.sc / nfold)
  
  fold = 1
  # for (fold in c(1:nfold)) {
    if (fold < nfold) {
      i.test = c( (fold.size*(fold-1)+1) : (fold.size*fold))
    } else {
      i.test = c((fold.size*(fold-1)+1) : n.sc)
    }
    i.train = c(1:n.sc)[-i.test]
    
    ## SOUP on training set
    soup.out = SOUP(expr[i.train, ], Ks=Ks, type=type)
    
    ## Predict remaining
    for (i.K in c(1:length(Ks))) {
      train.centers = soup.out$centers[[i.K]]
      test.theta = predictTheta(expr[i.test, ], 
                                centers=t(train.centers))
      test.pred = test.theta %*% train.centers
      cv.error[fold, i.K] = Matrix::norm(expr[i.test, ] - test.pred, "F")
    }
  # }
  
  return(list(cvm = colMeans(cv.error),
              cvsd = apply(cv.error, 2, sd)))
}