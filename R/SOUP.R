
#' SOUP
#'
#' A semi-soft clustering algorithm for single cells.
#' 
#' @param expr a cell-by-gene expression matrix, either the raw counts or log-transformed expressions. 
#' @param Ks number of clusters, can be a single integer or a list of integers.
#' @param type "log" if \code{expr} has been normalized and log-transformed (default),
#'     or "count" if \code{expr} contains the raw counts.
#'     It is recommended to use the log scale, which usually gives better results in practice.
#' @param i.pure (optional) the indices of the pure cells. By default is \code{NULL}, and SOUP will infer the pure list.
#' If the list is already known (for example, from previous runs), then providing it will reduce the computation time.
#' @param ext.prop (optional) the proportion of extreme neighbors for each cell, such that \code{ext.prop*n.cells} is roughly the number of pure cells \emph{per cluster}. 
#' By default, \code{ext.prop=0.1} for less than 1,000 cells, and \code{ext.prop=0.05} for larger datasets.
#' @param pure.prop (optional) the expected proportion of pure cells in the data. By default \code{pure.prop=0.5}.
#' @param nPC (optional) the number of principal components to be used, which is at least \code{K} (default). 
#' Sometimes it can be helpful to use a slightly larger \code{nPC}.
#'  
#' @return A list containing  \describe{
#'   \item{Ks}{the list of \code{K}, number of clusters, used by SOUP.}
#'   \item{memberships}{a list of cell-by-K membership matrices, one per \code{K} in \code{Ks}.}
#'   \item{centers}{a list of K-by-gene cluster center matrices, one per \code{K} in \code{Ks}.}
#'   \item{purity}{a vector containing the purity scores of all cells.}
#'   \item{i.pure}{the indices of pure cells with the highest purity scores.}
#' }
#' 
#' @export  
SOUP <- function(expr, Ks=3, 
                 type="log", 
                 i.pure=NULL, ext.prop=NULL, pure.prop=0.5,
                 nPC=NULL, nstart=50, verbose=FALSE) {
  if (! type %in% c("count", "log")) {
    stop("Data type must be eiter 'log' or 'count'.")
  }
  
  ## for raw counts, normalize by TPM
  if (type == "count") {
    expr = scaleRowSums(expr) * 10^6
  }
  
  ## find pure cells by ranking the purity score
  
  if (is.null(i.pure)) {
    if (verbose) {
      cat("Finding pure cells...\n")
    }
    A = getSimilarity(expr, type=type)
    purity.out = findPure(A=A,
                          ext.prop=ext.prop, 
                          pure.prop=pure.prop)
    i.pure = purity.out$i.pure
    purity = purity.out$purity
  } else {
    purity=NULL
  }
  
  ## SVD
  k.svd = max(nPC, max(Ks), na.rm=TRUE)
  G = RSpectra::svds(expr, k=k.svd, nu=k.svd, nv=0)$u
  
  ## SOUP
  memberships = list()
  centers = list()
  major.labels = list()
  pure.kms = list()
  for (K in Ks) {
    ## by default nPC=K, and should be at least K
    nPC = max(nPC, K, na.rm=TRUE)
    if (verbose) {
      cat("Clustering with K=", K, ", nPC=", nPC, "...\n",sep="")
    }
    
    ## pure cells partition by K-means
    km.pure = kmeans(expr[i.pure, ], centers=K, 
                     nstart=nstart, iter.max=50)
    
    ## solve for theta
    opt.out = getTheta(expr=expr,
                       i.pure=i.pure, 
                       pure.cluster=km.pure$cluster, 
                       G=G[, 1:nPC])
    
    ## major labels
    major = apply(opt.out$membership, 1, nnet::which.is.max)
    
    ## record
    pure.kms = c(pure.kms, km.pure)
    memberships = c(memberships, list(opt.out$membership))
    centers = c(centers, list(opt.out$center))
    major.labels = c(major.labels, list(major))
  }
  
  return(list(Ks=Ks,
              memberships=memberships,
              centers=centers,
              major.labels=major.labels,
              purity=purity,
              i.pure=i.pure,
              pure.kms=pure.kms))
  
}

#' Find pure cells
#' 
#' Find the list of pure cells with the highest purity scores.
#' 
#' @param A cell-by-cell similarity matrix. 
#' @param ext.prop (optional) the proportion of extreme neighbors for each cell, such that \code{ext.prop*n.cells} is roughly the number of pure cells \emph{per cluster}. 
#' By default, \code{ext.prop=0.1} for less than 1,000 cells, and \code{ext.prop=0.05} for larger datasets.
#' @param pure.prop (optional) the proportion of pure cells in the data. By default \code{pure.prop=0.5}.
#'
#' @return A list containing \describe{
#'   \item{i.pure}{indices for the top \code{pure.prop*n.cells} pure cells.}
#'   \item{purity}{purity scores of all cells.}
#' }
#' 
#' @export
findPure <- function(A, ext.prop = NULL, pure.prop = 0.5) {
  ## cell-cell similarity matrix
  n = nrow(A)
  
  ## proportion of extreme neighbors
  if (is.null(ext.prop)) {
    if (n < 1e3) {
      ext.prop = 0.1
    } else if (n < 2*1e3) {
      ext.prop = 0.05
    } else {
      ext.prop = 0.03
    }
  }
  
  ## for stable results, require at least 5 cells in the extreme neighbors
  n.ext = round(n * ext.prop)
  if (n.ext < 5) {
    warning("ext.prop is too small: at least 5 extreme neighbors are required. Use ext.prop=5/n.cells instead.\n")
    n.ext = 5
  }
  
  ## for each row i, compute its extreme value m_i
  ext.values = rep(0, n)
  for (i in c(1:n)) {
    ext.values[i] = max(abs(A[i, -i]))
  }
  
  ## purity score
  purity = rep(0, n)
  names(purity) = colnames(A)
  for (i in c(1:n)) {
    ## extreme neighbors (excluding itself)
    ext.neighbors = order(abs(A[i, ]), decreasing = TRUE)[1:n.ext]
    ext.neighbors = ext.neighbors[ext.neighbors != i]
    
    ## are extreme neighbors themselves extreme values?
    extreme.score = rep(0, length(ext.neighbors))
    for (index in c(1:length(ext.neighbors))) {
      j = ext.neighbors[index]
      extreme.score[index] = abs(A[i, j]) / ext.values[j]
    }
    
    ## purity score
    purity[i] = mean(extreme.score)
  }
  
  ## pure cells
  ## for stable results, require at least 10 pure cells
  n.pure = round(pure.prop * n)
  if (n.pure < 10) {
    warning("pure.prop is too small: at least 10 pure cells are required. Use pure.prop=10/n.cells instead.\n")
    n.pure = 10
  }
  i.pure = order(purity, decreasing = TRUE)[1:n.pure]
  
  return(list(i.pure=i.pure, purity=purity))
}

#' Get cluster centers
#' 
#' Compute the centers of clusters.
#' 
#' @param expr a cell-by-gene expression matrix, either the raw counts or log-transformed expressions. 
#' @param theta the cell-by-K membership matrix
#' 
#' @return A K-by-gene matrix representing the centers of K clusters.
#' @export
getCenters <- function(expr, theta) {
  n.gene = ncol(expr)
  K = ncol(theta)

  ## solve for X = theta %*% centers
  centers = matrix(0, nrow=K, ncol=n.gene)
  for (i in c(1:n.gene)) {
    centers[, i] = limSolve::lsei(A = theta, 
                                  B=expr[, i], 
                                  type=2)$X
  }
  colnames(centers) <- colnames(expr)

  return(centers)
}


#' Get the membership matrix and centers
#' 
#' Compute the membership matrix and cluster centers based on pure cells.
#' 
#' @param expr a cell-by-gene expression matrix, either the raw counts or log-transformed expressions. 
#' @param i.pure the indices of the pure cells..
#' @param pure.cluster a vector containing {1, ..., K}, indicating the cluster labels of pure cells.
#' @param G the cell-by-nPC matrix obtained from SVD of the similarity matrix. 
#' 
#' @return A list containing  \describe{
#'   \item{membership}{the cell-by-K membership matrix.}
#'   \item{centers}{the K-by-gene membership matrix.}
#' }
#' 
#' @export
getTheta <- function(expr, i.pure, pure.cluster, G) {
  if (length(i.pure) != length(pure.cluster)) {
    stop("Length of i.pure does not match with pure.cluster.")
  }
  
  K = length(table(pure.cluster))
  nPC = ncol(G)
  
  ## construct theta for the pure cells
  n.pure = length(i.pure)
  pure.theta = matrix(0, nrow=n.pure, ncol=K)
  for (i in c(1:n.pure)) {
    pure.theta[i, pure.cluster[i]] = 1
  }
  
  ## estimate the rotation matrix Q such that
  ## theta = G %*% Q
  Q = matrix(0, nrow=nPC, ncol=K)
  for (k in c(1:K)) {
    Q[, k] = limSolve::lsei(A = G[i.pure, ], 
                            B=pure.theta[, k], 
                            type=2)$X
  }
  theta = G %*% Q
  row.names(theta) = row.names(expr)

  ## get cluster centers
  centers = getCenters(expr, theta)
  
  ## post-process for a proper membership matrix
  membership = projMembership(theta)

  return(list(membership=membership, 
              centers=centers))
  
}

#' Clean up membership matrix
#' 
#' @param theta The estimated raw theta
#' 
#' @return The cleaned-up membership matrix.
projMembership <- function(theta) {
  membership = theta
  membership[membership < 0] = 0
  membership = scaleRowSums(membership)
  membership[is.na(membership)] = 0
  row.names(membership) = row.names(theta)
  
  return(membership)
}

#' Scale matrix by row sums
#' 
#' Scale matrix such that each row sums to 1.
#' 
#' @param x a matrix
#' @return the scaled matrix.
#' @export
#' 
scaleRowSums <- function(x) {
  return (t(scale(t(x), center=FALSE, scale=rowSums(x))))
}

#' Compute similarity matrix
#' 
#' Compute the cell-cell similarity matrix.
#' 
#' @param expr a cell-by-gene expression matrix, either the raw counts or log-transformed expressions. 
#' @param type "log" if \code{expr} has been normalized and log-transformed (default), 
#'            or "count" if \code{expr} contains the raw counts.
#' @param method the method to measure similarity between cells, one of \{"cov", "cor", "inner\},
#'
#' @return The cell-cell similarity matrix.
getSimilarity <- function(expr, type="log") {
  
  ## for raw counts, normalize and log-transform
  if (type == "count") {
    depths = rowSums(expr)
    A = expr %*% t(expr) - diag(depths)
  } else if (type == "log") {
    ## center by genes is helpful
    expr = scale(expr, center=TRUE, scale=FALSE)
    A = expr %*% t(expr)
  } else {
    stop("data type must be either 'log' or 'count'.\n")
  }

  return (A) 
}

