
#' SOUP
#'
#' A semi-soft clustering algorithm for single cells.
#' 
#' @param expr a cell-by-gene expression matrix, either the raw counts or log-transformed expressions. 
#' @param Ks number of clusters, can be a single integer or a list of integers.
#' @param type "count" (default) if \code{expr} contains the raw counts, or "log" if \code{expr} has been normalized and log-transformed.
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
#' @examples 
#' select.genes = zeisel$select.genes
#' counts = zeisel$counts[, colnames(counts) %in% select.genes]
#' soup.out = SOUP(counts, Ks=7, type="count")
#' 
#' @export  

SOUP <- function(expr, Ks=3, 
                 type="count", 
                 i.pure=NULL, ext.prop=NULL, pure.prop=0.5,
                 nPC=NULL, verbose=FALSE) {
  if (! type %in% c("count", "log")) {
    stop("Data type must be eiter 'count' or 'log'.")
  }
  
  ## for raw counts, normalize
  if (type == "count") {
    expr = scaleRowSums(expr) * 10^6
  }
  
  ## cell-cell similarity matrix
  if (verbose) {
    cat("Computing similarity matrix...\n")
  }
  A = getSimilarity(expr, type=type)
  
  ## find pure cells by ranking the purity score
  if (is.null(i.pure)) {
    if (verbose) {
      cat("Finding pure cells...\n")
    }
    purity.out = findPure(A=A,
                          ext.prop=ext.prop, 
                          pure.prop=pure.prop)
    i.pure = purity.out$i.pure
    purity = purity.out$purity
  } else {
    purity=NULL
  }
  
 
  ## SVD: only do it once
  if (type == "count") {
    G = RSpectra::eigs_sym(A, k=max(nPC, max(Ks), na.rm=TRUE))$vectors
  } else {
    G = RSpectra::svds(expr, k=max(nPC, max(Ks), na.rm=TRUE))$u
  }
  
  ## SOUP
  memberships = list()
  centers = list()
  for (K in Ks) {
    if (verbose) {
      cat("Clustering with K =", K, "...\n")
    }
    
    ## by default nPC=K, and should be at least K
    if (is.null(nPC)) {
      nPC = K
    }
    nPC = max(nPC, K, na.rm=TRUE)
    
    ## pure cells partition by K-means
    km.pure = kmeans(expr[i.pure, ], centers=K, 
                     nstart=20, iter.max=50)
    
    ## solve for theta
    opt.out = getTheta(expr=expr,
                       i.pure=i.pure, 
                       pure.cluster=km.pure$cluster, 
                       G=G[, 1:nPC])
    
    memberships = c(memberships, list(opt.out$membership))
    centers = c(centers, list(opt.out$center))
  }
  
  return(list(Ks=Ks,
              memberships=memberships,
              centers=centers,
              purity=purity,
              i.pure=i.pure))
  
}

#' Find pure cells
#' 
#' Find the list of pure cells with the highest purity scores.
#' 
#' @param A the cell-cell similarity matrix. 
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
  n = nrow(A)
  
  ## proportion of extreme neighbors
  if (is.null(ext.prop)) {
    if (n < 1e3) {
      ext.prop = 0.1
    } else {
      ext.prop = 0.05
    }
  }
  ## for stable results, require at least 5 cells in the extreme neighbors
  n.ext = round(n * ext.prop)
  n.ext = max(n.ext, 5)
  
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
  n.pure = round(pure.prop*nrow(A))
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

  ## get cluster centers
  centers = getCenters(expr, theta)
  
  ## post-process for a proper membership matrix
  membership = theta
  membership[membership < 0] = 0
  membership = scaleRowSums(membership)
  membership[is.na(membership)] = 0
  row.names(membership) = row.names(expr)
  
  return(list(membership=membership, 
              centers=centers))
  
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
#' @param type "count" (default) if \code{expr} contains the raw counts, or "log" if \code{expr} has been normalized and log-transformed.
#' @param method the method to measure similarity between cells, one of \{"cov", "cor", "inner\},
#'
#' @return The cell-cell similarity matrix.
getSimilarity <- function(expr, type="count",
                          method="inner") {
  
  ## for raw counts, normalize and log-transform
  if (type == "count") {
    depths = rowSums(expr)
    A = expr %*% t(expr) - diag(depths)
  } else if (type == "log") {
    ## center by genes is helpful
    expr = scale(expr, center=TRUE, scale=FALSE)
    A = expr %*% t(expr)
  } else {
    stop("data type must be 'count' or 'log'.\n")
  }
  
  # 
  # ## similarity matrix: covariance, correlation, or inner product?
  # ## For now, inner product works the best.
  # if (method == "cov") {
  #   A = cov(t(expr))
  # } else if (method == "cor") {
  #   A = cor(t(expr))
  # } else if (method == "inner") {
  #   ## center by genes is helpful
  #   expr = scale(expr, center=TRUE, scale=FALSE)
  #   A = expr %*% t(expr)
  # } else {
  #   stop("Similarity method must be 'cov', 'cor', or 'inner'.")
  # }
  return (A) 
}

