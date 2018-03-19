
#' Estimate Developmental Timeline
#' 
#' Estimate the developmental timeline for cells from SOUP soft membership.
#' 
#' @param membership An n.cell-by-K soft membership matrix
#' @param centers A K-by-n.gene cluster center matrix
#' @param k.start Which cluster is the starting point? 
#'     One, and only one of \code{k.start} and \code{k.end} must be specified.
#' @param k.end Which cluster is the ending point?
#'     One, and only one of \code{k.start} and \code{k.end} must be specified.
#'     
#' @return A vector of timepoints for all cells.
#' @export
getTimeline <- function(membership, centers, 
                        k.start=NULL, k.end=NULL) {
  
  K = nrow(centers)
  if (ncol(membership) != K) {
    stop("Error: Membership must have the same ncol as centers.\n")
  }
  
  ## re-order clusters
  cor.center = cor(t(centers))
  if (!is.null(k.start)) {
    k.order = getClusterOrder(cor.center, k.start)
  } else if (!is.null(k.end)) {
    ## start from end, then reverse
    k.order = getClusterOrder(cor.center, k.end)
    k.order = k.order[c(K:1)]
  }
  
  ## timeline
  timeline = rowSums(membership[, k.order] %*% diag(c(1:K)))
  names(timeline) <- row.names(membership)
  
  return(timeline)
}

#' Order Clusters
#' 
#' Re-order cluster labels, starting from \code{k.start}, 
#' and consecutively append the next one with the highest correlation in cluster centers.
#' 
#' @param cor.center The K-by-K correlation matrix among cluster centers
#' @param k.start The cluster to start with
#' 
#' @return A list with length K that re-orders the clusters. 
#' @export
getClusterOrder <- function(cor.center, ## K-by-K correlation matrix of centers
                            k.start=1 ## which cluster to start with
) {
  K = ncol(cor.center)
  ## order clusters
  colnames(cor.center) = c(1:K)
  k.order = c(k.start)
  k.cur = k.start
  while (length(k.order) < K) {
    k.next = names(which.max(cor.center[k.cur, -k.order]))
    k.next = as.numeric(k.next)
    k.order = c(k.order, k.next)
    k.cur = k.next
  }
  return (k.order)
}
