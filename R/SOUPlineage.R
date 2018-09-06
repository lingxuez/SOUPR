
#' Get lineages from cluster centers
#' 
#' Given cluster centers, construct a linearge using minimum spanning tree (MST) with Euclidean distances.
#' It is often suggested to use a low-dimensional space for identifying MST.
#' A starting or ending cluster can be specified if prior knowledge is available.
#' This code is adapted from the `getLineages` method in the R package `slingshot`.
#' 
#' @param centers K-by-p' matrix of cluster centers, where p' is usually for the low-dimensional space.
#' @param end.clust (optional) the end point
#' @param start.clust (optional) the starting point
#' 
#' @references Street K, Risso D, Fletcher RB, et al. (2017). 
#' "Slingshot: Cell lineage and pseudotime inference for single-cell transcriptomics."
#'  \emph{bioRxiv. 10.1101/128843}
#'  
#' @export

getClusterLineages <- function(centers,
                               end.clust=NULL,
                               start.clust=NULL,
                               dist.fun=NULL) {
  nclust = nrow(centers)
  if (is.null(row.names(centers))) {
    row.names(centers) = c(1:nclust)
  }
  
  ####################
  ## Forest
  ####################
  
  ## pairwise cluster distance; use Euclidean
  if (is.null(dist.fun)) {
    D = as.matrix(dist(centers, method="euclidean"))
  }
  
  ## add an artificial cluster OMEGA; let its distance larger than any other
  omega = max(D) + 1
  D = rbind(D, rep(omega, nclust))
  D = cbind(D, c(rep(omega, nclust), 0))
  row.names(D) = c(row.names(centers), "O")
  colnames(D) = c(row.names(centers), "O")
  
  
  ## draw MST on cluster centers + OMEGA
  ## exclude endpoint cluster(s) if specified
  if (!is.null(end.clust)) {
    i.end = which(row.names(D) %in% end.clust)
    mstree = ape::mst(D[-i.end, -i.end, drop=FALSE]) 
  } else {
    mstree = ape::mst(D)
  }
  
  ## add in endpoint cluster(s)
  if (!is.null(end.clust)) {
    forest = matrix(0, nrow=nrow(D), ncol=ncol(D))
    row.names(forest) = row.names(D)
    colnames(forest) = colnames(D)
    forest[-i.end, -i.end] = mstree
    ## append end point to closest one
    for (cl in end.clust) {
      i.cl = which(row.names(D) == cl)
      dists = D[-i.end, i.cl]
      closest = names(dists)[which.min(dists)]
      i.closest = which(row.names(D) == closest)
      forest[i.cl, i.closest] = 1
      forest[i.closest, i.cl] = 1
    }
  } else {
    forest = mstree
  }
  
  ## remove OMEGA
  forest = forest[1:nclust, 1:nclust, drop=FALSE]
  
  ###################
  ## Lineage
  ###################
  
  ## identify trees  
  unused = row.names(forest)
  trees = list()
  ntree = 0
  while(length(unused) > 0) {
    ntree = ntree + 1
    newtree = getDescendants(unused[1], forest)
    trees[[ntree]] = newtree
    unused = setdiff(unused, newtree)
  }
  trees <- trees[order(sapply(trees,length),decreasing = TRUE)]
  
  ## identify lineages (paths through trees)
  lineages = list()
  for (tree in trees) {
    ## isolated node
    if (length(tree) == 1) {
      lineages <- c(lineages, list(tree))
      next
    }
    
    ## tree
    i.tree = which(row.names(forest) %in% tree)
    tree.graph = forest[i.tree, i.tree, drop=FALSE]
    degree = rowSums(tree.graph)
    leaves = rownames(tree.graph)[degree == 1]
    g = igraph::graph.adjacency(tree.graph, mode="undirected")
    
    ## pick the root
    ## if you have starting cluster(s) in this tree then they are roots
    if (!is.null(start.clust) && sum(start.clust %in% tree) > 0) {
      starts = intersect(start.clust, tree)
    } else {
      ## else, pick the start (root) with higest average length (~parsimony)
      avg.length = sapply(leaves, function(l){
        ends = leaves[leaves != l]
        paths = igraph::shortest_paths(g, from=l, to=ends, mode="out", output="vpath")$vpath
        mean(sapply(paths, length))
      })
      starts = names(avg.length)[nnet::which.is.max(avg.length)]
    }
    
    ## draw a path from each start to each leaf
    ends = setdiff(leaves, starts)
    for (st in starts) {
      paths = igraph::shortest_paths(g, from=st, to=ends, mode="out", output="vpath")$vpath
      for (p in paths) {
        lineages = c(lineages, list(names(p)))
      }
    }
    
  }
  
  ## sort lineages by number of clusters included
  lineages = lineages[order(sapply(lineages, length), decreasing=TRUE)]
  names(lineages) = paste0("Lineage", 1:length(lineages))
  
  return(list(lineages=lineages,
              D=D[1:nclust, 1:nclust],
              forest=forest,
              trees=trees,
              start.clust=start.clust,
              end.clust=end.clust))

}

#' Decide cell pseudotime
#' 
#' Compute the pseudotime of each cell along each lineage.
#' 
#' @param lineages A list of vectors, one per lineage, 
#'              each containing the ordered cluster labels along the lineage
#' @param membership the n-by-K soft membership matrix
#' 
#' @return An n-by-L pseudotime matrix for L lineages. 
#'     Each column contains the pseudotime of cells along one lineage, and NA if a cell does not belong to it. 
#'     Cells are assigned to different lineages, potentially with overlaps, according to their major clusters.
#' 
#' @export
getLineageTime <- function(lineages, membership) {
  ## the pseudotime matrix, one column per lineage
  cell.lineage.time = matrix(NA, nrow=nrow(membership), ncol=length(lineages))
  colnames(cell.lineage.time) = names(lineages)
  
  if (is.null(colnames(membership))) {
    colnames(membership) <- c(1:ncol(membership))
  }
  
  ## compute pseudotime along each lineage
  cell.major.clust = apply(membership, 1, nnet::which.is.max)
  for (l in c(1:length(lineages))) {
    ## the cells belonging to this lineage, decided based on their major clusters
    i.clust = match(lineages[[l]], colnames(membership))
    i.cell = which(cell.major.clust %in% i.clust)
    
    ## psuedotime is computed from soft memberships
    scaled.membership = scaleRowSums(membership[i.cell, i.clust])
    l.time = rowSums(scaled.membership %*% diag(1:length(i.clust)))
    cell.lineage.time[i.cell, l] = l.time
  }
  
  return(cell.lineage.time)
}

#' Get descendants in a tree
#' 
#' Helper function to extract all descendants of `clus` in a tree, one node per cluster.
#' This code is adapted from the R package `slingshot`.
#' 
#' @param clust the cluster (node) to start with
#' @param forest the adjacency matrix between clusters (nodes)
#' @param parent the parent cluster (node) that will be excluded from the result
#' 
#' @references Street K, Risso D, Fletcher RB, et al. (2017). 
#' "Slingshot: Cell lineage and pseudotime inference for single-cell transcriptomics."
#'  \emph{bioRxiv. 10.1101/128843}
#'   
#' @return A vector containing the list of descendants
#' 
#' @export
getDescendants <- function(clus, forest, parent = NULL){
  ## children (excluding its parent)
  children.idx = which(forest[, clus] == 1)
  children = rownames(forest)[children.idx]
  if (!is.null(parent)) {
    children = children[children != parent]
  }
  
  ## recursively find all descendants
  out = clus
  for(child in children){
    out <- c(out, Recall(clus=child, forest=forest, parent=clus))
  }

  return(out)
}


#' Get lineage curves
#' 
#' This code is adapted from the `getCurves` method in the R package `slingshot`, modified to allow soft memberships.
#' 
#' @param expr.ld The n-by-p' data matrix, usually projected to a low-dimensional space with a small p'
#' @param lineages A list of vectors, one per lineage, 
#'              each containing the ordered cluster labels along the lineage
#' @param centers K-by-p' matrix of cluster centers in the low-dimensional space.
#'             The row names are the cluster labels, which should to be consistent with the names in `lineages`.
#' @param shrink (optional) logical or numeric between 0 and 1, determines whether and how 
#'   much to shrink branching lineages toward their average prior to the split.
#' @param smoother (optional) choice of scatter plot smoother. Same as 
#'   \code{\link{principal.curve}}, but \code{"lowess"} option is replaced with 
#'   \code{"loess"} for additional flexibility.
#' @param maxit (optional) maximum iteration
#' @param extend character, how to handle root and leaf clusters of lineages 
#'   when constructing the initial, piece-wise linear curve. Accepted values are
#'   \code{'y'} (default) and \code{'n'}. See 'Details' for more.
#' @param reweight logical, whether to allow cells shared between lineages to be
#'   reweighted during curve-fitting. If \code{TRUE}, cells shared between 
#'   lineages will be weighted by: distance to nearest curve / distance to
#'   curve.
#' @param drop.multi logical, whether to drop shared cells from lineages which 
#'   do not fit them well. If \code{TRUE}, shared cells with a distance to one 
#'   lineage above the 90th percentile and another below the 50th will be
#'   dropped from the further lineage.
#' @param shrink.method character denoting how to determine the appropriate 
#'   amount of shrinkage for a branching lineage. Accepted values are the same
#'   as for \code{kernel} in \code{\link{density}} (default is \code{"cosine"}),
#'   as well as \code{"tricube"} and \code{"density"}. See 'Details' for more.
#'   
#' @details When there is only a single lineage, the curve-fitting algorithm is 
#'   nearly identical to that of \code{\link{principal.curve}}. When there are 
#'   multiple lineages and \code{shrink == TRUE}, an additional step is added to
#'   the iterative procedure, forcing curves to be similar in the neighborhood 
#'   of shared points (ie., before they branch).
#'   
#' @details The \code{extend} argument determines how to construct the
#'   piece-wise linear curve used to initiate the recursive algorithm. The
#'   initial curve is always based on the lines between cluster centers and if
#'   \code{extend = 'n'}, this curve will terminate at the center of the
#'   endpoint clusters. Setting \code{extend = 'y'} will allow the first and
#'   last segments to extend beyond the cluster center to the orthogonal
#'   projection of the furthest point.
#'   These options typically have little to no impact on the final curve, but can
#'   occasionally help with stability issues.
#'   
#' @details When \code{shink == TRUE}, we compute a shrinkage curve,
#'   \eqn{w_l(t)}, for each lineage, a non-increasing function of pseudotime
#'   that determines how much that lineage should be shrunk toward a shared
#'   average curve. We set \eqn{w_l(0) = 1}, so that the curves will perfectly
#'   overlap the average curve at pseudotime \code{0}. The weighting curve
#'   decreases from \code{1} to \code{0} over the non-outlying pseudotime values
#'   of shared cells (where outliers are defined by the \code{1.5*IQR} rule).
#'   The exact shape of the curve in this region is controlled by
#'   \code{shrink.method}, and can follow the shape of any standard kernel
#'   function's cumulative density curve (or more precisely, survival curve,
#'   since we require a decreasing function). Different choices of
#'   \code{shrink.method} seem to have little impact on the final curves, in 
#'   most cases.
#'   
#'   
#' @references Hastie, T., and Stuetzle, W. (1989). "Principal Curves."
#'   \emph{Journal of the American Statistical Association}, 84:502–516.
#' @references Street K, Risso D, Fletcher RB, et al. (2017). 
#' "Slingshot: Cell lineage and pseudotime inference for single-cell transcriptomics."
#'  \emph{bioRxiv. 10.1101/128843}
#'  
#' @importFrom princurve project_to_curve
#' @export
getLineageCurves <- function(expr.ld, 
                             lineages, 
                             centers, 
                             membership,
                             shrink=TRUE, 
                             smoother = 'smooth.spline',
                             extend="y",
                             drop.multi = TRUE,
                             reweight=TRUE,
                             shrink.method = 'cosine',
                             stretch = 2,
                             thresh=0.001, maxit = 15) {
  
  ## Setup
  shrink = as.numeric(shrink)
  L = length(lineages)
  n = nrow(expr.ld); p = ncol(expr.ld)
  nclust = nrow(centers)
  if (is.null(row.names(centers))) {
    row.names(centers) = c(1:nclust)
  }
  if (is.null(colnames(membership))) {
    colnames(membership) = c(1:nclust)
  }
  if (!all.equal(row.names(centers), colnames(membership))) {
    stop("Cluster labels do not match between center matrix and membership matrix.\n")
  }
  clusters = row.names(centers)
  
  ## cells' weights in each lineage
  major.clust = apply(membership, 1, nnet::which.is.max)
  W = sapply(c(1:L), function(l){
    as.numeric(major.clust %in% lineages[[l]])
  })
  row.names(W) = row.names(membership)
  colnames(W) = names(lineages)
  
  W.orig = W
  ## distance from each cell to each lineage
  D = W; D[,] = NA
  
  ## smoother function
  smootherFcn <- switch(smoother, loess = function(lambda, xj, 
                                                   w = NULL, ...){
    loess(xj ~ lambda, weights = w, ...)$fitted
  }, smooth.spline = function(lambda, xj, w = NULL, ..., df = 5, 
                              tol = 1e-4){
    fit <- tryCatch({
      smooth.spline(lambda, xj, w = w, ..., df = df, 
                    tol = tol, keep.data = FALSE)
    }, error = function(e){
      smooth.spline(lambda, xj, w = w, ..., df = df, 
                    tol = tol, keep.data = FALSE, spar = 1)
    })
    predict(fit, x = lambda)$y
  })
  
  # determine curve hierarchy
  ## which clusters are in which lineage(s)
  C = as.matrix(sapply(lineages, function(lin) {
    sapply(clusters, function(clID) {
      as.numeric(clID %in% lin)
    })
  }))
  rownames(C) = clusters
  segmnts = unique(C[rowSums(C)>1,,drop = FALSE])
  segmnts = segmnts[order(rowSums(segmnts),decreasing = FALSE), ,
                     drop = FALSE]
  avg.order = list()
  for(i in seq_len(nrow(segmnts))){
    idx = (segmnts[i,] == 1)
    avg.order[[i]] = colnames(segmnts)[idx]
    new.col = rowMeans(segmnts[,idx, drop = FALSE])
    segmnts = cbind(segmnts[, !idx, drop = FALSE],new.col)
    colnames(segmnts)[ncol(segmnts)] = paste('average',i,sep='')
  }
  
  # initial curves are piecewise linear paths through the tree
  pcurves = list()
  for(l in seq_len(L)){
    ## which cells to use to for each lineage
    ## they will extend the initial path if extend == 'y'
    idx = which(W[,l] > 0.8)
    
    ## centers along the lineage
    line.initial <- centers[clusters %in% lineages[[l]], , 
                            drop = FALSE]
    line.initial <- line.initial[match(lineages[[l]],
                                       rownames(line.initial)),  ,
                                 drop = FALSE]
    K <- nrow(line.initial)
    
    
    # special case: single-cluster lineage
    if(K == 1){
      pca <- prcomp(expr.ld[idx, ,drop = FALSE])
      ctr <- line.initial
      line.initial <- rbind(ctr - 10*pca$sdev[1] * 
                              pca$rotation[,1], 
                            ctr, 
                            ctr + 10*pca$sdev[1] *
                              pca$rotation[,1])
      curve <- project_to_curve(expr.ld[idx, ,drop = FALSE], 
                                s = line.initial,
                                stretch = 9999)
      
      # do this twice because all points should have projections
      # on all lineages, but only those points on the lineage
      # should extend it
      pcurve <- project_to_curve(expr.ld, 
                                 s = curve$s[curve$ord, ,drop = FALSE],
                                 stretch=0)
      pcurve$dist <- abs(pcurve$dist) 
      # ^ force non-negative distances
      pcurve$lambda <- pcurve$lambda - min(pcurve$lambda, 
                                           na.rm=TRUE)
      # ^ force pseudotime to start at 0
      pcurve$w <- W[,l]
      pcurves[[l]] <- pcurve
      D[,l] <- abs(pcurve$dist)
      next
    }
    
    ## project in-branch cells to the initial piece-wise linear curve
    ## potentially extend the line in the end point
    if(extend == 'y'){
      curve <- project_to_curve(expr.ld[idx, ,drop = FALSE], 
                                s = line.initial,
                                stretch = 9999)
      curve$dist <- abs(curve$dist)
    }
    if(extend == 'n'){
      curve <- project_to_curve(expr.ld[idx, ,drop = FALSE], 
                                s = line.initial,
                                stretch = 0)
      curve$dist <- abs(curve$dist)
    }
    
    ## project *all* cells to the curve, without extending it
    pcurve <- project_to_curve(expr.ld, 
                               s = curve$s[curve$ord, ,drop = FALSE],
                               stretch=0)
    
    # force non-negative distances
    pcurve$dist <- abs(pcurve$dist)
    # force pseudotime to start at 0
    pcurve$lambda <- pcurve$lambda - min(pcurve$lambda, na.rm=TRUE) 
    pcurve$w <- W[,l]
    pcurves[[l]] <- pcurve
    D[,l] <- abs(pcurve$dist)
  }
  
  ## track distances between curves and data points to determine convergence
  dist.new <- sum(abs(D[W>0.8]), na.rm=TRUE)
  
  it <- 0
  hasConverged <- FALSE
  while (!hasConverged && it < maxit){
    it <- it + 1
    dist.old <- dist.new
    
    ## reweight cells by distance
    if(reweight){
      W[,] <- t(sapply(seq_len(nrow(W)),function(i){
        ds <- D[i,]
        out <- min(ds)/ds
        out[is.nan(out)] <- 1 # handle 0/0
        return(out)
      }))
      W[W > 1] <- 1
      W[W < 0] <- 0
      W[W.orig==0] <- 0
    }
    
    ## drop shared cells
    if(drop.multi){
      Z <- D; Z[,] <- NA
      for(l in seq_len(L)){
        idx <- W[,l] > 0
        Z[idx,l] <- rank(D[idx,l]) / sum(idx)
      }
      W[,] <- t(sapply(seq_len(nrow(W)),function(i){
        zs <- Z[i,]
        out <- W[i,]
        ## if a cell with a distance to one lineage abote 90th quantile,
        ## and another below 50th, it'll be dropped from the further lineage
        if(max(zs,na.rm=TRUE) > .9 && 
           min(zs,na.rm=TRUE) <= .5){
          out[!is.na(zs) & zs > .9] <- 0
        }
        return(out)
      }))
    }
    
    # predict each dimension as a function of lambda (pseudotime)
    for(l in seq_len(L)){
      pcurve <- pcurves[[l]]
      ## previous proj
      s <- pcurve$s 
      ord <- order(pcurve$lambda)
      ## smooth predict
      for(jj in seq_len(p)){
        s[, jj] <- smootherFcn(pcurve$lambda, expr.ld[,jj], 
                               w = pcurve$w)[ord]
      }
      ## new projection
      new.pcurve <- project_to_curve(expr.ld, s = s, stretch = stretch)
      new.pcurve$dist <- abs(new.pcurve$dist)
      new.pcurve$lambda <- new.pcurve$lambda - 
        min(new.pcurve$lambda, na.rm = TRUE)
      new.pcurve$w <- W[,l]
      pcurves[[l]] <- new.pcurve
    }
    D[,] <- sapply(pcurves, function(p){ p$dist })
    
    # shrink together lineages near shared clusters
    if(shrink > 0){
      if(max(rowSums(C)) > 1){
        
        segmnts <- unique(C[rowSums(C)>1,,drop=FALSE])
        segmnts <- segmnts[order(rowSums(segmnts),
                                 decreasing = FALSE),
                           , drop = FALSE]
        seg.mix <- segmnts
        avg.lines <- list()
        pct.shrink <- list()
        
        # determine average curves and amount of shrinkage
        for(i in seq_along(avg.order)){
          ns <- avg.order[[i]]
          to.avg <- lapply(ns,function(n){
            if(grepl('Lineage',n)){
              l.ind <- as.numeric(gsub('Lineage','',n))
              return(pcurves[[l.ind]])
            }
            if(grepl('average',n)){
              a.ind <- as.numeric(gsub('average','',n))
              return(avg.lines[[a.ind]])
            }
          })
          avg <- avg_curves(to.avg, expr.ld, stretch = stretch)
          avg.lines[[i]] <- avg
          common.ind <- rowMeans(sapply(to.avg,
                                        function(crv){
                                          crv$w > 0})) == 1
          pct.shrink[[i]] <- lapply(to.avg,function(crv){
            percent_shrinkage(crv, common.ind, 
                               method = shrink.method)
          })
          # check for degenerate case (if one curve won't be
          # shrunk, then the other curve shouldn't be,
          # either)
          new.avg.order <- avg.order
          all.zero <- sapply(pct.shrink[[i]], function(pij){
            return(all(pij == 0))
          })
          if(any(all.zero)){
            if(allow.breaks){
              new.avg.order[[i]] <- NULL
              message('Curves for', ns[1], 'and', ns[2],
                      'appear to be going in opposite 
                      directions. No longer forcing them
                      to share an initial point. To 
                      manually override this, set 
                      allow.breaks = FALSE.')
            }
            pct.shrink[[i]] <- lapply(pct.shrink[[i]], 
                                      function(pij){
                                        pij[] <- 0
                                        return(pij)
                                      })
            }
          }
        # do the shrinking in reverse order
        for(j in rev(seq_along(avg.lines))){
          ns <- avg.order[[j]]
          avg <- avg.lines[[j]]
          to.shrink <- lapply(ns,function(n){
            if(grepl('Lineage',n)){
              l.ind <- as.numeric(gsub('Lineage','',n))
              return(pcurves[[l.ind]])
            }
            if(grepl('average',n)){
              a.ind <- as.numeric(gsub('average','',n))
              return(avg.lines[[a.ind]])
            }
          })
          shrunk <- lapply(seq_along(ns),function(jj){
            crv <- to.shrink[[jj]]
            return(shrink_to_avg(crv, avg, 
                                  pct.shrink[[j]][[jj]], 
                                  expr.ld, stretch = stretch))
          })
          for(jj in seq_along(ns)){
            n <- ns[jj]
            if(grepl('Lineage',n)){
              l.ind <- as.numeric(gsub('Lineage','',n))
              pcurves[[l.ind]] <- shrunk[[jj]]
            }
            if(grepl('average',n)){
              a.ind <- as.numeric(gsub('average','',n))
              avg.lines[[a.ind]] <- shrunk[[jj]]
            }
          }
        }
        avg.order <- new.avg.order
      }
    }
    D[,] <- sapply(pcurves, function(p){ p$dist })
    
    ## update distance
    dist.new <- sum(D[W>0.8], na.rm=TRUE)
    hasConverged <- (abs((dist.old - 
                            dist.new)/dist.old) <= thresh)
  }
  
  ## final processing
  if(reweight){
    W[,] <- t(sapply(seq_len(nrow(W)),function(i){
      ds <- D[i,]
      out <- min(ds)/ds
      out[is.nan(out)] <- 1 # handle 0/0
      return(out)
    }))
    W[W > 1] <- 1
    W[W < 0] <- 0
    W[W.orig==0] <- 0
  }
  if(drop.multi){
    Z <- D; Z[,] <- NA
    for(l in seq_len(L)){
      idx <- W[,l] > 0
      Z[idx,l] <- rank(D[idx,l]) / sum(idx)
    }
    W[,] <- t(sapply(seq_len(nrow(W)),function(i){
      zs <- Z[i,]
      out <- W[i,]
      if(max(zs,na.rm=TRUE) > .9 && min(zs,na.rm=TRUE) <= .5){
        out[!is.na(zs) & zs > .9] <- 0
      }
      return(out)
    }))
  }
  
  for(l in seq_len(L)){
    class(pcurves[[l]]) <- 'principal.curve'
  }
  names(pcurves) <- paste('curve',1:length(pcurves),sep='')
  
  for(l in seq_len(L)){
    pcurve$pseudotime <- pcurve$lambda
    pcurve$w <- W[,l]
    pcurve$pseudotime[pcurve$w==0] <- NA
  }
  
  return(pcurves)
}


avg_curves <- function(pcurves, X, stretch = 2){
  p <- ncol(pcurves[[1]]$s)
  lambdas.all <- lapply(pcurves, function(pcv){pcv$lambda})
  lambdas.all <- unique(unlist(lambdas.all))
  max.shared.lambda <- min(sapply(pcurves, function(pcv){max(pcv$lambda)}))
  lambdas.all <- sort(lambdas.all[lambdas.all <= max.shared.lambda])
  pcurves.dense <- lapply(pcurves,function(pcv){
    sapply(seq_len(p),function(jj){
      interpolated <- approx(pcv$lambda, pcv$s[,jj], xout = lambdas.all)$y
      return(interpolated)
    })
  })
  avg <- sapply(seq_len(p),function(jj){
    dim.all <- sapply(1:length(pcurves.dense),function(i){
      pcurves.dense[[i]][,jj]
    })
    return(rowMeans(dim.all))
  })
  avg.curve <- princurve::project_to_curve(X, avg, stretch=stretch)
  avg.curve$w <- rowMeans(sapply(pcurves, function(p){ p$w }))
  return(avg.curve)
}

percent_shrinkage <- function(crv, share.idx, method = 'cosine'){
  pst <- crv$lambda
  if(method %in% eval(formals(density.default)$kernel)){
    dens <- density(0, bw=1, kernel = method)
    surv <- list(x = dens$x, y = (sum(dens$y) - cumsum(dens$y))/sum(dens$y))
    box.vals <- boxplot(pst[share.idx], plot = FALSE)$stats
    surv$x <- scaleAB(surv$x, a = box.vals[1], b = box.vals[5])
    if(box.vals[1]==box.vals[5]){
      pct.l <- rep(0, length(pst))
    }else{
      pct.l <- approx(surv$x, surv$y, pst, rule = 2)$y
    }
  }
  if(method == 'tricube'){
    tc <- function(x){ ifelse(abs(x) <= 1, (70/81)*((1-abs(x)^3)^3), 0) }
    dens <- list(x = seq(-3,3,length.out = 512))
    dens$y <- tc(dens$x)
    surv <- list(x = dens$x, y = (sum(dens$y) - cumsum(dens$y))/sum(dens$y))
    box.vals <- boxplot(pst[share.idx], plot = FALSE)$stats
    surv$x <- scaleAB(surv$x, a = box.vals[1], b = box.vals[5])
    if(box.vals[1]==box.vals[5]){
      pct.l <- rep(0, length(pst))
    }else{
      pct.l <- approx(surv$x, surv$y, pst, rule = 2)$y
    }
  }
  if(method == 'density'){
    bw1 <- bw.SJ(pst)
    bw2 <- bw.SJ(pst[share.idx])
    bw <- (bw1 + bw2) / 2
    d2 <- density(pst[share.idx], bw = bw, 
                  weights = crv$w[share.idx]/sum(crv$w[share.idx]))
    d1 <- density(pst, bw = bw, weights = crv$w/sum(crv$w))
    scale <- sum(crv$w[share.idx]) / sum(crv$w)
    pct.l <- (approx(d2$x,d2$y,xout = pst, yleft = 0, 
                     yright = 0)$y * scale) / 
      approx(d1$x,d1$y,xout = pst, yleft = 0, yright = 0)$y
    pct.l[is.na(pct.l)] <- 0
    pct.l <- cumMin(pct.l, pst)
  }
  return(pct.l)
}

shrink_to_avg <- function(pcurve, avg.curve, pct, X, stretch = 2){
  n <- nrow(pcurve$s)
  p <- ncol(pcurve$s)
  lam <- pcurve$lambda
  s <- vapply(seq_len(p),function(jj){
    orig.jj <- pcurve$s[,jj]
    avg.jj <- approx(x = avg.curve$lambda, y = avg.curve$s[,jj], xout = lam,
                     rule = 2)$y
    return(avg.jj * pct + orig.jj * (1-pct))
  }, rep(0,n))
  w <- pcurve$w
  pcurve <- princurve::project_to_curve(X, 
                                        as.matrix(s[pcurve$ord, ,drop = FALSE]), 
                                        stretch = stretch)
  pcurve$w <- w
  return(pcurve)
}


cumMin <- function(x,time){
  sapply(seq_along(x),function(i){ min(x[time <= time[i]]) })
}

scaleAB <- function(x,a=0,b=1){
  ((x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE)))*(b-a)+a
}
