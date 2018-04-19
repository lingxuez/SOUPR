
#' Get lineages from cluster centers
#' 
#' Given cluster centers, construct a linearge using minimum-spanning tree.
#' A starting or ending cluster can be specified.
#' This code is adapted from the `getLineages` method in R package `slingshot`.
#' 
#' @param centers K-by-p matrix of cluster centers
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
    newtree = get_connections(unused[1], forest)
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
#' @param lineages A list of cluster lineages
#' @param membership cell-by-K membership matrix
#' 
#' @return The pseudotime along each lineage, where cells are assigned to
#' different lineages according to their major clusters
#' 
#' @export
getLineageTime <- function(lineages, membership) {
  ## major cluster
  cell.major.clust = apply(membership, 1, nnet::which.is.max)
  cell.lineage.time = matrix(NA, nrow=nrow(membership), ncol=length(lineages))
  colnames(cell.lineage.time) = names(lineages)
  
  if (is.null(colnames(membership))) {
    colnames(membership) <- c(1:ncol(membership))
  }
  
  ## time along each lineage
  for (l in c(1:length(lineages))) {
    i.clust = match(lineages[[l]], colnames(membership))
    i.cell = which(cell.major.clust %in% i.clust)
    
    ## only use the largest 2 proportions
    # scaled.membership = t(apply(membership[i.cell, i.clust], 1, function(x){
    #   i.top2 = order(x, decreasing=TRUE)[1:2]
    #   x[-i.top2] = 0
    #   x = x / sum(x)
    #   return (x)
    # }))
    scaled.membership = scaleRowSums(membership[i.cell, i.clust])
    l.time = rowSums(scaled.membership %*% diag(1:length(i.clust)))
    cell.lineage.time[i.cell, l] = l.time
  }
  
  return(cell.lineage.time)
}

#' Helper function
#' 
#' Extract all descendants starting from `clus`
#' 
#' @param clust the cluster to start with
#' @param forest the adjacency matrix
#' @param parent parent that will be excluded from the result
#' 
#' @return The list of descendants
get_connections <- function(clus, forest, parent = NULL){
  children.idx <- which(forest[,clus] == 1)
  children <- rownames(forest)[children.idx]
  if(is.null(parent)){
    out <- clus
    for(child in children){
      out <- c(out, Recall(child, forest, clus))
    }
  }else{
    children <- children[children != parent]
    out <- clus
    for(child in children){
      out <- c(out, Recall(child, forest, clus))
    }
  }
  return(out)
}
