
#' Visualize SOUP hard assignments.
#' 
#' @param memberships A list of membership matrices of different K
#' @param Ks A list of different Ks 
#' @param cell.type The golden standard cell types for reference.
#' @param ref.lab Name to be shown for the reference cell types; default is "Reference"
#' 
#' @return A ggplot object.
#' 
#' @export
heatmapKseq <- function(memberships, Ks, cell.type, ref.lab="Reference") {
  
  assign.Kseq = getMajorMatrix(memberships, Ks, cell.type,
                              ref.lab=ref.lab)
  
  ## order cells by reference type and the last SOUP assignments
  i.order = order(cell.type, assign.Kseq[, paste0("K=", max(Ks))])
  dat3 <- assign.Kseq[i.order, ]

  ## Heatmap
  dat3$Cell = factor(c(1:nrow(dat3)))
  dat3 <- melt(dat3, id.var="Cell")
  dat3$value <- as.factor(dat3$value)
  g <- ggplot2::ggplot(dat3, aes(Cell, variable)) +
    geom_tile(aes(fill = value)) +
    labs(y="", x="Single Cells", fill="cluster") +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=15),
          axis.ticks=element_blank(),
          axis.title.x=element_text(size=15),
          axis.title.y=element_blank(),
          legend.position="right", 
          legend.direction="vertical",
          legend.title = element_text(size=15)
    )
  
  ## color scheme
  K = length(table(dat3$value))
  if (K <= 10) {
    mycols = RColorBrewer::brewer.pal(K+1, "Spectral")[-c((K+1)/2+1)]
    g <- g + scale_fill_manual(values=mycols)
  }
  
  return(g)
}


#' Get the SOUP majority matrix
#' 
#' @export
#' 
getMajorMatrix <- function(memberships, Ks, cell.type, ref.lab="Reference") {
  if (length(memberships) != length(Ks)) {
    stop("memberships and Ks must be two lists with same lengths.\n")
  }
  
  n.sc = nrow(memberships[[1]])
  n.K = length(Ks)
  
  ## Get a [n.sc, n.K] matrix, recording the SOUP hard assignments of each K
  assign.Kseq = matrix(nrow=n.sc, ncol=n.K)
  for (i in c(1:n.K)) {
    K <- Ks[i]
    soup.label <- apply(memberships[[i]], 1, nnet::which.is.max)
    ## re-order the labels to match
    i.permute = getPermute(soup.label, cell.type)
    assign.Kseq[, i] <- apply(memberships[[i]][, i.permute, drop=FALSE], 
                                1, nnet::which.is.max)
  }

  ## reference type
  assign.Kseq = data.frame(assign.Kseq)
  assign.Kseq$Reference = as.numeric(as.factor(cell.type))
  colnames(assign.Kseq) <- c(paste0("K=", Ks), ref.lab)
  
  return(assign.Kseq)
}


#' Permute labels
#' 
#' Permute the estimated labels to be consistent with the reference labels.
#' 
#' @param est.label a vector of estimated clusters, containing {1, ..., K}
#' @param true.label a vector of reference/true labels
#' 
#' @return A vector of length K, i.permute, 
#' such that the i-th reference cluster maps to i.permute[i] in estimated clusters. 
getPermute = function(est.label, true.label) {
  cont.table = table(est.label, true.label)
  K = nrow(cont.table)
  mapped.labels = rep(0, K)
  
  while (sum(!is.na(cont.table)) > 0) {
    ## find the next largest overlap
    i.max <- which.max(cont.table)
    max.row <- row(cont.table)[i.max]
    max.col <- col(cont.table)[i.max]
    mapped.labels[max.row] <- max.col
    
    ## remove
    cont.table[max.row, ] <- NA
    cont.table[, max.col] <- NA
  }
  
  ## remaining: randomly assign
  clust.remaining <- which(mapped.labels == 0)
  if (length(clust.remaining) > 0) {
    mapped.labels[mapped.labels==0] <- clust.remaining
  }
  
  i.permute = as.numeric(row.names(cont.table)[order(mapped.labels)])
  return(i.permute)
}

#' Visualize contingency table
#' 
#' @param est_label Estimated cluster assignments, a vector of characters or factors
#' @param true_label The true cell types, a vector of characters or factors
#' @param short.names (optional) If the true cell type names are too long, 
#'         you can supply abbreviated names to save space in visualizations
#' @param xlab (optional) The x-axis label; default is "Reference"
#' 
#' @return A ggplot object.
#' 
#' @export
plotContTable <- function(est_label, true_label, short.names=NULL, xlab="Reference") {
  if ("factor" %in% class(true_label)) {
    true_label = droplevels(true_label)
  }
  if ("factor" %in% class(est_label)) {
    est_label = droplevels(est_label)
  }
  if (is.null(short.names)) {
    short.names = levels(factor(true_label))
  }
  cont.table <- table(true_label, est_label)
  K <- ncol(cont.table)
  sub.clusters <- paste0("cluster ", colnames(cont.table))
  
  cont.table <- apply(as.matrix(cont.table), 2, as.integer)
  cont.table <- data.frame(cont.table)
  cont.table$Reference = factor(short.names, levels=short.names)
  colnames(cont.table) <- c(sub.clusters, "Reference")
  
  dat3 <- melt(cont.table, id.var="Reference")
  grid.labels = as.character(dat3$value)
  grid.labels[grid.labels == "0"] = ""
  
  
  g <- ggplot(dat3, aes(Reference, variable)) +
    geom_tile(aes(fill = value)) +
    geom_text(aes(label= grid.labels), size=4.5) +
    scale_fill_gradient(low="white", high="purple") +
    labs(y="SOUP", x=xlab) +
    theme(panel.background = element_blank(),
          axis.line=element_blank(),
          axis.text.x=element_text(size=13),
          axis.text.y=element_text(size=13),
          axis.ticks=element_blank(),
          axis.title.x=element_text(size=18),
          axis.title.y=element_text(size=18),
          legend.position="none"
    )
  
  return(g)
}