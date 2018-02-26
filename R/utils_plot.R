
#' Visualize SOUP hard assignments.
#' 
#' @param memberships A list of membership matrices of different K
#' @param Ks A list of different Ks 
#' @param cell.type The golden standard cell types for reference.
#' 
#' @return A ggplot object.
heatmapKseq <- function(memberships, Ks, cell.type) {
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
    assign.Kseq[, i] <- apply(memberships[[i]][, i.permute], 
                        1, nnet::which.is.max)
  }
  timeline = rowSums(assign.Kseq %*% diag(Ks))
  
  ## reference type
  assign.Kseq = data.frame(assign.Kseq)
  assign.Kseq$Reference = as.numeric(as.factor(cell.type))
  colnames(assign.Kseq) <- c(paste0("K=", Ks), "Reference")
  
  ## Heatmap
  dat3 <- assign.Kseq[order(cell.type, timeline), ]
  dat3$Cell = factor(c(1:n.sc))
  dat3 <- melt(dat3, id.var="Cell")
  dat3$value <- as.factor(dat3$value)
  
  K = length(table(dat3$value))
  mycols = RColorBrewer::brewer.pal(K+1, "Spectral")[-c((K+1)/2+1)]
  
  g <- ggplot2::ggplot(dat3, aes(Cell, variable)) +
    geom_tile(aes(fill = value)) +
    scale_fill_manual(values=mycols) +
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
  
  return(g)
}


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
  
  return(order(mapped.labels))
}