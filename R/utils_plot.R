
#' Visualize SOUP hard assignments.
#' 
#' @param memberships A list of membership matrices of different K
#' @param Ks A list of different Ks 
#' @param cell.type The golden standard cell types for reference.
#' 
#' @return A ggplot object.
heatmapKseq <- function(memberships, Ks, cell.type) {
  
  assign.out = getMajorMatrix(memberships, Ks, cell.type,
                                type="hard")
  ## order
  i.order = order(cell.type, assign.out$timeline)
  dat3 <- assign.out$assign.Kseq[i.order, ]

  ## Heatmap
  dat3$Cell = factor(c(1:nrow(dat3)))
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


#' Get the SOUP majority matrix
#' 
#' @export
#' 
getMajorMatrix <- function(memberships, Ks, cell.type,
                            type="hard") {
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
    i.na <- which(rowSums(memberships[[i]]) == 0)
    ## re-order the labels to match
    i.permute = getPermute(soup.label, cell.type)
    if (type == "hard") {
      assign.Kseq[, i] <- apply(memberships[[i]][, i.permute], 
                                1, nnet::which.is.max)
      assign.Kseq[i.na, i] <- NA
    } else {
      assign.Kseq[, i] <- apply(memberships[[i]][, i.permute], 
                                1, max)
    }
   
  }
  ## timeline
  last_membership = memberships[[n.K]][, i.permute]
  timeline = rowSums(last_membership  %*% diag(c(1:ncol(last_membership))))
  
  ## reference type
  assign.Kseq = data.frame(assign.Kseq)
  assign.Kseq$Reference = as.numeric(as.factor(cell.type))
  colnames(assign.Kseq) <- c(paste0("K=", Ks), "Reference")
  
  return(list(assign.Kseq=assign.Kseq,
              timeline=timeline,
              last_membership=last_membership,
              i.permute=i.permute))
}


#' Visualize SOUP major proportions
#' 
#' @param memberships A list of membership matrices of different K
#' @param Ks A list of different Ks 
#' @param cell.type The golden standard cell types for reference.
#' 
#' @return A ggplot object.
heatmapKseqProp <- function(memberships, Ks, cell.type) {
  assign.out = getMajorMatrix(memberships, Ks, cell.type,
                              type="soft")
  ## order
  i.order = order(cell.type, assign.out$timeline)
  dat3 <- assign.out$assign.Kseq[i.order, c(1:length(Ks))]
  
  ## Heatmap
  dat3$Cell = factor(c(1:nrow(dat3)))
  dat3 <- melt(dat3, id.var="Cell")
  
  mycol = RColorBrewer::brewer.pal(7, "Spectral")[7]
  g <- ggplot2::ggplot(dat3, aes(Cell, variable)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient(low="white", high=mycol) +
    labs(y="", x="Single Cells", fill="Major\nproportion") +
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
  
  return(as.numeric(row.names(cont.table)[order(mapped.labels)]))
}


plotContTable <- function(true_label, label_est, short.names=NULL, xlab="Reference") {
  if (is.null(short.names)) {
    short.names = levels(factor(true_label))
  }
  cont.table <- table(true_label, label_est)
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