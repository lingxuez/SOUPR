#' Plot lineages or curves
#' 
#' @export
plotLineageCurves <- function(x,
                              centers=NULL,
                          dims=c(1:2),
                          add=FALSE,
                          labels = NULL, 
                          lwd=2, 
                          col=1, 
                          lab.cex=2, 
                          cex=2, 
                          pos=2, 
                          offset=0.5,
                          type = "lineages") {
  
  curves <- FALSE
  lineages <- FALSE
  if(type %in% c('lineages','both')){
    lineages <- TRUE
  }
  if(type %in% c('curves','both')){
    curves <- TRUE
  }
  if(lineages & (length(x@lineages)==0)){
    stop('No lineages detected.')
  }
  if(curves & (length(x@curves)==0)){
    stop('No curves detected.')
  }
  
  ## info
  if(lineages){
    X <- reducedDim(x)
    clusterLabels <- clusterLabels(x)
    connectivity <- slingAdjacency(x)
    clusters <- rownames(connectivity)
    nclus <- nrow(connectivity)
    if (is.null(centers)) {
      centers <- t(vapply(clusters,function(clID){
        w <- clusterLabels[,clID]
        return(apply(X, 2, weighted.mean, w = w))
      }, rep(0,ncol(X))))
    }
    
    rownames(centers) <- clusters
    X <- X[rowSums(clusterLabels) > 0, , drop = FALSE]
    clusterLabels <- clusterLabels[rowSums(clusterLabels) > 0, ,drop = FALSE]
    linC <- slingParams(x)
  }
  
  ## basic plot
  if(!add){
    xs <- NULL
    ys <- NULL
    if(lineages){
      xs <- c(xs, centers[,dims[1]])
      ys <- c(ys, centers[,dims[2]])
    }
    if(curves){
      xs <- c(xs, as.numeric(vapply(slingCurves(x), 
                                    function(cr){ cr$s[,dims[1]] }, 
                                    rep(0,nrow(reducedDim(x))))))
      ys <- c(ys, as.numeric(vapply(slingCurves(x), 
                                    function(cr){ cr$s[,dims[2]] }, 
                                    rep(0,nrow(reducedDim(x))))))
    }
    plot(x = NULL, y = NULL, #asp = asp,
         xlim = range(xs), ylim = range(ys),
         xlab = colnames(reducedDim(x))[dims[1]],
         ylab = colnames(reducedDim(x))[dims[2]])
  }
  
  ## lineage
  if(lineages){
    for(i in 1:(nclus-1)){
      for(j in (i+1):nclus){
        if(connectivity[i,j]==1){
          lines(centers[c(i,j), dims], lwd = lwd, col = col)
        }
      }
    }
    
    points(centers[, dims], cex = cex, pch = 16)
    text(centers[, dims], 
         labels=labels, pos=pos, cex=lab.cex, offset = offset)
  }
  
  if(curves){
    for(ii in seq_along(slingCurves(x))){
      cr <- slingCurves(x)[[ii]]
      lines(cr$s[cr$ord, dims], lwd = lwd, col = col)
    }
  }
}
