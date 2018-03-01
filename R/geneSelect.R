#' Select highly variable genes
#' 
#' Select highly variable genes for SOUP clustering, using SPCA and/or DESCEND.
#' 
#' @param expr a cell-by-gene expression matrix, either the raw counts or log-transformed expressions. 
#' @param type "count" (default) if \code{expr} contains the raw counts, or "log" if \code{expr} has been normalized and log-transformed.
#' @param SPCA boolean, whether to use SPCA or not.
#' @param DESCEND boolean, whether to use DESCEND or not.
#' @param n.cores the number of cores used for parallel computing of DESCEND. 
#' DESCEND can be slow so parallelization is highly recommended.
#' @param threshold the threshold for Gini index of DESCEND. Higer threshold leads to fewer selected genes.
#' @param sumabs a measurement of sparsity of genes in SPCA, between \code{1/sqrt(n.gene)} and 1.
#'  Smaller values result in sparser results, hence fewer selected genes. 
#' @param nPC the number of sparse singular vectors to use in SPCA.
#' 
#' @return A list of gene names that are selected.
#' 
#' @export
selectGenes <- function(expr, type="count",
                        SPCA=TRUE, DESCEND=TRUE,
                        n.cores=1, threshold=3,
                        sumabs=0.05, nPC=3) {
  if (! type %in% c("count", "log")) {
    stop("Data type must be eiter 'count' or 'log'.")
  }
  
  ## remove genes with too few expression
  n = nrow(expr)
  if (n > 1000) {
    min.count = 10
  } else if (n > 500) {
    min.count = 5
  } else {
    min.count = 2
  }
  
  num.express = colSums(expr > 0)
  cat("Removed", sum(num.express < min.count), 
      "genes that are expressed in less than", 
      min.count,
      "cells.\n")
  expr = expr[, num.express >= min.count]
  cat("Selecting from remaining", ncol(expr), "genes...\n")
  
  ## Gene selection
  select.genes = c()
  results = list()
  
  if (SPCA) {
    cat("SPCA selection...\n")
    spca.out = SPCAselect(expr=expr, type=type, 
                          sumabs=sumabs, nPC=nPC)
    spca.genes = spca.out$select.genes
    select.genes = c(select.genes, spca.genes)
    results = c(results, 
                list(spca.vectors=spca.out$vectors,
                     spca.genes=spca.genes))
  }
  
  if (DESCEND) {
    cat("DESCEND selection...\n")
    ## DESCEND requires raw counts, so make "pseudo" counts if input is log scale
    if (type == "log") {
      expr=2^expr - 1
      expr=round(expr)
    }
    descend.out = DESCENDselect(counts=expr, n.cores=n.cores, 
                                threshold=threshold)
    descend.genes = descend.out$select.genes
    select.genes = c(select.genes, descend.genes)
    results = c(results, 
                list(descend.scores=descend.out$scores,
                     descend.genes=descend.genes))
  }
  
  ## final list of genes
  select.genes = unique(select.genes)
  results$select.genes = select.genes
  
  return(results)
}

#' DESCEND gene selection
#' 
#' Select highly variable genes for clustering using DESCEND.
#' 
#' @param counts the cell-by-gene expression counts. Note that the data should \emph{not} be normalized or log-transformed.
#' @param n.cores the number of cores used for parallel computing. DESCEND can be slow so parallelization is highly recommended.
#' @param threshold the threshold for Gini index. Higer threshold leads to fewer selected genes.
#' 
#' @return A list containing \describe{
#'   \item{select.genes}{the names of selected genes, ordered by decreasing scores.}
#'   \item{scores}{a score matrix containing the Gini scores of all genes.}
#' }
#' 
#' @export
DESCENDselect <- function(counts, n.cores=1, threshold=3) {
  ## DESCEND requires gene-by-cell raw count as input
  result <- runDescend(t(as.matrix(counts)), 
                       n.cores = n.cores,
                       do.LRT.test = FALSE,
                       family = "Poisson")
  
  hvg <- findHVG(result, criteria="Gini")
  
  return(list(select.genes = hvg$HVG.genes,
              scores=hvg$score.mat))
}

#' SPCA gene selection
#' 
#' Select highly variable genes for clustering using Sparse Principal Component Analysis (SPCA).
#' 
#' @param expr a cell-by-gene expression matrix, either the raw counts or log-transformed expressions. 
#' @param type "count" (default) if \code{expr} contains the raw counts, or "log" if \code{expr} has been normalized and log-transformed.
#' @param sumabs a measurement of sparsity for SPCA, between \code{1/sqrt(n.gene)} and 1.
#'  Smaller values result in sparser results, hence fewer selected genes. 
#' @param nPC the number of sparse singular vectors to look into.
#' 
#' @return A list containing \describe{
#'   \item{select.genes}{the names of selected genes, ordered by decreasing importance.}
#'   \item{vectors}{a gene-by-nPC matrix of the sparse eigen vectors.}
#' }
#' 
#' @references 
#' Witten, DM and Tibshirani, R and T Hastie (2009) 
#' A penalized matrix decomposition, with applications to 
#' sparse principal components and canonical correlation analysis.
#' \emph{Biostatistics}.
#' 
#' @export
SPCAselect <- function(expr, type="log", sumabs=0.05, nPC=3) {
  ## for raw counts, normalize and log-transform
  if (type == "count") {
    expr = scaleRowSums(expr) * 10^6
    expr = log2(expr + 1)
  } else if (type != "log") {
    stop("Data type must be eiter 'count' or 'log'.")
  }
  
  ## SPCA
  n.sc = nrow(expr)
  n.gene = ncol(expr)
  ## sumabsv must be between 1 and sqrt(n.gene)
  sumabsv = max(1, sumabs*sqrt(n.gene))
  spca.out <- PMA::PMD(scale(expr), ## scale by genes is helpful 
                      type="standard", 
                      sumabs=NULL, 
                      sumabsu=sqrt(n.sc),  ## no sparsity for cells
                      sumabsv=sumabsv, ## sparsity for genes
                      K=nPC,
                      center=TRUE,
                      trace=FALSE)
  
  ## selected genes
  gene.scores = rowSums((spca.out$v)^2)
  select.genes = colnames(expr)[order(gene.scores, decreasing = TRUE)[1:sum(gene.scores > 0)]]
  vectors=spca.out$v
  row.names(vectors) = colnames(expr)
  
  return(list(select.genes = select.genes,
              vectors=vectors))
}