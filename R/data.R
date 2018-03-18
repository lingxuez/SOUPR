#' The mouse brain single-cell dataset
#'
#' An RNA-seq dataset published in
#' Zeisel \emph{et al.} (2015),
#' containing 3,005 mouse brain single cells. 
#' Cells are grouped into 7 major cell types, 
#' and we retain the 16,450 genes that are expressed in at least 10 cells.
#'
#' @format A list with three items:
#' \describe{
#'   \item{cell.info}{a dataframe with 3,005 rows and two columns, cell.id and cell.type.}
#'   \item{counts}{a matrix of RNA-seq counts, with 3,005 rows and 16,450 columns.}
#'   \item{select.genes}{names of the 2818 genes being selected by DESCEND and SPCA.}
#' }
#' 
#' @source 
#' Zeisel \emph{et al.} (2015) 
#' Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq.
#' \emph{Science}. (\url{http://linnarssonlab.org/cortex/})
#' 
#' @examples 
#' cell.info = zeisel$cell.info
#' counts = zeisel$counts
#' select.genes = zeisel$select.genes
"zeisel"


#' The fetal brain single-cell dataset.
#'
#' An RNA-seq dataset published in
#' Camp \emph{et al.} (2015),
#' containing 220 fetal brain single cells. 
#' Cells are grouped into 7 cell types, 
#' and we retain the 12,694 genes that are expressed in at least 2 cells.
#'
#' @format A list with three items:
#' \describe{
#'   \item{cell.info}{a dataframe with 220 rows and two columns, cell.id and cell.type.}
#'   \item{counts}{a matrix of RNA-seq counts, with 220 rows and 12,694 columns.}
#'   \item{select.genes}{names of the 430 genes being selected by DESCEND and SPCA.}
#' }
#' 
#' @source 
#' Camp \emph{et al.} (2015) 
#' Human cerebral organoids recapitulate gene expression programs of fetal neocortex development.
#' \emph{PNAS}.
#' 
#' @examples 
#' cell.info = camp$cell.info
#' counts = camp$counts
#' select.genes = camp$select.genes
"camp"
