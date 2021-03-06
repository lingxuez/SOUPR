% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{zeisel}
\alias{zeisel}
\title{The mouse brain single-cell dataset}
\format{A list with three items:
\describe{
  \item{cell.info}{a dataframe with 3,005 rows and two columns, cell.id and cell.type.}
  \item{counts}{a matrix of RNA-seq counts, with 3,005 rows and 16,450 columns.}
  \item{select.genes}{names of the 2818 genes being selected by DESCEND and SPCA.}
}}
\source{
Zeisel \emph{et al.} (2015) 
Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq.
\emph{Science}. (\url{http://linnarssonlab.org/cortex/})
}
\usage{
zeisel
}
\description{
An RNA-seq dataset published in
Zeisel \emph{et al.} (2015),
containing 3,005 mouse brain single cells. 
Cells are grouped into 7 major cell types, 
and we retain the 16,450 genes that are expressed in at least 10 cells.
}
\examples{
cell.info = zeisel$cell.info
counts = zeisel$counts
select.genes = zeisel$select.genes
}
\keyword{datasets}
