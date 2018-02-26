library(mclust)
library(devtools)
library(ggplot2)

## set working dir to this folder
## setwd("SOUP/")
devtools::load_all()

############################
## Example on Zeisel data
############################

## load data: 3,005 mouse brain single cells
counts = zeisel$counts
cell.info = zeisel$cell.info
dim(counts)
## 3005 16450
table(cell.info$cell.type)
## 7 major cell types

## 1. Select genes using DESCEND and SPCA.
## This step can be slow, and please consider using parallelization.
## The selected results has been provided, so you can skip this step.

## Input: a cell-by-gene count matrix.
# select.out = selectGenes(counts, type="count", n.cores=10)
# select.genes = select.out$select.genes

## 2. SOUP on 2,819 selected genes
select.genes = zeisel$select.genes
sel.counts = counts[, colnames(counts) %in% select.genes]
system.time({
  ## Input: a cell-by-gene count matrix.
  ## You can also use a sequence of Ks, e.g., Ks=c(5:7)
  soup.out = SOUP(sel.counts, Ks=7, type="count")
})

membership = soup.out$memberships[[1]]
soup.label = apply(membership, 1, nnet::which.is.max)
cat("ARI =", 
    mclust::adjustedRandIndex(soup.label, cell.info$cell.type), 
    "\n")
## ARI = 0.8804754

############################
## Example on Camp data
############################

## load data: 220 fetal brain single cells
counts = camp$counts
cell.info = camp$cell.info
dim(counts)
## 220 12858
table(cell.info$cell.type)
## 7 cell types

## 1. Select genes using DESCEND and SPCA.
## This step can be slow, and please consider using parallelization.
## The selected results has been provided, so you can skip this step.

## Input: a cell-by-gene count matrix.
# select.out = selectGenes(counts, type="count", n.cores=10)
# select.genes = select.out$select.genes

## 2. SOUP on 425 selected genes
select.genes = camp$select.genes
sel.counts = counts[, colnames(counts) %in% select.genes]
system.time({
  ## Input: a cell-by-gene count matrix.
  soup.out = SOUP(sel.counts, Ks=c(3:7), type="count")
})

## visualize the SOUP assignments versus the reference labels
g <- heatmapKseq(memberships=soup.out$memberships, 
                 Ks=soup.out$Ks,
                 cell.type=cell.info$cell.type)
g
