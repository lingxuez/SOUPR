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

## 2. SOUP
## log-transformation recommended
log.expr = log2(scaleRowSums(counts) * 1e6 + 1)

## selected genes
select.genes = zeisel$select.genes
## 2818 selected genes
select.expr = log.expr[, colnames(log.expr) %in% select.genes]
dim(select.expr)
## 3005 2818

## SOUP with a single K
## takes ~1.5 minutes
system.time({
  soup.out = SOUP(select.expr, Ks=7, type="log")
})

## Treat SOUP as hard clustering
membership = soup.out$memberships[[1]]
soup.label = apply(membership, 1, nnet::which.is.max)
cat("ARI =", 
    mclust::adjustedRandIndex(soup.label, cell.info$cell.type), 
    "\n")
## ARI = 0.8944478 

## visualize the contingency table
g <- plotContTable(true_label=cell.info$cell.type,
              est_label=soup.label,
              xlab="Zeisel")
g <- g + theme(axis.text.x=element_text(size=12, angle=90))
g

############################
## Example on Camp data
############################

## load data: 220 fetal brain single cells
counts = camp$counts
cell.info = camp$cell.info
dim(counts)
## 220 12694
table(cell.info$cell.type)
## 7 cell types

## 1. Select genes using DESCEND and SPCA.
## This step can be slow, and please consider using parallelization.
## The selected results has been provided, so you can skip this step.

## Input: a cell-by-gene count matrix.
# select.out = selectGenes(counts, type="count", n.cores=10)
# select.genes = select.out$select.genes

## 2. SOUP
## log-transformation recommended
log.expr = log2(scaleRowSums(counts) * 1e6 + 1)

## selected genes
select.genes = camp$select.genes
## 430 selected genes
select.expr = log.expr[, colnames(log.expr) %in% select.genes]
dim(select.expr)
## 220 430

## SOUP with a sequence of Ks
system.time({
  soup.out = SOUP(select.expr, Ks=c(2:7), type="log")
})

## visualize the SOUP assignments of different Ks
## compared tp the reference labels
heatmapKseq(memberships=soup.out$memberships, 
                 Ks=soup.out$Ks,
                 cell.type=cell.info$cell.type)

## get developmental trajectory
