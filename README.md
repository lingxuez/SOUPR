# SOUP

R package for SOUP.

SOUP, for Semi-sOft clUstering with Pure cells, is a single-cell clustering method
that can classify both pure and transitional cells. 
SOUP performs well at the hard-clustering problem for pure cell types 
and excels at identifying transitional cells with soft memberships.
To find pure cells, SOUP uses the special block structure the $K$ cell types form in a similarity matrix, 
devised by pairwise comparison of the gene expression profiles of individual cells. 
Once pure cells are identified, they provide the key information from which the membership matrix can be computed. 
SOUP is applicable to general clustering problems as well, as long as the unrestrictive modeling assumptions hold. 

## Citation
As part of the gene selection procedure, we utilized the R package [DESCEND](https://github.com/jingshuw/descend).
To avoid compatibility issues, we directly included the three source R scripts, 
`deconvSingle.R`, `descend.R`, and `g_model.R`, 
forked from its [Github Repo](https://github.com/jingshuw/descend/tree/master/R). 



## Installation
This package can be installed through `devtools` in R:
```{r}
install.packages("devtools") ## if not installed
library("devtools")
devtools::install_github("lingxuez/SOUP")
```

## Examples
Please follow the vignette for an example of using this R package on a fetal brain dataset, 
where we re-produce the results in our paper.
