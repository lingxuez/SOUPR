# SOUP

R package for SOUP.

SOUP, for Semi-sOft clUstering with Pure cells, is a clustering method for single cell data.
SOUP reveals the clustering structure for both pure cells, which belong to one single cluster, 
as well as transitional cells with soft memberships.
SOUP is applicable to general clustering problems as well, as long as the unrestrictive modeling assumptions hold. 
A preprint of the manuscript is available [here](https://www.biorxiv.org/content/early/2018/03/19/285056.full.pdf).

## Citation
Pease cite our paper if it helps your research:
```
@article{Zhu2018,
    title={Semi-soft Clustering of Single Cell Data},
    author={Zhu, Lingxue and Lei, Jing and Devlin, Bernie and Roeder, Kathryn},
    year={2018},
    journal={bioRxiv},
    doi={10.1101/285056}
}
```

In this package, we utilized the R package [DESCEND](https://github.com/jingshuw/descend).
To avoid compatibility issues, we directly included the three source R scripts, 
`deconvSingle.R`, `descend.R`, and `g_model.R`, 
forked from the [Github Repo](https://github.com/jingshuw/descend/tree/master/R). 



## Installation
This package can be installed through `devtools` in R:
```{r}
library("devtools")
devtools::install_github("lingxuez/SOUP")
```

## Examples
Please follow the [vignette](https://github.com/lingxuez/SOUP/blob/master/vignettes/SOUP-vignette.pdf) for an example of using this R package on a fetal brain dataset, which re-produces the results in our paper.
