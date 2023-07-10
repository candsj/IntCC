
<!-- README.md is generated from README.Rmd. Please edit that file -->

# weightedCC

<!-- badges: start -->

[![R-CMD-check](https://github.com/candsj/weightedCC/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/candsj/weightedCC/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->


## Installation

You can install the development version of weightedCC from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("candsj/weightedCC")
```

# Introduction
Analysis of omic datasets has become more and more important recently. It could help to define cancer subtypes and reveal novel discoveries. How to integrate omic datasets raised much-growing attention. The Cancer Genome Atlas Program began in 2006 and generated over 2.5 petabytes of genomic, epigenomic, transcriptomic, and proteomic data. This greatly improved the development of omic datasets research. Various integrative clustering methods and tools have been proposed. We proposed two-layer weighted consensus clustering method to deal with some limitations existed.

# Usage of weighted consensus clustering

The R package **wcc** (**w**eighted **c**onsensus **c**lustering) introduces weighted consensus clustering functions inspired by consensus clustering methods.

The main function in this package is `WeightedConsensusCluster`. Users input a list of datasets, clusteing methods, the number of clutsers for separate datasets and global, criterions and so on. Then consensus clustering will be applied to each dataset to derive consensus matrix. Depends on one layer or two-layer, consensus matrix will be weighted and combined for each dataset and final clustering. Finally user could choose pam or hierarchical clustering method to derive the clustering results. Some useful functions are also included in this package: 

* `criterion`. This function takes a list of consensus matrix, a list of corresponding clustering label results and a vector of corresponding number of clusters. It could use Silhouette, widestGap, dunns and dunn2s index as criterion to choose the best number of clusters.

* `weightcal`. This function takes a list of consensus matrix and clustering label results. Then weights are calculated based on ratio of in-cluster proportion to out-of-cluster proportion using the cluster estimated by the algorithm itself. It could help to estimate weights for different methods of each methods or weights for consensus matrix of each dataset.

The other function needed for this is `consensusCluster`. This function can be found in the R package `coca` and is used to perform consensus clustering on one dataset and obtain a co-clustering matrix  (Monti et al. 2003). We made some changes to this function to allow users to define more distance metric and sample percentage for features.


# How wcc works 

First, we generate four different type datasets with the same clustering structure (3 clusters of equal size) and different levels of noise.
