
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

```r
library("weightedCC")
```


## Introduction
Analysis of omic datasets has become more and more important recently. It could help to define cancer subtypes and reveal novel discoveries. How to integrate omic datasets raised much-growing attention. The Cancer Genome Atlas Program began in 2006 and generated over 2.5 petabytes of genomic, epigenomic, transcriptomic, and proteomic data. This greatly improved the development of omic datasets research. Various integrative clustering methods and tools have been proposed. We proposed two-layer weighted consensus clustering method to deal with some limitations existed.

## Usage of weighted consensus clustering

The R package **wcc** (**w**eighted **c**onsensus **c**lustering) introduces weighted consensus clustering functions inspired by consensus clustering methods.

The main function in this package is `WeightedConsensusCluster`. Users input a list of datasets, clusteing methods, the number of clutsers for separate datasets and global, criterions and so on. Then consensus clustering will be applied to each dataset to derive consensus matrix. Depends on one layer or two-layer, consensus matrix will be weighted and combined for each dataset and final clustering. Finally user could choose pam or hierarchical clustering method to derive the clustering results. Some useful functions are also included in this package: 

* `criterion`. This function takes a list of consensus matrix, a list of corresponding clustering label results and a vector of corresponding number of clusters. It could use Silhouette, widestGap, dunns and dunn2s index as criterion to choose the best number of clusters.

* `weightcal`. This function takes a list of consensus matrix and clustering label results. Then weights are calculated based on ratio of in-cluster proportion to out-of-cluster proportion using the cluster estimated by the algorithm itself. It could help to estimate weights for different methods of each methods or weights for consensus matrix of each dataset.

The other function needed for this is `consensusCluster`. This function can be found in the R package `coca` and is used to perform consensus clustering on one dataset and obtain a co-clustering matrix  (Monti et al. 2003). We made some changes to this function to allow users to define more distance metric and sample percentage for features.


## How wcc works 

First, we generate four different type datasets with the same clustering structure (3 clusters of equal size) and different levels of noise.

```r
  load(system.file("extdata", "exampleData.RData", package = "weightedCC"))
  normData=exampleData[[1]]
  binomData=exampleData[[2]]
  poisData=exampleData[[3]]
  multData=exampleData[[4]]
  lopoisData=log(poisData+1)
```

Now we can use the `consensusclustering` function to compute a consensus matrix for each dataset.

```r
  normRes=consensuscluster(normData,K=3,B=1000,pItem = 0.8,pFeature = 0.8 ,clMethod ="kmeans",finalclmethod="pam")
  binomRes=consensuscluster(binomData,K=3,B=1000,pItem = 0.8,pFeature = 0.8 ,clMethod ="kmeans",finalclmethod="pam")
  poisRes=consensuscluster(lopoisData,K=3,B=1000,pItem = 0.8,pFeature = 0.8 ,clMethod ="kmeans",finalclmethod="pam")
  multRes=consensuscluster(multData,K=3,B=1000,pItem = 0.8,pFeature = 0.8 ,clMethod ="kmeans",finalclmethod="pam")

  normRes1=consensuscluster(normData,K=3,B=1000,pItem = 0.8,pFeature = 0.8 ,clMethod ="hclust",finalclmethod="pam")
  binomRes1=consensuscluster(binomData,K=3,B=1000,pItem = 0.8,pFeature = 0.8 ,clMethod ="hclust",dist = "binary",finalclmethod="pam")
  poisRes1=consensuscluster(poisData,K=3,B=1000,pItem = 0.8,pFeature = 0.8 ,clMethod ="hclust",dist = "jaccard",finalclmethod="pam")
  multRes1=consensuscluster(multData,K=3,B=1000,pItem = 0.8,pFeature = 0.8 ,clMethod ="hclust",dist = "hamming",finalclmethod="pam")
```

Then we could perform one layer weighted consensus clustering.

```r
  res.all <- c(list(normRes),list(binomRes),list(poisRes),list(multRes))
  weights <-weightcal(res.all)
  wcm=weights[1]*res.all[[1]]$consensusMatrix+weights[2]*res.all[[2]]$consensusMatrix+
    weights[3]*res.all[[3]]$consensusMatrix+weights[4]*res.all[[4]]$consensusMatrix
  distances <- stats::as.dist(1 - wcm)
  weight_clusterLabels <- pam(distances,3,diss = TRUE,metric = "euclidean" )$clustering
  ari.weight<- adjustedRandIndex(true.class, weight_clusterLabels)
```

Or we could perform two-layer weighted consensus clustering.

```r
  #norm
  res.norm <- c(list(normRes),list(normRes1))
  weight.norm<- weightcal(res.norm)
  wcm=weight.norm[1]*res.norm[[1]]$consensusMatrix+weight.norm[2]*res.norm[[2]]$consensusMatrix
  distances <- stats::as.dist(1 - wcm)
  weight_clusterLabels <-pam(distances,3,diss = TRUE,metric = "euclidean" )$clustering
  weightnormRes=list(consensusMatrix=wcm,class=weight_clusterLabels)
  
  #binom
  res.binom <- c(list(binomRes),list(binomRes1))
  weight.binom <- weightcal(res.binom)
  wcm=weight.binom[1]*res.binom[[1]]$consensusMatrix+weight.binom[2]*res.binom[[2]]$consensusMatrix
  distances <- stats::as.dist(1 - wcm)
  weight_clusterLabels <- pam(distances,3,diss = TRUE,metric = "euclidean" )$clustering
  weightbinomRes=list(consensusMatrix=wcm,class=weight_clusterLabels)
  
  #poisson
  res.pois <- c(list(poisRes),list(poisRes1))
  weight.pois <- weightcal(res.pois)
  wcm=weight.pois[1]*res.pois[[1]]$consensusMatrix+weight.pois[2]*res.pois[[2]]$consensusMatrix
  distances <- stats::as.dist(1 - wcm)
  weight_clusterLabels <-pam(distances,3,diss = TRUE,metric = "euclidean" )$clustering
  weightpoisRes=list(consensusMatrix=wcm,class=weight_clusterLabels)
  
  #mult
  res.mult <- c(list(multRes),list(multRes1))
  weight.mult <- weightcal(res.mult)
  wcm=weight.mult[1]*res.mult[[1]]$consensusMatrix+weight.mult[2]*res.mult[[2]]$consensusMatrix
  distances <- stats::as.dist(1 - wcm)
  weight_clusterLabels <- pam(distances,3,diss = TRUE,metric = "euclidean" )$clustering
  weightmultRes=list(consensusMatrix=wcm,class=weight_clusterLabels)
  
  #2 layer wcm
  res.wcm <- c(list(weightnormRes),list(weightbinomRes),list(weightpoisRes),list(weightmultRes))
  weight.wcm <- weightcal(res.wcm)
  wcm=weight.wcm[1]*res.wcm[[1]]$consensusMatrix+weight.wcm[2]*res.wcm[[2]]$consensusMatrix+
    weight.wcm[3]*res.wcm[[3]]$consensusMatrix+weight.wcm[4]*res.wcm[[4]]$consensusMatrix
  distances <- stats::as.dist(1 - wcm)
  weight_clusterLabels <- pam(distances,3,diss = TRUE,metric = "euclidean" )$clustering
  end.time3=Sys.time()
  ari.wcm <- adjustedRandIndex(true.class, weight_clusterLabels)
```


The same can be done simply using the function `WeightedConsensusCluster`

```r
#one layer
  wcc <- WeightedConsensusCluster(exampleData, method="one layer", individualK = rep(3, 4),
                                   globalK = 3, pFeature = 0.8 ,ccClMethods = "kmeans",
                                    ccDistHCs = "euclidean",hclustMethod = "average",finalclmethod="hclust",
                                    finalhclustMethod = "average",Silhouette=TRUE)

#two-layer
  data_for_wcc <- list()
  data_for_wcc[[1]]=normData
  data_for_wcc[[2]]=binomData
  data_for_wcc[[3]]=lopoisData
  data_for_wcc[[4]]=multData
  data_for_wcc[[5]]=normData
  data_for_wcc[[6]]=binomData
  data_for_wcc[[7]]=poisData
  data_for_wcc[[8]]=multData
  wcc <- WeightedConsensusCluster(data_for_wcc, method="two-layer", individualK = rep(3, 8),
                     globalK = 3, pFeature = 0.8 ,ccClMethods = c("kmeans","kmeans","kmeans","kmeans",
                                                                  "hclust","hclust","hclust","hclust"),
                     ccDistHCs = c("euclidean","euclidean","euclidean","euclidean",
                                   "euclidean","binary","jaccard","hamming"),hclustMethod = "average",finalclmethod="hclust",
                                    finalhclustMethod = "average",Silhouette=TRUE)

```

Heatmap could be generated by plot_CM function.

```r
#Heatmap for consensus matrix
plot_CM(result=normRes)
#Heatmap for consensus matrix and dataset if dataset provided
plot_CM(data=normData,result=normRes)
```

## Session info
```{r}
sessionInfo()
```


## References 

Cabassi, A. and Kirk, P. D. W. (2019). Multiple kernel learning for integrative
consensus clustering of 'omic datasets. arXiv preprint. arXiv:1904.07701.
