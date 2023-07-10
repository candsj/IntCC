#' Consensus clustering
#'
#' This function allows to perform consensus clustering using the k-means
#' clustering algorithm, for a fixed number of clusters. We consider the number
#' of clusters K to be fixed.
#' @param data N X P data matrix
#' @param result Results derived from consensus clustering.
#' It consists of a list of several consensus matrix and corresponding clustering label results.
#' @return heatmap for consensus matrix from the result or heatmap for dataset
#' @export
#'
#' @examples library(weightedCC)
#' load(system.file("extdata", "exampleData.RData", package = "weightedCC"))
#' normData=exampleData[[1]]
#' normRes=consensuscluster(normData,K=3,B=1000,pItem = 0.8,pFeature = 0.8,
#' clMethod ="kmeans",finalclmethod="pam")
#' plot_CM(result=normRes)
#' plot_CM(data=normData,result=normRes)


plot_CM <-
  function(data = NULL,result=NULL) {
    if (is.null(result)) {stop("Result is not provided!")}
    CM <- result[[1]]
    clusterID <- result[[2]]
    cluster_order=unlist(split(seq_along(clusterID), clusterID))
    pheatmap::pheatmap(CM[cluster_order,cluster_order],treeheight_row = 0,treeheight_col = 0)
    if (!is.null(data)) {pheatmap::pheatmap(data[cluster_order,],treeheight_row = 0,treeheight_col = 0)}
  }
