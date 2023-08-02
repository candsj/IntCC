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
#' @examples library(intCC)
#' load(system.file("extdata", "exampleData.RData", package = "intCC"))
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
    k=length(unique(clusterID))
    cluster_order=unlist(split(seq_along(clusterID), clusterID))
    sample_annotation=data.frame(Sample=rep(paste0("Cluster ",1:k),lengths(split(seq_along(clusterID), clusterID))))
    CM_order=CM[cluster_order,cluster_order]
    rownames(CM_order)=rownames(sample_annotation)=cluster_order
    pheatmap::pheatmap(CM_order,annotation_row = sample_annotation,cluster_rows=FALSE,cluster_cols=FALSE,
                       show_rownames=FALSE,treeheight_row = 0,treeheight_col = 0)
    if (!is.null(data)) {
      data_order=data[cluster_order,]
      rownames(data_order)=cluster_order
      pheatmap::pheatmap(data_order,annotation_row = sample_annotation,cluster_rows=FALSE,cluster_cols=FALSE,
                         show_rownames=FALSE,treeheight_row = 0,treeheight_col = 0)
      }
    }

