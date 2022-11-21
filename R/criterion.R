#' Title
#'
#' @param CM A list of consensus matrix.
#' @param clLabels A matrix consists of clustering result corresponding to each consensus matrix.
#' @param K A vector denotes number of clusters for each consensus matrix.
#' @param Silhouette Boolean. If TRUE, compute Silhouette to choose
#' best number of clusters. Default is TRUE.
#' @param widestGap Boolean. If TRUE, compute also widest gap index to choose
#' best number of clusters. Default is FALSE.
#' @param dunns Boolean. If TRUE, compute also Dunn's index to choose best
#' number of clusters. Default is FALSE.
#' @param dunn2s Boolean. If TRUE, compute also alternative Dunn's index to
#' choose best number of clusters. Default is FALSE.
#'
#' @return best final K and criterion used
#' @export
#' @examples library(weightedCC)
#' simulation=data_generation(20,5,30)
#' normData=simulation[[1]]
#' tempCM <- array(NA, c(60, 60, 4))
#' clLabels <- matrix(NA, 4, 60)
#' for (k in 2:5){
#' CCoutput=consensuscluster(normData,K=k,B=1000,pItem = 0.8,pFeature = 0.8 ,
#' clMethod ="kmeans",finalclmethod="pam")
#'tempCM[, , k - 1] <- CCoutput[[1]]
#' clLabels[k - 1, ] <- CCoutput[[2]]
#'
#' }
#' chooseK <- criterion(tempCM, clLabels,K = 2:5, Silhouette = TRUE)


criterion <- function(CM, clLabels, K, Silhouette = Silhouette, widestGap = widestGap,
                      dunns = dunns,
                      dunn2s = dunn2s) {
  output <- list()
  n <- length(K)
  criterion <- rep(NA, n)
  if (Silhouette == TRUE) {
    for (i in 1:n) {
      criterion[i] <- mean(cluster::silhouette(clLabels[i, ], stats::as.dist(1 - CM[, , i]))[, "sil_width"])
    }
  } else {
    if (widestGap == TRUE) {
      for (i in 1:n) {
        criterion[i] <- fpc::cluster.stats(stats::as.dist(1 - CM[, , i]), clLabels[i, ])$widestgap
      }
    }

    if (dunns == TRUE) {
      for (i in 1:n) {
        criterion[i] <- fpc::cluster.stats(stats::as.dist(1 - CM[, , i]), clLabels[i, ])$dunn
      }
    }

    if (dunn2s == TRUE) {
      for (i in 1:n) {
        criterion[i] <- fpc::cluster.stats(stats::as.dist(1 - CM[, , i]), clLabels[i, ])$dunn2
      }
    }
  }

  output$finalK <- K[which.max(criterion)]
  output$criterion <- criterion

  return(output)
}
