#' Consensus clustering
#'
#' This function allows to perform consensus clustering using the k-means
#' clustering algorithm, for a fixed number of clusters. We consider the number
#' of clusters K to be fixed.
#' @param data N X P data matrix
#' @param K Number of clusters.
#' @param B Number of iterations.
#' @param pItem Proportion of items sampled at each iteration.
#' @param pFeature Proportion of features sampled at each iteration.
#' @param clMethod Clustering algorithm. Can be "hclust" for hierarchical
#' clustering, "kmeans" for k-means clustering, "pam" for partitioning around
#' medoids, "sparse-kmeans" for sparse k-means clustering or "sparse-hclust"
#' for sparse hierarchical clustering. Default is "hclust". However, if the data
#' contain at least one covariate that is a factor, the default clustering
#' algorithm is "pam".
#' @param dist Distance used for hierarchical clustering. Can be "pearson" (for
#' 1 - Pearson correlation), "spearman" (for 1- Spearman correlation), any of
#' the distances provided in stats::dist() (i.e. "euclidean", "maximum",
#' "manhattan", "canberra", "binary" or "minkowski"), or a matrix containing the
#' distances between the observations.
#' @param hclustMethod Hierarchical clustering method. Default is "average". For
#' more details see \code{?hclust}.
#' @param finalclmethod Final clustering method. Could be "pam" or "hclust".
#' @param finalhclustMethod  Hierarchical clustering method for final clustering. Default is "average". For
#' more details see \code{?hclust}.
#' @param sparseKmeansPenalty If the selected clustering method is
#' "sparse-kmeans", this is the value of the parameter "wbounds" of the
#' "KMeansSparseCluster" function. The default value is the square root of the
#' number of variables.
#' @param maxIterKM Number of iterations for the k-means clustering algorithm.
#' @param nstart nstart parameter for k-means if used.
#' @return consensus matrix and label results
#' @export
#'
#' @examples
consensuscluster <-
  function(data = NULL,
           K = 2,
           B = 100,
           pItem = 0.8, pFeature = 1,
           clMethod = "hclust",
           dist = "euclidean",
           hclustMethod = "average", finalclmethod = "hclust", finalhclustMethod = "average",
           sparseKmeansPenalty = NULL,
           maxIterKM = 1000, nstart = 20) {
    containsFactors <- 0
    if (!is.null(data)) {
      # Save number of observations
      N <- dim(data)[1]
      # Save number of covariates
      P <- dim(data)[2]
      # Check wheter there are factors among the covariates
      for (i in seq_len(P)) {
        containsFactors <-
          as.numeric(is.factor(data[, i])) + containsFactors
      }
    } else if (is.double(dist)) {
      N <- dim(dist)[1]
    } else {
      stop("If the data matrix is not provided, `dist` must be a symmetric
        matrix of type double providing the distances between each pair of
        observations.")
    }

    # Create vector of data indices that will be used to update the
    # co-clustering matrix
    dataIndices <- seq_len(N)
    # Initialise empty co-clustering matrix and an auxiliar
    coClusteringMatrix <- indicatorMatrix <- matrix(0, N, N)

    # For each step of the algorithm
    if (pFeature == 1) {
      for (b in seq_len(B)) {
        # Sample a proportion pItem of the observations without replacement
        items <- sample(N, ceiling(N * pItem), replace = FALSE)

        nUniqueDataPoints <- 0
        if (!is.null(data)) {
          # If there are more unique data points than clusters
          uniqueData <- unique(data[items, ])
          nUniqueDataPoints <- nrow(uniqueData)
        }

        if (nUniqueDataPoints > K | is.null(data)) {
          # If the chosen clustering method is PAM or there is at least one
          # covariate that is a factor
          if (clMethod == "pam" | containsFactors) {
            if (is.double(dist)) {
              distances <- stats::as.dist(dist[items, items])
            } else if (dist == "euclidean") {
              distances <- stats::dist(data[items, ], method = dist)
            } else if (dist == "cor") {
              distances <- stats::as.dist(1 - stats::cor(t(data[items, ])))
            } else if (dist == "binary") {
              distances <- stats::dist(data[items, ], method = dist)
            } else if (dist == "gower") {
              distances <- cluster::daisy(data[items, ], metric = "gower")
            } else {
              stop("Distance not recognized. If method is `pam`, distance
                must be one of `cor`, `binary`, `gower` or the symmetric
                matrix of distances.")
            }
            # Apply pam to the subsample and extract cluster labels
            cl <- cluster::pam(distances, K)$clustering
          } else if (clMethod == "kmeans" & !is.null(data)) {
            # Apply k-means to the subsample and extract cluster labels
            cl <- stats::kmeans(
              data[items, ], K,
              iter.max = maxIterKM, nstart = nstart
            )$cluster
          } else if (clMethod == "sparse-kmeans" & !is.null(data)) {
            if (is.null(sparseKmeansPenalty)) {
              km.perm <- KMeansSparseCluster.permute(data[items, ], K = K, nperms = 50)
            }

            # Apply sparse k-means to the subsample and extract cluster labels
            cl <- KMeansSparseCluster(
              data[items, ], K,
              wbounds = km.perm$bestw
            )[[1]]$Cs
          } else if (clMethod == "kmodes" & !is.null(data)) {
            # Apply k-means to the subsample and extract cluster labels
            cl <- NULL
            while (is.null(cl)) {
              try(cl <- klaR::kmodes(moc[items, ], 6, iter.max = 1000)$cluster, silent = TRUE)
            }
          } else if (clMethod == "hclust" | clMethod == "sparse-hclust") {
            if (dist == "pearson" | dist == "spearman") {
              pearsonCor <- stats::cor(t(data[items, ]), method = dist)
              distances <- stats::as.dist(1 - pearsonCor)
            } else if (is.double(dist)) {
              distances <- stats::as.dist(dist[items, items])
            } else if (dist == "hamming") {
              distances <- suppressMessages(stats::as.dist(e1071::hamming.distance(data[items, ])))
            } else if (dist == "binary") {
              distances <- stats::dist(data[items, ], method = dist)
            } else {
              # Calculate pairwise distances between observations
              distances <- suppressMessages(stats::as.dist(philentropy::distance(data[items, ], method = dist)))
            }
            # Apply hierarchical clustering to the subsample
            if (clMethod == "hclust") {
              hClustering <- stats::hclust(distances, method = hclustMethod)
            } else {
              hc.perm <- HierarchicalSparseCluster.permute(data[items, ], dissimilarity = "squared.distance", nperms = 50)
              hClustering <- HierarchicalSparseCluster(
                x = data[items, ], method = "average", dissimilarity = "squared.distance", wbound = hc.perm$bestw
              )$hc
            }
            cl <- stats::cutree(hClustering, K)
          } else {
            stop("Clustering algorithm name not recognised. Please choose
                     one of `kmeans`, `hclust`, `pam`, `sparse-kmeans`,
                     `sparse-hclust`.")
          }

          # Update matrix containing counts of number of times each pair has
          # been sampled together
          indicatorMatrix <- indicatorMatrix + crossprod(
            t(as.numeric(dataIndices %in% items))
          )

          # For each cluster
          for (k in seq_len(K)) {
            # Update counts of number of times each pair has been put in the
            # same cluster
            coClusteringMatrix[items, items] <-
              coClusteringMatrix[items, items] +
              crossprod(t(as.numeric(cl == k)))
          }
        }
      }
    } else {
      for (b in seq_len(B)) {
        # Sample a proportion pItem of the observations without replacement
        items <- sample(N, ceiling(N * pItem), replace = FALSE)
        features <- sample(P, ceiling(P * pFeature), replace = FALSE)
        nUniqueDataPoints <- 0
        if (!is.null(data)) {
          # If there are more unique data points than clusters
          uniqueData <- unique(data[items, ])
          nUniqueDataPoints <- nrow(uniqueData)
        }

        if (nUniqueDataPoints > K | is.null(data)) {
          # If the chosen clustering method is PAM or there is at least one
          # covariate that is a factor
          if (clMethod == "pam" | containsFactors) {
            if (is.double(dist)) {
              distances <- stats::as.dist(dist[items, items])
            } else if (dist == "euclidean") {
              distances <- stats::dist(data[items, features], method = dist)
            } else if (dist == "cor") {
              distances <- stats::as.dist(1 - stats::cor(t(data[items, features])))
            } else if (dist == "binary") {
              distances <- stats::dist(data[items, features], method = dist)
            } else if (dist == "gower") {
              distances <- cluster::daisy(data[items, features], metric = "gower")
            } else {
              stop("Distance not recognized. If method is `pam`, distance
                must be one of `cor`, `binary`, `gower` or the symmetric
                matrix of distances.")
            }
            # Apply pam to the subsample and extract cluster labels
            cl <- cluster::pam(distances, K)$clustering
          } else if (clMethod == "kmeans" & !is.null(data)) {
            # Apply k-means to the subsample and extract cluster labels
            cl <- stats::kmeans(
              data[items, features], K,
              iter.max = maxIterKM, nstart = nstart
            )$cluster
          } else if (clMethod == "sparse-kmeans" & !is.null(data)) {
            if (is.null(sparseKmeansPenalty)) {
              km.perm <- KMeansSparseCluster.permute(data[items, features], K = K, nperms = 50)
            }

            # Apply sparse k-means to the subsample and extract cluster labels
            cl <- KMeansSparseCluster(
              data[items, features], K,
              wbounds = km.perm$bestw
            )[[1]]$Cs
          } else if (clMethod == "kmodes" & !is.null(data)) {
            # Apply k-means to the subsample and extract cluster labels
            cl <- NULL
            while (is.null(cl)) {
              try(cl <- klaR::kmodes(moc[items, ], 6, iter.max = 1000)$cluster, silent = TRUE)
            }
          } else if (clMethod == "hclust" | clMethod == "sparse-hclust") {
            if (dist == "pearson" | dist == "spearman") {
              pearsonCor <- stats::cor(t(data[items, features]), method = dist)
              distances <- stats::as.dist(1 - pearsonCor)
            } else if (is.double(dist)) {
              distances <- stats::as.dist(dist[items, items])
            } else if (dist == "hamming") {
              distances <- suppressMessages(stats::as.dist(e1071::hamming.distance(data[items, features])))
            } else if (dist == "binary") {
              distances <- stats::dist(data[items, features], method = dist)
            } else {
              # Calculate pairwise distances between observations
              distances <- suppressMessages(stats::as.dist(philentropy::distance(data[items, features], method = dist)))
            }
            # Apply hierarchical clustering to the subsample
            if (clMethod == "hclust") {
              hClustering <- stats::hclust(distances, method = hclustMethod)
            } else {
              hc.perm <- HierarchicalSparseCluster.permute(data[items, features], dissimilarity = "squared.distance", nperms = 50)
              hClustering <- HierarchicalSparseCluster(
                x = data[items, features], method = "average", dissimilarity = "squared.distance", wbound = hc.perm$bestw
              )$hc
            }
            cl <- stats::cutree(hClustering, K)
          } else {
            stop("Clustering algorithm name not recognised. Please choose
                     one of `kmeans`, `hclust`, `pam`, `sparse-kmeans`,
                     `sparse-hclust`.")
          }

          # Update matrix containing counts of number of times each pair has
          # been sampled together
          indicatorMatrix <- indicatorMatrix + crossprod(
            t(as.numeric(dataIndices %in% items))
          )

          # For each cluster
          for (k in seq_len(K)) {
            # Update counts of number of times each pair has been put in the
            # same cluster
            coClusteringMatrix[items, items] <-
              coClusteringMatrix[items, items] +
              crossprod(t(as.numeric(cl == k)))
          }
        }
      }
    }

    if (!sum(indicatorMatrix) == 0) {
      # Normalise co-clustering matrix
      consensusMatrix <- coClusteringMatrix / indicatorMatrix
    } else {
      consensusMatrix <- indicatorMatrix
      warning(paste("Consensus matrix is empty for K =", K, "because there are
                      less than", K, "distinct data points", sep = ""))
    }

    distances <- stats::as.dist(1 - consensusMatrix)
    if (finalclmethod == "hclust") {
      hClustering <- stats::hclust(distances, method = finalhclustMethod)
      clusterLabels <- stats::cutree(hClustering, K)
    } else {
      clusterLabels <- pam(distances, K, diss = TRUE, metric = "euclidean")$clustering
    }
    res <- alist()
    res[[1]] <- consensusMatrix
    res[[2]] <- clusterLabels
    names(res) <- c("consensusMatrix", "class")
    return(res)
  }
