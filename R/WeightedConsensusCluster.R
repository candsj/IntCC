#' Title
#'
#' @param data List of M datasets, each of size N X P_m, m = 1, ..., M.
#' @param method One layer or two-layer. Default is one layer.
#' @param individualK Vector containing the number of clusters in each dataset.
#' Default is NULL. If the number of clusters is not provided, then all the
#' possible values between 2 and individualMaxK are considered and the best
#' value is chosen for each dataset by maximising the silhouette.
#' @param individualMaxK Maximum number of clusters considered for the
#' individual data. Default is 10.
#' @param globalK Number of global clusters. Default is NULL. If the number of
#' clusters is not provided, then all the number of clusters for each dataset will be used and decided
#' according to criterion defined.
#' @param B Number of iterations for consensus clustering. Default is 1000.
#' @param pItem Proportion of items sampled at each iteration.
#' @param pFeature Proportion of features sampled at each iteration.
#' @param ccClMethods The i-th element of this vector goes into the
#' "clMethod" argument of consensusCluster() for the i-th dataset. If only
#' one string is provided, then the same method is used for all datasets.
#' @param ccDistHCs Distance used for each clustering method. Default is Euclidean distance.
#' @param hclustMethod Method used for hclust method only worked when hclust method used in consensus clustering.
#' @param finalclmethod Method used for final clustering. It could be pam or hclust.
#' @param finalhclustMethod Method used for hclust method only worked when hclust method used in final clustering.
#' @param Silhouette Boolean. If TRUE, compute Silhouette to choose
#' best number of clusters. Default is TRUE.
#' @param widestGap Boolean. If TRUE, compute also widest gap index to choose
#' best number of clusters. Default is FALSE.
#' @param dunns Boolean. If TRUE, compute also Dunn's index to choose best
#' number of clusters. Default is FALSE.
#' @param dunn2s Boolean. If TRUE, compute also alternative Dunn's index to
#' choose best number of clusters. Default is FALSE.
#' @return clustering results, weights and weighted consensus matrix
#' @export
#' @examples
#' library(weightedCC)
#' load(system.file("extdata", "exampleData.RData", package = "weightedCC"))
#' wcc <- WeightedConsensusCluster(exampleData, method="one layer", individualK = rep(3,4),
#' globalK = 3, pFeature = 0.8 ,ccClMethods = "kmeans",
#' ccDistHCs = "euclidean",hclustMethod = "average",finalclmethod="hclust",
#' finalhclustMethod = "average",Silhouette=TRUE)
WeightedConsensusCluster <- function(data, method = NULL, individualK = NULL,
                                     individualMaxK = 10,
                                     globalK = NULL,
                                     B = 1000, pItem = 0.8,
                                     pFeature = 1, ccClMethods = "kmeans",
                                     ccDistHCs = "euclidean", hclustMethod = "average", finalclmethod = "hclust",
                                     finalhclustMethod = "average",
                                     Silhouette = TRUE,
                                     widestGap = FALSE,
                                     dunns = FALSE,
                                     dunn2s = FALSE) {
  if (is.null(data)) stop("data should be provided")
  N <- dim(data[[1]])[1]
  M <- length(data)
  for (i in 1:M) {
    if (dim(data[[i]])[1] != N) {
      stop("All datasets must have the same number of rows.")
    }
  }
  if (is.null(method)) stop("method should be provided")
  method <- match.arg(method, c("one layer", "two-layer"))

  res.all <- list()
  output <- list()
  output$bestK <- rep(NA, M)
  CM <- array(NA, c(N, N, M))

  if (is.null(individualK)) {
    if (method == "one layer") {
      if (length(ccClMethods) == 1) {
        ccClMethods <- rep(ccClMethods, M)
      } else if (length(ccClMethods) != M) {
        stop("Please specify a method for each dataset by passing a vector of length", M, "to ccClMethods.")
      }

      if (length(ccDistHCs) == 1) {
        ccDistHCs <- rep(ccDistHCs, M)
      } else if (length(ccDistHCs) != M) {
        stop("Please specify a distance metric for each dataset by passing a vector of length", M, "to ccDistHCs.")
      }

      output$bestK <- rep(NA, M)
      # Initialise empty consensus matrices for all possible numbers of
      # clusters
      tempCM <- array(NA, c(N, N, individualMaxK - 1))
      # Initialise empty cluster labels for all possible numbers of clusters
      clLabels <- matrix(NA, individualMaxK - 1, N)

      for (i in 1:M) {
        dataset_i <- data[[i]]
        ccClMethod_i <- ccClMethods[i]
        ccDistHC_i <- ccDistHCs[i]
        for (j in 2:individualMaxK) {
          CCoutput <-
            consensuscluster(dataset_i, j, B, pItem, pFeature,
              clMethod = ccClMethod_i,
              dist = ccDistHC_i, finalclmethod
            )

          tempCM[, , j - 1] <- CCoutput[[1]]

          clLabels[j - 1, ] <- CCoutput[[2]]
        }
        chooseK <- criterion(tempCM, clLabels,
          K = 2:individualMaxK, Silhouette = Silhouette, widestGap = widestGap,
          dunns = dunns, dunn2s = dunn2s
        )
        # If there is more than one, choose smallest number of clusters
        # among the ones that maximise the silhouette
        bestK <- output$bestK[i] <- chooseK$finalK
        res.all[[i]] <- list(tempCM[, , bestK - 1], clLabels[bestK - 1, ])
        names(res.all[[i]]) <- c("consensusMatrix", "class")
      }

      weights <- weightcal(res.all)
      wcm <- 0
      for (l in 1:M) {
        wcm <- wcm + weights[l] * res.all[[l]]$consensusMatrix
      }
      distances <- stats::as.dist(1 - wcm)

      if (is.null(globalK)) {
        Ks <- sort(unique(output$bestK))
        if (length(Ks) == 1) {
          if (finalclmethod == "pam") {
            weight_clusterLabels <- cluster::pam(distances, Ks, diss = TRUE, metric = "euclidean")$clustering
          } else {
            hClustering <- stats::hclust(distances, method = finalhclustMethod)
            weight_clusterLabels <- stats::cutree(hClustering, Ks)
          }
          output$globalK <- Ks
          # Save chosen cluster labels
          output$globalClusterLabels <- weight_clusterLabels
          # Save chosen weights
          output$weights <- weights
          # Save chosen consensus matrix
          output$weightedKM <- wcm
        } else {
          number.Ks <- length(Ks)
          tempCM <- array(NA, c(N, N, number.Ks))
          # Initialise empty cluster labels for all possible numbers of clusters
          clLabels <- matrix(NA, number.Ks, N)

          for (k in 1:number.Ks) {
            if (finalclmethod == "pam") {
              weight_clusterLabels <- cluster::pam(distances, Ks[k], diss = TRUE, metric = "euclidean")$clustering
              tempCM[, , k] <- wcm
              clLabels[k, ] <- weight_clusterLabels
            } else {
              hClustering <- stats::hclust(distances, method = finalhclustMethod)
              weight_clusterLabels <- stats::cutree(hClustering, Ks[k])
            }
            tempCM[, , k] <- wcm
            clLabels[k, ] <- weight_clusterLabels
          }
          chooseK <- criterion(tempCM, clLabels,
            K = Ks, Silhouette = Silhouette, widestGap = widestGap,
            dunns = dunns, dunn2s = dunn2s
          )
          globalK <- output$globalK <- chooseK$finalK
          # Save chosen cluster labels
          output$globalClusterLabels <- clLabels[which(Ks == globalK), ]
          # Save chosen weights
          output$weights <- weights
          # Save chosen consensus matrix
          output$weightedKM <- wcm
        }
      } else {
        if (finalclmethod == "pam") {
          weight_clusterLabels <- cluster::pam(distances, globalK, diss = TRUE, metric = "euclidean")$clustering
        } else {
          hClustering <- stats::hclust(distances, method = finalhclustMethod)
          weight_clusterLabels <- stats::cutree(hClustering, globalK)
        }
        output$globalK <- globalK
        # Save chosen cluster labels
        output$globalClusterLabels <- weight_clusterLabels
        # Save chosen weights
        output$weights <- weights
        # Save chosen consensus matrix
        output$weightedKM <- wcm
      }
    } else {
      if (length(ccClMethods) == 1) {
        ccClMethods <- vector("list", M)
        for (i in 1:M) {
          ccClMethods[[i]] <- c("kmeans", "hclust")
        }
      } else if (length(ccClMethods) != M) {
        stop("Please specify methods for each dataset by passing a list of length", M, "to ccClMethods.")
      }

      if (length(ccDistHCs) == 1) {
        ccDistHCs <- rep(ccDistHCs, M)
        for (i in 1:M) {
          ccDistHCs[[i]] <- c("euclidean", "euclidean")
        }
      } else if (length(ccDistHCs) != M) {
        stop("Please specify distance metrics for each dataset by passing a list of length", M, "to ccDistHCs.")
      }

      number.methods <- rep(NA, M)
      for (i in 1:M) {
        if (length(ccDistHCs[[i]]) != length(ccClMethods[[i]])) {
          stop("All datasets must have the same number of methods and distance metrics.")
        }
        number.methods[i] <- length(ccDistHCs[[i]])
      }


      output$bestK <- vector("list", M)
      # Initialise empty consensus matrices for all possible numbers of
      # clusters
      tempCM <- array(NA, c(N, N, individualMaxK - 1))

      # Initialise empty cluster labels for all possible numbers of clusters
      clLabels <- matrix(NA, individualMaxK - 1, N)


      for (i in 1:M) {
        dataset_i <- data[[i]]
        if (number.methods[i] == 1) {
          ccClMethod_i <- ccClMethods[[i]]
          ccDistHC_i <- ccDistHCs[[i]]
          for (j in 2:individualMaxK) {
            CCoutput <-
              consensuscluster(dataset_i, j, B, pItem, pFeature,
                clMethod = ccClMethod_i,
                dist = ccDistHC_i, finalclmethod
              )

            tempCM[, , j - 1] <- CCoutput[[1]]

            clLabels[j - 1, ] <- CCoutput[[2]]
          }
          chooseK <- criterion(tempCM, clLabels,
            K = 2:individualMaxK, Silhouette = Silhouette, widestGap = widestGap,
            dunns = dunns, dunn2s = dunn2s
          )
          # If there is more than one, choose smallest number of clusters
          # among the ones that maximise the silhouette
          bestK <- output$bestK[[i]] <- chooseK$finalK
          output$weights.methods[[i]] <- 1
          output$weightedKM.methods[[i]] <- tempCM[, , bestK - 1]
          res.all[[i]] <- list(tempCM[, , bestK - 1], clLabels[bestK - 1, ])
          names(res.all[[i]]) <- c("consensusMatrix", "class")
        } else {
          weight.res <- vector("list", number.methods[i])
          for (l in 1:number.methods[i]) {
            ccClMethod_i <- ccClMethods[[i]][l]
            ccDistHC_i <- ccDistHCs[[i]][l]
            for (j in 2:individualMaxK) {
              CCoutput <-
                consensuscluster(dataset_i, j, B, pItem, pFeature,
                  clMethod = ccClMethod_i,
                  dist = ccDistHC_i, finalclmethod
                )

              tempCM[, , j - 1] <- CCoutput[[1]]

              clLabels[j - 1, ] <- CCoutput[[2]]
            }
            chooseK <- criterion(tempCM, clLabels,
              K = 2:individualMaxK, Silhouette = Silhouette, widestGap = widestGap,
              dunns = dunns, dunn2s = dunn2s
            )
            # If there is more than one, choose smallest number of clusters
            # among the ones that maximise the silhouette
            bestK <- output$bestK[[i]][l] <- chooseK$finalK
            weight.res[[l]] <- list(tempCM[, , bestK - 1], clLabels[bestK - 1, ])
          }
          weights <- weightcal(weight.res)
          wcm <- 0
          for (k in 1:number.methods[i]) {
            wcm <- wcm + weights[k] * weight.res[[k]]$consensusMatrix
          }
          distances <- stats::as.dist(1 - wcm)
          # Save chosen consensus matrix
          output$weightedKM.methods[[i]] <- wcm

          Ks <- sort(output$bestK[[i]])
          if (length(Ks) == 1) {
            if (finalclmethod == "pam") {
              weight_clusterLabels <- cluster::pam(distances, Ks, diss = TRUE, metric = "euclidean")$clustering
            } else {
              hClustering <- stats::hclust(distances, method = finalhclustMethod)
              weight_clusterLabels <- stats::cutree(hClustering, Ks)
            }

            output$bestK[[i]] <- Ks
            output$weights.methods[[i]] <- weights
            res.all[[i]] <- list(wcm, weight_clusterLabels)
            names(res.all[[i]]) <- c("consensusMatrix", "class")
          } else {
            number.Ks <- length(Ks)
            tempCM <- array(NA, c(N, N, number.Ks))
            # Initialise empty cluster labels for all possible numbers of clusters
            clLabels <- matrix(NA, number.Ks, N)

            for (k in 1:number.Ks) {
              if (finalclmethod == "pam") {
                weight_clusterLabels <- cluster::pam(distances, Ks[k], diss = TRUE, metric = "euclidean")$clustering
                tempCM[, , k] <- wcm
                clLabels[k, ] <- weight_clusterLabels
              } else {
                hClustering <- stats::hclust(distances, method = finalhclustMethod)
                weight_clusterLabels <- stats::cutree(hClustering, Ks[k])
              }
              tempCM[, , k] <- wcm
              clLabels[k, ] <- weight_clusterLabels
            }
            chooseK <- criterion(tempCM, clLabels,
              K = Ks, Silhouette = Silhouette, widestGap = widestGap,
              dunns = dunns, dunn2s = dunn2s
            )
            output$bestK[[i]] <- chooseK$finalK
            output$weights.methods[[i]] <- weights



            res.all[[i]] <- list(wcm, clLabels[which(Ks == output$bestK[[i]]), ])
            names(res.all[[i]]) <- c("consensusMatrix", "class")
          }
        }
      }

      if (is.null(globalK)) {
        Ks <- sort(unique(unlist(output$bestK)))
        if (length(Ks) == 1) {
          if (finalclmethod == "pam") {
            weight_clusterLabels <- cluster::pam(distances, Ks, diss = TRUE, metric = "euclidean")$clustering
          } else {
            hClustering <- stats::hclust(distances, method = finalhclustMethod)
            weight_clusterLabels <- stats::cutree(hClustering, Ks)
          }
          output$globalK <- Ks
          # Save chosen cluster labels
          output$globalClusterLabels <- weight_clusterLabels
          # Save chosen weights
          output$weights <- weights
          # Save chosen consensus matrix
          output$weightedKM <- wcm
        } else {
          number.Ks <- length(Ks)
          tempCM <- array(NA, c(N, N, number.Ks))
          # Initialise empty cluster labels for all possible numbers of clusters
          clLabels <- matrix(NA, number.Ks, N)

          for (k in 1:number.Ks) {
            if (finalclmethod == "pam") {
              weight_clusterLabels <- cluster::pam(distances, Ks[k], diss = TRUE, metric = "euclidean")$clustering
              tempCM[, , k] <- wcm
              clLabels[k, ] <- weight_clusterLabels
            } else {
              hClustering <- stats::hclust(distances, method = finalhclustMethod)
              weight_clusterLabels <- stats::cutree(hClustering, Ks[k])
            }
            tempCM[, , k] <- wcm
            clLabels[k, ] <- weight_clusterLabels
          }
          chooseK <- criterion(tempCM, clLabels,
            K = Ks, Silhouette = Silhouette, widestGap = widestGap,
            dunns = dunns, dunn2s = dunn2s
          )
          globalK <- output$globalK <- chooseK$finalK
          # Save chosen cluster labels
          output$globalClusterLabels <- clLabels[which(Ks == globalK), ]
          # Save chosen weights
          output$weights <- weights
          # Save chosen consensus matrix
          output$weightedKM <- wcm
        }
      } else {
        if (finalclmethod == "pam") {
          weight_clusterLabels <- cluster::pam(distances, globalK, diss = TRUE, metric = "euclidean")$clustering
        } else {
          hClustering <- stats::hclust(distances, method = finalhclustMethod)
          weight_clusterLabels <- stats::cutree(hClustering, globalK)
        }
        output$globalK <- globalK
        # Save chosen cluster labels
        output$globalClusterLabels <- weight_clusterLabels
        # Save chosen weights
        output$weights <- weights
        # Save chosen consensus matrix
        output$weightedKM <- wcm
      }
    }
  } else if (method == "one layer") {
    if (length(individualK) == 1) {
      individualK <- rep(individualK, M)
    } else if (length(individualK) != M)
    {
      stop("Please specify the cluster number for each dataset by passing a vector of length", M, "to individualK.")
    }

    if (length(ccClMethods) == 1) {
      ccClMethods <- rep(ccClMethods, M)
    } else if (length(ccClMethods) != M) {
      stop("Please specify a method for each dataset by passing a vector of length", M, "to ccClMethods.")
    }

    if (length(ccDistHCs) == 1) {
      ccDistHCs <- rep(ccDistHCs, M)
    } else if (length(ccDistHCs) != M) {
      stop("Please specify a distance metric for each dataset by passing a vector of length", M, "to ccDistHCs.")
    }

    output$bestK <- rep(NA, M)
    # Initialise empty consensus matrices for all possible numbers of
    # clusters
    tempCM <- array(NA, c(N, N, individualMaxK - 1))
    # Initialise empty cluster labels for all possible numbers of clusters
    clLabels <- matrix(NA, individualMaxK - 1, N)

    for (i in 1:M) {
      dataset_i <- data[[i]]
      ccClMethod_i <- ccClMethods[i]
      ccDistHC_i <- ccDistHCs[i]

      res.all[[i]] <- consensuscluster(dataset_i, individualK[i], B, pItem, pFeature,
        clMethod = ccClMethod_i,
        dist = ccDistHC_i, finalclmethod
      )
      names(res.all[[i]]) <- c("consensusMatrix", "class")
    }

    weights <- weightcal(res.all)
    wcm <- 0
    for (l in 1:M) {
      wcm <- wcm + weights[l] * res.all[[l]]$consensusMatrix
    }
    distances <- stats::as.dist(1 - wcm)

    if (is.null(globalK)) {
      Ks <- sort(unique(individualK))
      if (length(Ks) == 1) {
        if (finalclmethod == "pam") {
          weight_clusterLabels <- cluster::pam(distances, Ks, diss = TRUE, metric = "euclidean")$clustering
        } else {
          hClustering <- stats::hclust(distances, method = finalhclustMethod)
          weight_clusterLabels <- stats::cutree(hClustering, Ks)
        }
        output$globalK <- Ks
        # Save chosen cluster labels
        output$globalClusterLabels <- weight_clusterLabels
        # Save chosen weights
        output$weights <- weights
        # Save chosen consensus matrix
        output$weightedKM <- wcm
      } else {
        number.Ks <- length(Ks)
        tempCM <- array(NA, c(N, N, number.Ks))
        # Initialise empty cluster labels for all possible numbers of clusters
        clLabels <- matrix(NA, number.Ks, N)

        for (k in 1:number.Ks) {
          if (finalclmethod == "pam") {
            weight_clusterLabels <- cluster::pam(distances, Ks[k], diss = TRUE, metric = "euclidean")$clustering
            tempCM[, , k] <- wcm
            clLabels[k, ] <- weight_clusterLabels
          } else {
            hClustering <- stats::hclust(distances, method = finalhclustMethod)
            weight_clusterLabels <- stats::cutree(hClustering, Ks[k])
          }
          tempCM[, , k] <- wcm
          clLabels[k, ] <- weight_clusterLabels
        }
        chooseK <- criterion(tempCM, clLabels,
          K = Ks, Silhouette = Silhouette, widestGap = widestGap,
          dunns = dunns, dunn2s = dunn2s
        )
        globalK <- output$globalK <- chooseK$finalK
        # Save chosen cluster labels
        output$globalClusterLabels <- clLabels[which(Ks == globalK), ]
        # Save chosen weights
        output$weights <- weights
        # Save chosen consensus matrix
        output$weightedKM <- wcm
      }
    } else {
      if (finalclmethod == "pam") {
        weight_clusterLabels <- cluster::pam(distances, globalK, diss = TRUE, metric = "euclidean")$clustering
      } else {
        hClustering <- stats::hclust(distances, method = finalhclustMethod)
        weight_clusterLabels <- stats::cutree(hClustering, globalK)
      }
      output$globalK <- globalK
      # Save chosen cluster labels
      output$globalClusterLabels <- weight_clusterLabels
      # Save chosen weights
      output$weights <- weights
      # Save chosen consensus matrix
      output$weightedKM <- wcm
    }
  } else {
    number.methods <- rep(NA, M)

    if (length(individualK) == 1) {
      individualK <- rep(individualK, M)
    } else if (length(individualK) != M)
    {
      stop("Please specify the cluster number for each dataset by passing a vector of length", M, "to individualK.")
    }

    for (i in 1:M) {
      if (length(ccDistHCs[[i]]) != length(ccClMethods[[i]])) {
        stop("All datasets must have the same number of methods and distance metrics.")
      }
      number.methods[i] <- length(ccDistHCs[[i]])
    }

    if (length(ccClMethods) == 1) {
      ccClMethods <- vector("list", M)
      for (i in 1:M) {
        ccClMethods[[i]] <- c("kmeans", "hclust")
      }
    } else if (length(ccClMethods) != M) {
      stop("Please specify methods for each dataset by passing a list of length", M, "to ccClMethods.")
    }

    if (length(ccDistHCs) == 1) {
      ccDistHCs <- rep(ccDistHCs, M)
      for (i in 1:M) {
        ccDistHCs[[i]] <- c("euclidean", "euclidean")
      }
    } else if (length(ccDistHCs) != M) {
      stop("Please specify distance metrics for each dataset by passing a list of length", M, "to ccDistHCs.")
    }




    output$bestK <- rep(NA, M)
    # Initialise empty consensus matrices for all possible numbers of
    # clusters
    tempCM <- array(NA, c(N, N, individualMaxK - 1))

    # Initialise empty cluster labels for all possible numbers of clusters
    clLabels <- matrix(NA, individualMaxK - 1, N)


    for (i in 1:M) {
      dataset_i <- data[[i]]
      if (number.methods[i] == 1) {
        ccClMethod_i <- ccClMethods[[i]]
        ccDistHC_i <- ccDistHCs[[i]]

        bestK <- output$bestK[i] <- individualK[i]
        output$weights.methods[[i]] <- 1

        res.all[[i]] <- consensuscluster(dataset_i, individualK[i], B, pItem, pFeature,
          clMethod = ccClMethod_i,
          dist = ccDistHC_i, finalclmethod
        )
        names(res.all[[i]]) <- c("consensusMatrix", "class")
        output$weightedKM.methods[[i]] <- res.all[[i]]$consensusMatrix
      } else {
        weight.res <- vector("list", number.methods[i])
        for (l in 1:number.methods[i]) {
          ccClMethod_i <- ccClMethods[[i]][l]
          ccDistHC_i <- ccDistHCs[[i]][l]

          bestK <- output$bestK[i] <- individualK[i]
          weight.res[[l]] <- consensuscluster(dataset_i, individualK[i], B, pItem, pFeature,
            clMethod = ccClMethod_i,
            dist = ccDistHC_i, finalclmethod
          )
        }
        weights <- weightcal(weight.res)
        wcm <- 0
        for (k in 1:number.methods[i]) {
          wcm <- wcm + weights[k] * weight.res[[k]]$consensusMatrix
        }
        distances <- stats::as.dist(1 - wcm)
        # Save chosen consensus matrix
        output$weightedKM.methods[[i]] <- wcm


        if (finalclmethod == "pam") {
          weight_clusterLabels <- cluster::pam(distances, individualK[i], diss = TRUE, metric = "euclidean")$clustering
        } else {
          hClustering <- stats::hclust(distances, method = finalhclustMethod)
          weight_clusterLabels <- stats::cutree(hClustering, individualK[i])
        }

        output$weights.methods[[i]] <- weights
        res.all[[i]] <- list(wcm, weight_clusterLabels)
        names(res.all[[i]]) <- c("consensusMatrix", "class")
      }
    }

    if (is.null(globalK)) {
      Ks <- sort(unique(unlist(output$bestK)))
      if (length(Ks) == 1) {
        if (finalclmethod == "pam") {
          weight_clusterLabels <- cluster::pam(distances, Ks, diss = TRUE, metric = "euclidean")$clustering
        } else {
          hClustering <- stats::hclust(distances, method = finalhclustMethod)
          weight_clusterLabels <- stats::cutree(hClustering, Ks)
        }
        output$globalK <- Ks
        # Save chosen cluster labels
        output$globalClusterLabels <- weight_clusterLabels
        # Save chosen weights
        output$weights <- weights
        # Save chosen consensus matrix
        output$weightedKM <- wcm
      } else {
        number.Ks <- length(Ks)
        tempCM <- array(NA, c(N, N, number.Ks))
        # Initialise empty cluster labels for all possible numbers of clusters
        clLabels <- matrix(NA, number.Ks, N)

        for (k in 1:number.Ks) {
          if (finalclmethod == "pam") {
            weight_clusterLabels <- cluster::pam(distances, Ks[k], diss = TRUE, metric = "euclidean")$clustering
            tempCM[, , k] <- wcm
            clLabels[k, ] <- weight_clusterLabels
          } else {
            hClustering <- stats::hclust(distances, method = finalhclustMethod)
            weight_clusterLabels <- stats::cutree(hClustering, Ks[k])
          }
          tempCM[, , k] <- wcm
          clLabels[k, ] <- weight_clusterLabels
        }
        chooseK <- criterion(tempCM, clLabels,
          K = Ks, Silhouette = Silhouette, widestGap = widestGap,
          dunns = dunns, dunn2s = dunn2s
        )
        globalK <- output$globalK <- chooseK$finalK
        # Save chosen cluster labels
        output$globalClusterLabels <- clLabels[which(Ks == globalK), ]
        # Save chosen weights
        output$weights <- weights
        # Save chosen consensus matrix
        output$weightedKM <- wcm
      }
    } else {
      if (finalclmethod == "pam") {
        weight_clusterLabels <- cluster::pam(distances, globalK, diss = TRUE, metric = "euclidean")$clustering
      } else {
        hClustering <- stats::hclust(distances, method = finalhclustMethod)
        weight_clusterLabels <- stats::cutree(hClustering, globalK)
      }
      output$globalK <- globalK
      # Save chosen cluster labels
      output$globalClusterLabels <- weight_clusterLabels
      # Save chosen weights
      output$weights <- weights
      # Save chosen consensus matrix
      output$weightedKM <- wcm
    }
  }
  return(output)
}
