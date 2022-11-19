#' Title
#'
#' @param result Results derived from consensus clustering.
#' It consists of a list of several consensus matrix and corresponding clustering label results.
#'
#' @return weights for each consensus matrix
#' @export
#'
weightcal <- function(result) {
  n <- length(result)
  clusterID <- sapply(result, "[", 2)
  CCmat <- sapply(result, "[", 1)

  CCmat.sort <- lapply(1:n, FUN = function(n) {
    CCmat[[n]][order(clusterID[[n]]), order(clusterID[[n]])]
  })
  clusterID.sort <- lapply(1:n, FUN = function(n) {
    sort(clusterID[[n]])
  })
  weights <- rep(NA, n)

  for (i in 1:n) {
    K <- length(unique(result[[i]]$class))
    weights[i] <- mean(unlist(lapply(1:K, FUN = function(k) {
      mean(CCmat.sort[[i]][clusterID.sort[[i]] == k, clusterID.sort[[i]] == k])
    }))) /
      mean(unlist(lapply(1:K, FUN = function(k) {
        mean(CCmat.sort[[i]][clusterID.sort[[i]] == k, clusterID.sort[[i]] != k])
      })))
  }
  weights <- weights / sum(weights)
  return(weights)
}
