#' Title
#'
#' @param clustersize A number or a vector denotes the number of samples for each cluster.
#' @param info_features Number of informative features.
#' @param total_features Number of total features.
#' @return simulation datasets
#' @export
#'
#' @examples
#' library(weightedCC)
#' simulation=data_generation(20,5,30)

data_generation <- function(clustersize, info_features, total_features) {
  p <- total_features
  if (length(clustersize) == 1) {
    n1 <- n2 <- n3 <- clustersize
    n <- 3 * clustersize
  } else {
    n1 <- clustersize[1]
    n2 <- clustersize[2]
    n3 <- clustersize[3]
    n <- n1 + n2 + n3
  }

  if (length(info_features) == 1) {
    p1 <- p2 <- p3 <- info_features
  } else {
    p1 <- info_features[1]
    p2 <- info_features[2]
    p3 <- info_features[3]
  }

  ################# normal distribution ####################
  c1 <- matrix(stats::rnorm(n1 * p1, mean = 2.2, sd = 1.7), ncol = p1, nrow = n1)
  c2 <- matrix(stats::rnorm(n2 * p2, mean = 1.6, sd = 1.7), ncol = p2, nrow = n2)
  c3 <- matrix(stats::rnorm(n3 * p3, mean = 1, sd = 1.7), ncol = p3, nrow = n3)

  normData <- matrix(stats::rnorm(n * p, mean = 0, sd = 1.7), nrow = n, ncol = p)
  normData[1:n1, 1:p1] <- c1
  normData[(n1 + 1):(n1 + n2), (p1 + 1):(p1 + p2)] <- c2
  normData[(n1 + n2 + 1):n, (p1 + p2 + 1):(p1 + p2 + p3)] <- c3

  ######## binomial data #########################
  c1 <- matrix(stats::rbinom(n1 * p1, size = 1, prob = 0.54), ncol = p1, nrow = n1)
  c2 <- matrix(stats::rbinom(n2 * p2, size = 1, prob = 0.42), ncol = p2, nrow = n2)
  c3 <- matrix(stats::rbinom(n3 * p3, size = 1, prob = 0.3), ncol = p3, nrow = n3)

  binomData <- matrix(stats::rbinom(n * p, size = 1, prob = 0.1), nrow = n, ncol = p)
  binomData[1:n1, 1:p1] <- c1
  binomData[(n1 + 1):(n1 + n2), (p1 + 1):(p1 + p2)] <- c2
  binomData[(n1 + n2 + 1):n, (p1 + p2 + 1):(p1 + p2 + p3)] <- c3

  ######### poisson data ############################
  c1 <- matrix(stats::rpois(n1 * p1, lambda = 1.9), ncol = p1, nrow = n1)
  c2 <- matrix(stats::rpois(n2 * p2, lambda = 1.5), ncol = p2, nrow = n2)
  c3 <- matrix(stats::rpois(n3 * p3, lambda = 1), ncol = p3, nrow = n3)

  poisData <- matrix(stats::rpois(n * p, lambda = 0.6), nrow = n, ncol = p)
  poisData[1:n1, 1:p1] <- c1
  poisData[(n1 + 1):(n1 + n2), (p1 + 1):(p1 + p2)] <- c2
  poisData[(n1 + n2 + 1):n, (p1 + p2 + 1):(p1 + p2 + p3)] <- c3


  ################## multinomial distribution ##############
  c1 <- matrix(sample(1:3, size = n1 * p1, replace = TRUE, prob = c(0.8, 0.15, 0.05)), ncol = p1, nrow = n1)
  c2 <- matrix(sample(1:3, size = n2 * p2, replace = TRUE, prob = c(0.3, 0.4, 0.3)), ncol = p2, nrow = n2)
  c3 <- matrix(sample(1:3, size = n3 * p3, replace = TRUE, prob = c(0.05, 0.15, 0.8)), ncol = p3, nrow = n3)

  multData <- matrix(sample(1:3, size = n * p, replace = TRUE, prob = c(0.33, 0.33, 0.33)), nrow = n, ncol = p)
  multData[1:n1, 1:p1] <- c1
  multData[(n1 + 1):(n1 + n2), (p1 + 1):(p1 + p2)] <- c2
  multData[(n1 + n2 + 1):n, (p1 + p2 + 1):(p1 + p2 + p3)] <- c3


  # biID = seq(1,p,2)
  # multData[,biID] = binomData[,biID]

  ### some variables only have one category; need to remove these variables
  ncat <- function(x) {
    length(unique(x))
  }
  len <- apply(binomData, 2, ncat)
  unicat <- which(len == 1)
  if (length(unicat) > 0) {
    binomData <- binomData[, -unicat]
  }

  data_for_klic <- list()
  data_for_klic[[1]] <- normData
  data_for_klic[[2]] <- binomData
  data_for_klic[[3]] <- poisData
  data_for_klic[[4]] <- multData
  return(data_for_klic)
}
