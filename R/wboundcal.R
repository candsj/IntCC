#' Title
#'
#' @param data
#' @param K
#'
#' @return
#' @export
#'
#' @examples
wboundcal <- function(data = NULL, K = k) {
  # k-means
  res <- stats::kmeans(data, K, iter.max = 1000, nstart = 20)
  # between sum of squares
  average_bs <- res$betweenss / length(unique(res$cluster))
  # regression predict
}
