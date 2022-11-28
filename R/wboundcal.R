#' @keywords internal


wboundcal <- function(data = NULL, K = 2, pFeature=1) {
  # permute
  res <- KMeansSparseCluster.permute.ver(data,K,nperms=100,wbounds=exp(seq(log(1.2), log(sqrt(dim(data)[2])*1.5), len=50)))$bestw
  # regression predict
  if(pFeature!=1){
    res=res*pFeature^0.4454441
  }
  return(res)
}
