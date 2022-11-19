#' @keywords internal
UpdateWs.ver <- function(x, Cs, l1bound){
wcss.perfeature <- utils::getFromNamespace("GetWCSS","sparcl")(x, Cs)$wcss.perfeature
tss.perfeature <- utils::getFromNamespace("GetWCSS","sparcl")(x, rep(1, nrow(x)))$wcss.perfeature
lam <- utils::getFromNamespace("BinarySearch","sparcl")(-wcss.perfeature+tss.perfeature, l1bound)
ws.unscaled <- utils::getFromNamespace("soft","sparcl")(-wcss.perfeature+tss.perfeature,lam)
ws.unscaled[ws.unscaled<0] <- 0 ## added this line to ensure that negative estimates are set to zero
return(ws.unscaled/utils::getFromNamespace("l2n","sparcl")(ws.unscaled))
}

