#' Compute FDR threshold based on Benjamini-Hochberg (BH) algorithm
#'
#' @param p a vector of p-values
#' @param q False Discovery Rate level
#'
#' @return p-value threshold based on independence or positive dependence
#' @export
#' @examples
#' fdrBH(seq(0.01,0.1,0.01),q=0.1)
#'
fdrBH = function(p, q){
  p <- sort(p)
  V <- length(p)
  I <- 1:V
  cVID <- 1
  # cVN = sum(1/I)
  pID <- ifelse(all(p>I/V*q/cVID),0,max(p[p<=I/V*q/cVID]))
  # pN = max(p[p<=I/V*q/cVN])
  return(pID)
}
