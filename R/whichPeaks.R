#' Find local maxima and local minima of data sequence
#'
#' @param x numerical vector contains local maxima (minima)
#' @param partial logical value indicating if the two endpoints will be considered
#' @param decreasing logical value indicating whether to find local minima
#'
#' @return a vector of locations of local maxima or minima
#' @export
#'
#' @examples
#' a = 100:1
#' which.peaks(a*sin(a/3))
which.peaks = function(x, partial=FALSE, decreasing=FALSE){
  if (decreasing){
    if (partial){
      which(diff(c(FALSE,diff(x)>0,TRUE))>0)
    }else {
      which(diff(diff(x)>0)>0)+1
    }
  }else {
    if (partial){
      which(diff(c(TRUE,diff(x)>=0,FALSE))<0)
    }else {
      which(diff(diff(x)>=0)<0)+1
    }
  }
}
