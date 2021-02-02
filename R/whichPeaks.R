#' Find local maxima or minima
#'
#' @param x a vector with maxima or minima
#' @param partial if 'true' endpoints will be considered
#' @param decreasing find local minima if 'true', ortherwise local maxima
#'
#' @return a vector of positions of local maxima or minima
#' @export
#'
#' @examples
#' a = 100:1
#' which.peaks(a*sin(a/3))
which.peaks = function(x,partial=FALSE,decreasing=FALSE){
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
