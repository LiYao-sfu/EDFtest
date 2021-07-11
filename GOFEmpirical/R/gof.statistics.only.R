#' gof.statistics.only
#'
#' @param pit
#' @param a2
#' @param w2
#' @param u2
#'
#' @return
#' @export
#'
#' @examples
gof.statistics.only=function(pit,a2=T,w2=T,u2=T){
  #
  # This function takes in a vector of probability integral transforms
  #  and computes whichever of the three statistics is asked for.
  # By default A^2,W^2, and U^2 are computed
  #
  if(any(pit<0))stop("Negative Probability Integral Transforms Not Allowed")
  if(any(pit>1))stop("Probability Integral Transforms More than 1 Not Allowed")
  p=sort(pit)
  out=list()
  if(a2){
    ind = ((p>=1)||(p <= 0))
    if(any(ind)){
      if(any(p >= 1)){ warning("A2 not computed; some pit is 1")}
      if(any(p <= 0)){ warning("A2 not computed; some pit is 0")}
    }
    else out$A2 = AD(p)
  }
  if(w2){
    out$W2=CvM(p)
  }
  if(u2){
    out$U2=usq(p)
  }
  out
}
