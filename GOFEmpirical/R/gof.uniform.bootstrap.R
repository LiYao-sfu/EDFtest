#' gof.uniform.bootstrap
#'
#' @param x
#' @param stat
#' @param m
#'
#' @return
#' @export
#'
#' @examples
gof.uniform.bootstrap=function(x,stat=y,m=10000){
  Z=sort(x)
  n=length(Z)
  r=matrix(runif(n*m),nrow=m)
  r = t(apply(r,1,sort))
  print(dim(r))
  if(stat==1)
  {
    Wsq=CvM(Z)
    Wsq_r=apply(r,1,CvM)
    pval=sum(Wsq_r>=Wsq)/m
  }
  else if(stat==2)
  {
    Usq=usq(Z)
    Usq_r=apply(r,1,usq)
    pval=sum(Usq_r>=Usq)/m
  }
  else
  {
    Asq=AD(Z)
    Asq_r=apply(r,1,AD)
    pval=sum(Asq_r>=Asq)/m
  }
  return(pval)
}
