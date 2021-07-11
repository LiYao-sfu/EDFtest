#' MLE for Logistic Distribution
#'
#'
#'
#' @param x
#' @param eps
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
#'
estimate.logistic <- function(x,eps=1e-7,verbose=FALSE){
    hessianscore.logistic=function(x,theta){
      #
      # This is the hessian and score function for the two parameter logistic distribution
      #
      scale=theta[2]
      location=theta[1]
      y=(x-location)/scale
      e=exp(-y)
      r=(1-e)/(1+e)
      ell = sum(-log(scale) + y-2*log(1+exp(y)))
      s.scale= y*r/scale -1/scale
      s.location= r/scale
      score = c(sum(s.location),sum(s.scale))
      h.den=(scale*(1+e))^2
      r.mm = -2*sum(e/h.den)
      r.ms = sum((e^2-2*y*e-1)/h.den)
      r.ss= sum(((2*y+1)*e^2-2*(y^2-1)*e-(2*y-1))/h.den)
      Hessian = matrix(c(r.mm,r.ms,r.ms,r.ss),nrow=2)
      ev=eigen(-Hessian)$values
      if(any(ev <= 0)) cat("Warning: Hessian contains some eigenvalue of the wrong sign\n")
      list(loglikelihood = ell,score=score,Hessian=Hessian)
    } # for Newton-Raphson

    m1 <- mean(x)
    m2 <- var(x)
    iterates=0
    binit = sqrt(m2*3/pi^2)
    ainit=m1
    thetaold=c(ainit,binit)
    ders = hessianscore.logistic(x,theta=thetaold)
    ell.old=ders$loglikelihood
    if(verbose){
      cat("Initial Estimates and Likelihood\n")
      cat(thetaold)
      cat(ders$loglikelihood)
      cat("\n Score\n")
      cat(ders$score)
    }
    step = solve(ders$Hessian,ders$score)

    thetanew = thetaold-step
    while(thetanew[2]<0){
      if(verbose) cat("Step too far so \n New theta ",thetanew," step size ",step,"\n\n")
      step[2]=step[2]/2
      thetanew = thetaold-step
    }
    ders = hessianscore.logistic(x,theta=thetanew)
    iter=0
    while( ders$loglikelihood < ell.old){
      iter=iter+1
      step=step/2
      thetanew = thetaold-step
      if(verbose){
        cat("Likelihood didn't increase \n")
        cat(" was ",ell.old," now ",ders$loglikelihood,"\n")
        cat(" so \n New theta ",thetanew," step size ",step,"\n\n")
      }
      ders = hessianscore.logistic(x,theta=thetanew)
      if(verbose){
        cat(ders$score)
        cat("\n Iteration ",iter," log likelihood ",ders$loglikelihood,"\n")
      }
    }
    delta = sqrt(sum((thetanew-thetaold)^2))
    if(verbose){
      print(cbind(thetaold, thetanew, ders$score,abs(thetaold-thetanew),ders$loglikelihood))
    }
    while ( abs(delta) > eps){
      if(verbose){
        cat("Looping\n")
      }
      thetaold=thetanew
      ders = hessianscore.logistic(x,theta=thetanew)
      step =solve(ders$Hessian,ders$score)
      thetanew = thetaold-step
      delta = sqrt(sum((thetanew-thetaold)^2))
      if(verbose){
        print(cbind(thetaold, thetanew, ders$score, abs(thetaold-thetanew),ders$loglikelihood))
      }
      iterates=iterates+1
      if(iterates>50)break()
    }
    thetanew
  }
