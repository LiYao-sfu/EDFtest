#' Cramer Von Mises Statistics for the Poisson Distribution
#'
#' @param x
#' @param eps how much prob you are willing to omit
#'
#' @return
#' @export
#'
#' @examples
CvM.poisson = function(x,eps=10^(-9)){
    n=length(x)
    s=mean(x)
    obs=tabulate(x+1)
    #    print(obs)
    m=max(x)
    p=dpois(0:m,s)
    tw=p
    #
    # The Cramer von Mises statistic uses weights corresponding to dF
    #  in the continuous case.  There are at least 2 reasonable possible
    #  formulas for the weights.  We use our original choice here
    #
    # The formulas below are for weights which would give the same statistic if the #  cells were reversed in order.
    #
    #    p1=dpois(1:(m+1),s)
    #
    #    tw=(p+p1)/2
    #
    e = n*p
    z=cumsum(obs-e)
    #
    # The statistic is technically an infinite sum which we truncate
    #  in the code below to miss only eps of the tail probability.
    #  If the more elaborate weights above are used the code below
    #  needs to change, too.
    #
    m2=max(qpois(1-eps,s),100)
    p2=dpois((m+1):m2,s)
    q2=ppois(m:(m2-1),s,lower.tail=F)
    tail=n*sum(q2^2*p2)
    sum(z^2*tw)/n+tail
}
