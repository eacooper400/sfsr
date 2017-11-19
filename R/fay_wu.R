### Calculate Fay and Wu's H from an SFS
### Uses the normalized statistic from 
### Zheng et al. (2006) Genetics 174(3)

#' Calculate Normalize Fay and Wu's H
#'
#' @param sfs an SFS
#' @return a single numeric value for the estimate of H
#'
#' @export
fayWu_H <- function(sfs) {
    n=length(sfs)-1
    x=mask_corners(t(as.matrix(sfs)))
    s=sum(x)
    thetaPi = theta_pi(sfs)
    thetaW = theta_w(sfs)
    thetaL=((n-1)/n)*thetaW
    H=thetaPi - thetaL
    varH = variance.H(n, s, thetaW)
    H = H/(sqrt(varH))
    return(H)
}

#' Calculate a sub n, a component of the variance
#'
#' @param n an integer value for the number of chromosomes
#' @return a single numeric value for the estimate of a1
#'
#' @export
a1 <- function(n) {
    return(sum(1/seq_len(n-1)))
}

#' Calculate b sub n, a component of the variance
#'
#' @param n an integer value for the number of chromosomes
#' @return a single numeric value for the estimate of b1
#'
#' @export
b1 <- function(n) {
    return(sum(1/seq_len(n-1)^2))
}

#' Calculate the variance of H
#'
#' @param n an integer value for the number of chromosomes
#' @param s an integer value for the number of segregating sites
#' @param theta an estimate of theta (best to use thetaW)
#'
#' @return a single number value for the variance
#'
#' @export
variance.H <- function(n,s,theta) {
    thetaSQ = (s*(s-1))/(((a1(n))**2)+b1(n))
    f1 = ((n-2)/(6*(n-1))) * theta
    num1 = 18 * (n**2) * ((3*n)+2) * b1(n+1)
    num2 = (88 * (n**3)) + (9*(n**2)) - (13*n) + 6
    denom = (9*n) * ((n-1)**2)
    f2 = ((num1 - num2)/denom) * thetaSQ
    var = f1+f2
    return(var)
}

