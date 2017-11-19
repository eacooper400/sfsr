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
    thetaH = (2*thetaL) - thetaPi
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
    num1 = (((36*(n^2) * ((2*n)+1) * ((1/n)+sum(1/((seq(1,n))^2))))*2)*thetaSQ)
    num2=(2*(n^3)*thetaSQ)*116
    num3 = ((((n^2)*9) + (2*n) + 3)*2)*thetaSQ
    top=(num1-(num2 + num3))
    v=top/((9*n) * ((n-1)^2))
    return(v)
}
