dist.from.cov <-
function(Sigma){
    ## function to compute the resistance distance matrix
    ## (or equivalently the variogram matrix)
    ## from a covariance matrix
    ## See Hanks and Hooten (2013)
    d=diag(Sigma)
    n=nrow(Sigma)
    one=matrix(1,n,1)
    R=d%*%t(one) + one%*%t(d) -2*Sigma
    R
}
