cov.from.dist <-
function(R){
    ## function to compute a covariance matrix from
    ## a resistance distance (or variogram) matrix
    ## using the centering method of Gower 1996
    n=nrow(R)
    one=matrix(1,n,1)
    sum.all=sum(R)/n^2
    row.sum.matrix= ( apply(R,1,mean) )%*%t(one)
    col.sum.matrix= one%*%t( apply(R,2,mean) )
    Sigma = 1/2* (-R + row.sum.matrix + col.sum.matrix - sum.all)
    Sigma
}
