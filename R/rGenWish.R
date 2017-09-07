rGenWish <-
function(Sigma,df){
    n=nrow(Sigma)
    Y=t(rmvnorm(df,sigma=Sigma))
    S=Y%*%t(Y)
    onevec=matrix(1,n,1)
    D.gen=-2*S+diag(S)%*%t(onevec)+onevec%*%t(diag(S))
    D.gen
    }
