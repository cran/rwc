dGenWish <-
function(Dobs,Sigma, df,log=FALSE){
    ## calculate generalized wishart log-likelihood
    ## -D ~ GW(1,df=nu,Sigma)
    ## Sigma = covariance matrix, nu=degrees of freedom
    ##
    ## ll= df/2*log|Sigma.inv*A| + 1/4*tr(Sigma.inv*A*D)
    ##
    L=Diagonal(nrow(Sigma)-1,1)
    L=cbind(L,-1)
    LDL= L%*%(-Dobs)%*%t(L)
    W= L%*%(2*Sigma)%*%t(L)
    EE=eigen(W,symmetric=TRUE)
    W.inv=EE$vect%*%diag(1/EE$val)%*%t(EE$vect)
    WLDL=W.inv%*%LDL
    ll=-df/2*sum(log(EE$val))-sum(diag(WLDL))
    if(!log){
        ll=exp(ll)
    }
    ll
}
