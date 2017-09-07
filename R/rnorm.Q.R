rnorm.Q <-
function(Q,mu=rep(0,nrow(Q)),X=Matrix(1, nrow = nrow(Q), ncol = 1),zero.constraint=FALSE,canon=FALSE,diag.adjust=.Machine$double.eps*10){

    n=nrow(Q)
    
    diag(Q)=diag(Q)*(1+diag.adjust)
    L=t(chol(Q))
    
    if(canon==FALSE){
        z=rnorm(n)
        v=solve(t(L),z)
        eta.star=mu+v
    }
    if(canon==TRUE){
        w=solve(L,mu)
        q=solve(t(L),w)
        z=rnorm(length(mu))
        v=solve(t(L),z)
        eta.star=q+v
    }   

   qsolve=function(L,B){
       ## solve QX=B where Q=LL^T
       V=solve(L,B)
       X=solve(t(L),V)
       X
   }
    
    if(zero.constraint){
        if(is.na(X)[1]){
            X=Matrix(1,nrow=n,ncol=1)
        }
        V=qsolve(L,X)
        W=t(X)%*%V
        W=forceSymmetric(W)
        ## L.w.try=try(t(chol(W)))
        ## sink()
        ## if(length(dim(L.w.try))==0){
        ##     diag(W)=diag(W)+10^-12
        ##     L.w.try=t(chol(W))
        ## }else{
        ##     L.w=L.w.try
        ## }
        ##L.w=t(chol(W))
        ##W=nearPD(W)
        L.w=t(chol(W))
        U=qsolve(L.w,t(V))
        U=solve(W)%*%t(V)
        c=t(X)%*%eta.star
        eta=eta.star-t(U)%*%c
    }else{
        eta=eta.star
    }

    eta
}
