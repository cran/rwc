get.Phi <-
function(Q,obs.idx){
    if(length(obs.idx)<nrow(Q)){
        Q11=Q[obs.idx,obs.idx]
        Q22=Q[-obs.idx,-obs.idx]
        Q21=Q[-obs.idx,obs.idx]
        L=t(chol(Q22))
        V=solve(L,Q21)
        Phi=Q11-t(V)%*%V
    }else{
        Phi=Q
    }
}
