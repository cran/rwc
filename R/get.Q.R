get.Q <-
function(TL,beta){
    p=length(TL)
    M=nrow(TL[[1]])
    A=Matrix(0,nrow=M,ncol=M,sparse=T)
    for(i in 1:p){
        A=A+beta[i]*TL[[i]]
    }
    A@x=exp(A@x)
    m=rowSums(A)
    Q=Diagonal(M,m)-A
    Q=1/2*(Q+t(Q))
    rm(A)
    Q
}
