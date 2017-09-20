mcmc.wish.icar <- function(Dobs, TL, obs.idx, df=1,
                           beta.start = rep(0, length(TL)),
                           beta.prior.mean = rep(0, length(TL)),
                           beta.prior.cov = diag(10, length(TL)),
                           tau.start = 0.1, tau.prior.var = 1,
                           theta.tune = diag(10^-4,length(TL)+1),
                           n.mcmc=100, adapt.max=10000, adapt.int=100,
                           print.iter=FALSE, output.trace.plot=FALSE){


    ## ## wrapper function to create covariance matrix of observations
    get.Psi=function(TL,beta,tau,cells=cells,K=K,L=L){
        ## get precision matrix of whole graph
        Q=get.Q(TL,beta)
        ## get precision matrix of observed nodes
        max.diag=max(diag(Q))
        Q=Q/max.diag
        Phi=get.Phi(Q,cells)
        ## get covariance matrix of observations
        Sigma.nodes=ginv(as.matrix(Phi))
        Sigma.nodes=Sigma.nodes/max.diag
        Psi=K%*%Sigma.nodes%*%t(K)+tau*diag(nrow(K))
        Psi
    }

    ## Preliminaries----------------------------------------
    ## make L and W=L'(-Dobs)L ~ Wishart(df,L'(2Sigma)L)
    n=nrow(Dobs)
    L=diag(n-1)
    L=cbind(L,-1)
    L=t(L)
    W=t(L)%*%(-Dobs)%*%L

    ## number of params
    p <- length(TL)
    ## number of observations
    n.obs <- length(obs.idx)
    ## unique indexes of observations
    cells.idx <- unique(obs.idx)
    n.cells <- length(cells.idx)
    ## matrix for creating psi covariance matrix (above equation 11)
    K <- matrix(0, nrow = n.obs, ncol = n.cells)
    for (i in 1:n.obs){
        K[i, which(cells.idx == obs.idx[i])] <- 1
    }

    ## starting values
    beta=beta.start
    tau=tau.start
    logtau=log(tau)

    theta=c(logtau,beta)

    ## Starting value for ll
    Psi=get.Psi(TL,beta,tau,cells.idx,K,L)
    ll=dGenWish(Dobs,Psi,df,log=TRUE)
    ## MCMC ---------------------------------------------------
    beta.save=matrix(NA,n.mcmc,p)
    tau.save=rep(NA,n.mcmc)
    ll.save=rep(NA,n.mcmc)
    logprior.save=rep(NA,n.mcmc)
    dic.sum=0
    accept=0

    for(iter in 1:n.mcmc){
        if(print.iter)    cat(iter," ",beta,tau,"\n")
        if(iter%%100==0) cat (iter," ",beta,tau,"\n")
        if(iter%%100==0 & output.trace.plot){
            pdf("traceOut.pdf",width=12,height=5)
            matplot(beta.save[1:(iter),],type="l")
            legend("topleft",legend=1:length(beta),lwd=1,col=1:length(beta),bg="white")
            abline(h=0,col="yellow",lwd=2)
            dev.off()
        }

        ## propose
        pd=0
        while(pd==0){
            theta.star=rmvnorm(1,theta,theta.tune)
            beta.star=theta.star[-1]
            tau.star=exp(theta.star[1])

            ## get ll for gen wish
            Psi.star=try(get.Psi(TL,beta.star,tau.star,cells.idx,K,L),silent=TRUE)
            ll.star=try(dGenWish(Dobs,Psi.star,df,log=TRUE),silent=TRUE)
            if(is.numeric(ll.star)){
                pd=1
            }
        }

        ## MH step
        mh1=ll.star+dmvnorm(beta.star,beta.prior.mean,beta.prior.cov,log=TRUE)+dnorm(tau.star,0,sd=sqrt(tau.prior.var),log=TRUE)+log(tau.star) ## log(tau.star) is from jacobian for change of variables
        mh2=ll+dmvnorm(beta,beta.prior.mean,beta.prior.cov,log=TRUE)+dnorm(tau,0,sd=sqrt(tau.prior.var),log=TRUE)+log(tau) ## log(tau) is from jacobian for change of variables

        if(runif(1)<exp(as.numeric(mh1-mh2))){
            theta=theta.star
            tau=tau.star
            beta=beta.star
            ll=ll.star
            accept=accept+1
        }

        ## DIC update: (at end of chain, dic.sum=Dbar)
        ##   done by taking second half of chain
        if(iter>(n.mcmc/2)){
            dic.sum=dic.sum+1/(n.mcmc/2)*(-2*ll)
        }

        ## Save values---------------------------------------------

        beta.save[iter,] <- beta
        tau.save[iter] <- tau
        ll.save[iter] <- ll
        logprior.save[iter] <- dmvnorm(beta,beta.prior.mean,beta.prior.cov,log=TRUE)+dnorm(tau,0,sd=sqrt(tau.prior.var),log=TRUE)+log(tau)

        ## Adaptive MCMC for beta ---------------------------------
        if (iter < adapt.max + 1 & iter / adapt.int == round(iter/adapt.int)){
            theta.tune <- 2.4^2 /length(theta) *
                var(cbind(log(tau.save[1:iter]),beta.save[1:iter,]))+
                diag(10^-10,length(theta))
        }
    }

    ## DIC --------------------------------------------------------
    ##browser()
    Dbar=dic.sum
    if(ncol(beta.save)==1) {
        beta.hat=mean(beta.save[-c(1:n.mcmc/2),])
    } else {
        beta.hat=apply(beta.save[-c(1:n.mcmc/2),],2,mean)
    }

    tau.hat=mean(tau.save[-c(1:n.mcmc/2)])

    Psi=get.Psi(TL,beta.hat,tau.hat,cells.idx,K,L)
    ll.hat=dGenWish(Dobs,Psi,df,log=TRUE)
    
    Dhat=-2*ll.hat

    DIC=2*Dbar-Dhat
    DIC2=Dhat + var(-2*ll.save[-c(1:n.mcmc/2)])

    list(beta=beta.save,tau=tau.save,n.mcmc=n.mcmc,ll=ll.save,logprior=logprior.save,
         accept=accept,DIC=DIC,DIC2=DIC2,Dbar=Dbar,Dhat=Dhat)
}
