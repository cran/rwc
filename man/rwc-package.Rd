\name{rwc-package}
\alias{rwc-package}
\alias{rwc}
\docType{package}
\title{
Random Walk Covariance Models
}
\description{
Code to facilitate simulation and inference of connectivity defined by random walks.
}
\details{
This package contains code to simulate (rGenWish) and evaluate the
likelihood of (dGenWish) distance matrices from the Generalized Wishart
distribution.  It also contains helper functions to create and manage
spatial covariance models created from landscape grids with resistance
or conductance defined by landscape features.
}
\author{
Ephraim M. Hanks

Maintainer: Ephraim M. Hanks
}
\references{
McCullagh 2009.  Marginal likelihood for distance matrices.
Statistica Sinica 19: 631-649.

Hanks and Hooten 2013.  Circuit theory and model-based inference for
  landscape connectivity.  Journal of the American Statistical
  Association.  108(501), 22-33.

  Hanks 2017. Modeling spatial covariance using the limiting
  distribution of spatio-temporal random walks.  Journal of the American
  Statistical Association.  
}%% \seealso{
%% ~~ Optional links to other man pages, e.g. ~~
%% ~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%% }
\examples{

\dontrun{
## The following code conducts a simulation example by
## first simulating a heterogeneous landscape, then
## simulating pairwise distance data on the landscape
## and finally making inference on model parameters.

library(rwc)
library(MASS)

## source("GenWishFunctions_20170901.r")

##
## specify 2-d raster
##


ras=raster(nrow=30,ncol=30)
extent(ras) <- c(0,30,0,30)
values(ras) <- 1
plot(ras,asp=1)
ras

int=ras
cov.ras=ras



## simulate "resistance" raster
TL.int=get.TL(int)
Q.int=get.Q(TL.int,1)
set.seed(1248)
## values(cov.ras) <- as.numeric(rnorm.Q(Q.int%*%Q.int,zero.constraint=TRUE))
values(cov.ras) <- as.numeric(rnorm.Q(Q.int,zero.constraint=TRUE))
plot(cov.ras)


stack=stack(int,cov.ras)
plot(stack)
TL=get.TL(stack)


## Create "barrier" raster which has no effect, to test model selection

barrier=int
values(barrier) <- 0
barrier[,10:11] <- 1

plot(barrier)

TL.all=get.TL(stack(int,cov.ras,barrier))


##
## sampling locations
##

nsamps=150
set.seed(4567)
samp.xy=cbind(runif(nsamps,0,30),runif(nsamps,0,30))
## samp.xy=samp.xy[rep(1:nsamps,times=5),]
samp.locs=cellFromXY(int,samp.xy)
samp.cells=unique(samp.locs)
nsamps=nrow(samp.xy)

plot(cov.ras)
points(samp.xy)

K=matrix(0,nrow=nsamps,ncol=length(samp.cells))
for(i in 1:nsamps){
    K[i,which(samp.cells==samp.locs[i])]=1
}
image(K)

##
## beta values
##


betas=c(-2,-1)
tau=.01


Q=get.Q(TL,betas)
Phi=get.Phi(Q,samp.cells)



## simulate from ibr model
D.rand.ibr=rGenWish(df=20,Sigma=K\%*\%ginv(as.matrix(Phi))\%*\%t(K) + diag(nsamps)*tau)
image(D.rand.ibr)


## crude plot of geographic distance vs genetic distance

plot(as.numeric(as.matrix(dist(xyFromCell(ras,samp.locs)))),as.numeric(D.rand.ibr))

###################################
##
##
## fitting using optim
## 
##


nll.gen.wish.icar <- function(theta,D,df,TL,obs.idx){
    ## get K
    cells.idx=unique(obs.idx)
    n.cells=length(cells.idx)
    n.obs=length(obs.idx)
    K <- matrix(0, nrow = n.obs, ncol = n.cells)
    for (i in 1:n.obs){
        K[i, which(cells.idx == obs.idx[i])] <- 1
    }
    ## get precision matrix of whole graph
    tau=exp(theta[1])
    beta=theta[-1]
    Q=get.Q(TL,beta)
    ## get precision matrix of observed nodes
    max.diag=max(diag(Q))
    Q=Q/max.diag
    Phi=get.Phi(Q,cells.idx)
    ## get covariance matrix of observations
    Sigma.nodes=ginv(as.matrix(Phi))
    Sigma.nodes=Sigma.nodes/max.diag
    Psi=K\%*\%Sigma.nodes\%*\%t(K)+tau*diag(nrow(K))
    ## get nll
    nll=-dGenWish(D,Psi,df,log=TRUE)
    nll
}


fit=optim(c(log(tau),betas),nll.gen.wish.icar,D=D.rand.ibr,df=20,TL=TL,
    obs.idx=samp.locs,control=list(trace=10),hessian=TRUE)

fit.all=optim(c(log(tau),betas,0),nll.gen.wish.icar,D=D.rand.ibr,df=20,TL=TL.all,
    obs.idx=samp.locs,control=list(trace=10),hessian=FALSE)

fit.ibd=optim(c(log(tau),0),nll.gen.wish.icar,D=D.rand.ibr,df=20,TL=TL.int,
    obs.idx=samp.locs,control=list(trace=10),hessian=FALSE)


## model selection using AIC

aic.ibr=2*fit$value + 2*length(fit$par)
aic.all=2*fit.all$value + 2*length(fit.all$par)
aic.ibd=2*fit.ibd$value + 2*length(fit.ibd$par)

aic.ibr
aic.all
aic.ibd

## se's for best fit

str(fit)
## get standard errors from optim
H=fit$hessian
S=solve(H)
se=sqrt(diag(S))
se

## CI's for best fit

CImat=matrix(NA,3,4)
rownames(CImat) <- c("log(tau)","intercept","resistance")
colnames(CImat) <- c("truth","estimate","lower95CI","upper95CI")
CImat[,1]=c(log(tau),betas)
CImat[,2]=fit$par
CImat[,3]=fit$par-1.96*se
CImat[,4]=fit$par+1.96*se

CImat
}
}
