\name{mcmc.wish.icar}
\alias{mcmc.wish.icar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Markov chain Monte Carlo sampler for Generalized Wishart distance
  matrix data arising from an ICAR spatial model.
%%  ~~function to do ... ~~
}
\description{Constructs and runs an MCMC sampler to estimate resistance
  parameters of different landscape features.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mcmc.wish.icar(Dobs, TL, obs.idx, df=1,
               beta.start = rep(0, length(TL)),
               beta.prior.mean = rep(0, length(TL)),
               beta.prior.cov = diag(10, length(TL)),
               tau.start = 0.1, tau.prior.var = 1,
               theta.tune = diag(10^-4,length(TL)+1),
               n.mcmc=100, adapt.max=10000, adapt.int=100,
               print.iter=FALSE, output.trace.plot=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Dobs}{A square symmetric matrix of observed pairwise distances.
    For example, a genetic distance matrix.}
  \item{TL}{A list of transition matrices for different covariate raster
    layers, created by get.TL}
  \item{obs.idx}{A vector of unique indices of observed nodes in the
    graph defined by the raster grid. }
  \item{df}{An integer specifying the degrees of freedom of Dobs.}
  \item{beta.start}{Vector of initial values for conductance parameters
    beta.  Default is a vector of zeros.}
  \item{beta.prior.mean}{Vector of prior mean values for conductance parameters
     beta.  Default is a vector of zeros.}
  \item{beta.prior.cov}{Matrix of the prior covariance matrix for conductance parameters
     beta.  Default is a diagonal matrix with diagonal entries equal to
     10.}
   \item{tau.start}{Numeric starting value for the nugget variance tau.
   Default is 0.1}
   \item{tau.prior.var}{Variance of the half-normal prior for tau.
     Default is 1.}
   \item{theta.tune}{Covariance matrix for the random walk MH sampler
     for all parameters.  Default is a diagonal matrix with variance
     10^-4.}
   \item{n.mcmc}{Integer number of iterations of the MCMC sampler to
     run.}
   \item{adapt.max}{Integer number (or INF) specifying the last iteration
     at which the covariance matrix of the proposal distribution will be
     adapted.  Default is 10^5.}
   \item{adapt.int}{Interval at which the covariance matrix of the proposal
     distribution is adapted.  Default is every 100 iterations.}
   \item{print.iter}{Logical.  If TRUE, then the current state of the
     system will be printed to the console every 100 iterations.}
   \item{output.trace.plot}{Logical.  If TRUE, then the trace plots of the
     sampler will be saved to "traceOut.pdf" every 100 iterations.}
   
}
\details{
  Runs an MCMC sampler to draw samples from the posterior distribution
  of model parameters (beta,tau) from the following model for an
  observed squared distance matrix Dobs:

  -Dobs ~ GenWish(2*Sigma,df)

  Sigma = K(Psi)K'+tau*I

  where Psi is the covariance matrix of the observed nodes of a graph
  described by the transition list TL.  That is, the total graph has
  Laplacian (precision) matrix Q, with off-diagonal entries of Q given
by

Q_ij = exp( b0 + b1 x_1ij + ... + bp x_pij ), where beta=(b1,b2,...,bp) is
  a vector of log-conductance values of the covariates {x_kij}.  Each
  x_kij is equal to (x_ki+x_kj)/2.

  The prior on beta is N(beta.prior.mean,beta.prior.cov), and the prior
  on tau is tau~Half_Normal(0,tau.prior.var).
}
\value{ A list containing output from the MCMC sampler.
  \item{beta}{Posterior samples for conductance parameters beta.}
}
\references{
  Hanks and Hooten 2013.  Circuit theory and model-based inference for
  landscape connectivity.  Journal of the American Statistical
  Association.  108(501), 22-33.
%% ~put references to the literature/web site here ~
}
\author{Ephraim M. Hanks
%%  ~~who you are~~
}
%% \note{
%% %%  ~~further notes~~
%% }

%% %% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
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

##
## fitting using MCMC
## 

fit=mcmc.wish.icar(D.rand.ibr,TL,samp.locs,df=20,output.trace.plot=TRUE,
                   adapt.int=100,adapt.max=100000,n.mcmc=10000)

str(fit)

}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
