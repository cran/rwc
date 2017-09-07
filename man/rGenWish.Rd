\name{rGenWish}
\alias{rGenWish}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate realizations from the Generalized Wishart distribution
%%  ~~function to do ... ~~
}
\description{Simulates Wishart random variables, then computes the
  induced distance of the simulated Wishart random variables.  The
  result is a random matrix distributed as a Generalized Wishart random variable.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rGenWish(Sigma, df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Sigma}{The covariance parameter of the Generalized Wishart.
%%     ~~Describe \code{Sigma} here~~
}
  \item{df}{An integer specifying the degrees of freedom.
%%     ~~Describe \code{nu} here~~
}
}
%% \details{
%% %%  ~~ If necessary, more details than the description above ~~
%% }
\value{A matrix of dimension equal to the dimension of Sigma.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  McCullagh 2009.  Marginal likelihood for distance matrices.
  Statistica Sinica 19: 631-649.

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
ras=raster(nrow=30,ncol=30)
extent(ras) <- c(0,30,0,30)
values(ras) <- 1
int=ras
cov.ras=ras

## get precision matrix of entire graph
B.int=get.TL(int)
Q.int=get.Q(B.int,1)

## get precision at a few nodes
Phi=get.Phi(Q.int,obs.idx=1:20)
S=ginv(as.matrix(Phi))

## simulate distance matrix
Dsim=rGenWish(df=20,Sigma=S)
image(Dsim)

## calculate log-likelihood
ll=dGenWish(Dsim,S,df=20,log=TRUE)
ll

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
