\name{dGenWish}
\alias{dGenWish}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Density of the (singular) Generalized Wishart distribution
%%  ~~function to do ... ~~
}
\description{Computes the density of the (possibly singular) Generalized
  Wishart distribution with null-space equal to the space spanned by the
  "one" vector.  This corresponds to the case considered by McCullagh
  (2009) and Hanks and Hooten (2013).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dGenWish(Dobs, Sigma, df,log=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Dobs}{An observed squared-distance matrix.
%%     ~~Describe \code{Dobs} here~~
}
  \item{Sigma}{The covariance parameter of the Generalized Wishart.
%%     ~~Describe \code{Sigma} here~~
}
  \item{df}{An integer specifying the degrees of freedom.
%%     ~~Describe \code{nu} here~~
}
  \item{log}{Logical.  If True, then the log-likelihood is computed.
%%     ~~Describe \code{nu} here~~
}
}
\details{Following McCullagh (2009), the likelihood can be computed by considering any contrast matrix L of full rank, and with n-1 rows and n columns,
  where n is the number of columns of 'Dobs'.  If

  Dobs ~ GenWish(Sigma,df,1)

  is distributed as a generalized Wishart distribution with kernel (null
  space) equal to the one vector, and df degrees of freedom, then the
  likelihood can be computed by computing the likelihood of

  L(-Dobs)L' ~ Wishart(L(2*Sigma)L',df)

  Additionally, following Srivastava (2003), this likelihood holds (up
  to a proportionality constant) in the singular case where df<n.

  Following this formulation, the log-likelihood computed here (up to an
  additive constant) is

  -df/2*log|L(2*Sigma)L'| -1/2*tr( (L(2*Sigma)L')^-1 L(-D)L' )
  
%%  ~~ If necessary, more details than the description above ~~
}
\value{A numeric likelihood or log-likelihood
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  McCullagh 2009.  Marginal likelihood for distance matrices.
  Statistica Sinica 19: 631-649.

  Srivastava 2003.  Singular Wishart and multivariate beta
  distributions.  The Annals of Statistics.  31(5), 1537-1560.

  Hanks and Hooten 2013.  Circuit theory and model-based inference for
  landscape connectivity.  Journal of the American Statistical
  Association.  108(501), 22-33.
%% ~put references to the literature/web site here ~
}
\author{
  Ephraim M. Hanks
%%  ~~who you are~~
}
%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

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
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
