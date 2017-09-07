\name{cov.from.dist}
\alias{cov.from.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create covariance matrix from a distance matrix
%%  ~~function to do ... ~~
}
\description{ This computes a covariance matrix from a squared-distance
  matrix using the centering method of Gower (1996).  When the
  squared-distance matrix is a resistance distance matrix, or a
  variogram matrix from a spatial model, the resulting covariance matrix
  is the spatial covariance matrix corresponding to a random walk model
  for connectivity as in Hanks and Hooten (2013).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cov.from.dist(R)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R}{A negative semi-definite matrix of squared differences.
%%     ~~Describe \code{R} here~~
}
}
%% \details{
%% %%  ~~ If necessary, more details than the description above ~~
%% }
\value{A positive semi-definite covariance matrix, for which the variogram
  (or resistance distance) is equal to the input R matrix.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  Hanks and Hooten 2013.  Circuit theory and model-based inference for
  landscape connectivity.  Journal of the American Statistical
  Association.  108(501), 22-33.

  Gower 1996.  Some distance properties of latent root and vector
  methods used in multivariate analysis.  Biometrika 53(3), 325-338.
%% ~put references to the literature/web site here ~
}
\author{Ephraim M. Hanks
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
## create a Wishart covariance matrix with independent structure
Z=matrix(rnorm(10*20),ncol=20,nrow=10)
W=Z\%*\%t(Z)

## convert to resistance distance matrix
D=dist.from.cov(W)

## convert back to covariance matrix
C=cov.from.dist(D)

## compare C and W
max(abs(C-W))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
