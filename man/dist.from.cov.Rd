\name{dist.from.cov}
\alias{dist.from.cov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute a squared distance matrix from a
  covariance matrix.
%%  ~~function to do ... ~~
}
\description{ This computes a squared distance matrix from a covariance
  matrix, or other positive semi-definite matrix.  The resulting squared
  distance matrix is the variogram matrix or the resistance distance
  matrix under a random walk model for connectivity as in Hanks and
  Hooten (2013).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dist.from.cov(Sigma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Sigma}{ A symmetric positive definite matrix.
%%     ~~Describe \code{Sigma} here~~
}
}
%% \details{
%% %%  ~~ If necessary, more details than the description above ~~
%% }
\value{ A negative definite matrix of the same dimensions as Sigma.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Hanks and Hooten 2013.  Circuit theory and model-based inference for
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
## create a Wishart covariance matrix with independent structure

Z=matrix(rnorm(10*20),ncol=20,nrow=10)
W=Z \%*\% t(Z)

## convert to resistance distance matrix
D=dist.from.cov(W)

## convert back to covariance matrix
C=cov.from.dist(D)

## compare C and W
max(abs(C-W))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
