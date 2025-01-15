\name{get.Q}
\alias{get.Q}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create a precision matrix from a transition list and a set of
  log-conductance rates.
%%  ~~function to do ... ~~
}
\description{Creates a precision matrix Q, with off diagonal entries
  equal to exp(b1*x_1ij + ... + bp*x_pij), where beta=(b1,b2,...,bp) is
  a vector of log-conductance values of the covariates x_kij.  Each
  x_kij is equal to (x_ki+x_kj)/2.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get.Q(TL, beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TL}{A transition list from TL.from.stack
%%     ~~Describe \code{B} here~~
}
  \item{beta}{A vector of log-conductance rates with length equal to the
    length of TL.
%%     ~~Describe \code{beta} here~~
}
}
%% \details{
%% %%  ~~ If necessary, more details than the description above ~~
%% }
\value{A precision matrix, as a sparse matrix of class 'dgCMatrix',
  with dimension equal to n^2 by n^2, where n is the number of nodes in
  the raster stack used to compute TL.
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
int=raster(nrow=30,ncol=30)
values(int)=1
B.int=get.TL(int)
Q.int=get.Q(B.int,1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
