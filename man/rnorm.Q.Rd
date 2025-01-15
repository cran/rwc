\name{rnorm.Q}
\alias{rnorm.Q}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sample random normal variables with precision matrix Q.
%%  ~~function to do ... ~~
}
\description{General function to make use of sparse matrix methods to
  efficiently simulate random multivariate normal random variables with
  sparse precision matrices.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rnorm.Q(Q, mu = rep(0, nrow(Q)), X = Matrix(1, nrow = nrow(Q), ncol =
1),
zero.constraint = FALSE, canon = FALSE, diag.adjust = .Machine$double.eps * 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Q}{Precision matrix, defined as a sparse Matrix object.
%%     ~~Describe \code{Q} here~~
}
  \item{mu}{Mean vector of dimension equal to the dimension of Q.
%%     ~~Describe \code{mu} here~~
}
  \item{X}{Matrix of vectors which should be orthogonal to the simulated
    random variable.  
%%     ~~Describe \code{X} here~~
}
  \item{zero.constraint}{If TRUE, then the simulated random variable is
    orthogonal to the columns of X.
%%     ~~Describe \code{zero.constraint} here~~
}
  \item{canon}{If TRUE, then draw from the 'canonical form'.
%%     ~~Describe \code{canon} here~~
}
  \item{diag.adjust}{Numeric value to be added to the diagonal of Q to
    make it positive definite.
%%     ~~Describe \code{diag.adjust} here~~
}
}
\details{In the 'canonical form', the variable is drawn from:

  v~N(Q^-1 mu, Q^-1)

  In the non-canonical form, the variable is drawn from

  v~N( mu, Q^-1)

  
%%  ~~ If necessary, more details than the description above ~~
}
\value{A vector of the resulting random variable.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Rue and Held 2005.  Gaussian Markov Random Fields: theory
  and applications.  Chapman and Hall.
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
ras=raster(nrow=30,ncol=30)
extent(ras) <- c(0,30,0,30)
values(ras) <- 1

int=ras
cov.ras=ras


## simulate "resistance" raster
B.int=get.TL(int)
Q.int=get.Q(B.int,1)
values(cov.ras) <- as.numeric(rnorm.Q(Q.int,zero.constraint=TRUE))
plot(cov.ras)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
