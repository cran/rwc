\name{get.Phi}
\alias{get.Phi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute the precision matrix Phi of observed nodes on a graph.
%%  ~~function to do ... ~~
}
\description{Given a Gaussian Markov random field defined by a precision
  matrix Q, this returns Phi, which is the precision matrix of the nodes
  indexed by 'obs.idx', computed using the Schur complement.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get.Phi(Q, obs.idx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Q}{A precision matrix
%%     ~~Describe \code{Q} here~~
}
  \item{obs.idx}{A vector of unique indices of observed nodes in the
    graph defined by Q.  
%%     ~~Describe \code{obs.idx} here~~
}
}
%% \details{
%% %%  ~~ If necessary, more details than the description above ~~
%% }
\value{A square matrix of dimension equal to the length of 'obs.idx'.
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
TL.int=get.TL(int)
Q.int=get.Q(TL.int,1)
## get precision matrix of only the first 5 nodes
get.Phi(Q.int,1:5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
