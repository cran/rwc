\name{get.TL}
\alias{get.TL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Construct a transition list from a raster or raster stack
%%  ~~function to do ... ~~
}
\description{This computes a list of log-transition matrices as a
  preliminary step to creating a precision matrix from covariate rasters.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get.TL(rast.stack)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rast.stack}{ A raster layer or raster stack object.
%%     ~~Describe \code{rast.stack} here~~
}
}
%% \details{ 
%% %%  ~~ If necessary, more details than the description above ~~
%% }
\value{ A list of length equal to the number of raster
  layers in rast.stack.  Each element in the list is a sparse Matrix of
  class 'dgCMatrix'.
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
image(Q.int)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
