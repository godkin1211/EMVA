\name{IMPdata}
\alias{IMPdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
The constructor of IMPData object.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
IMPdata(x.ori, x.imputed, miss.idx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.ori}{
%%     ~~Describe \code{x.ori} here~~
}
  \item{x.imputed}{
%%     ~~Describe \code{x.imputed} here~~
}
  \item{miss.idx}{
%%     ~~Describe \code{miss.idx} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x.ori, x.imputed, miss.idx) 
{
    new("IMPData", ori.data = x.ori, imputed.data = x.imputed, 
        miss.idx = miss.idx)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
