\name{NRMSEcal}
\alias{NRMSEcal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
This function is for NRMSE calculation.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
NRMSEcal(mat.com, mat.imp, miss.idx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat.com}{
%%     ~~Describe \code{mat.com} here~~
}
  \item{mat.imp}{
%%     ~~Describe \code{mat.imp} here~~
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
function (mat.com, mat.imp, miss.idx) 
{
    if (sum(is.na(mat.com)) != 0) 
        stop("There exists NAs in ans!\n")
    if (sum(is.na(mat.imp)) != 0) 
        stop("There exists NAs in imputed data!\n")
    if (!is.logical(miss.idx)) 
        stop("Problems in miss.idx!\n")
    if (length(mat.com[miss.idx]) != length(mat.imp[miss.idx])) 
        stop("Length inequivalent!\n")
    if (!is.matrix(mat.com)) 
        mat.com <- as.matrix(mat.com)
    if (!is.matrix(mat.imp)) 
        mat.com <- as.matrix(mat.com)
    nrmse <- sqrt(mean((mat.com[miss.idx] - mat.imp[miss.idx])^2)/var(mat.com[miss.idx]))
    return(nrmse)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
