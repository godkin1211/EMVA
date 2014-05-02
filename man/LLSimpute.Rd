\name{LLSimpute}
\alias{LLSimpute}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
LLS
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LLSimpute(x, k = 50, sim.method = "euclidean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{sim.method}{
%%     ~~Describe \code{sim.method} here~~
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
function (x, k = 50, sim.method = "euclidean") 
{
    missidx <- is.na(x)
    missRowIdx <- which(rowSums(missidx) != 0)
    x.completePart <- x[-missRowIdx, ]
    x.missingPart <- x[missRowIdx, ]
    imputed <- t(apply(x.missingPart, 1, function(i) {
        missColIdx <- which(is.na(i))
        dist.list <- similarityCal(i[-missColIdx], x.completePart[, 
            -missColIdx], sim.method)
        neighborIdx <- order(dist.list, decreasing = T)[1:k]
        A <- x.completePart[neighborIdx, -missColIdx, drop = FALSE]
        b <- x.completePart[neighborIdx, missColIdx, drop = FALSE]
        Cp <- ginv(t(A))
        w <- i[-missColIdx, drop = FALSE]
        X <- Cp \%*\% w
        ans <- t(b) \%*\% X
        i[missColIdx] <- ans
        return(i)
    }))
    x[missRowIdx, ] <- imputed
    cat("LLS imputation completed!\n")
    return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
