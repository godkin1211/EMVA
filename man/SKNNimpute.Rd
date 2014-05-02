\name{SKNNimpute}
\alias{SKNNimpute}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
SKNNimpute.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SKNNimpute(x, k = 10, sim.method = "euclidean")
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
function (x, k = 10, sim.method = "euclidean") 
{
    rowNum <- nrow(x)
    colNum <- ncol(x)
    missIdx <- is.na(x)
    miss.RowIdx <- which(rowSums(missIdx) != 0)
    x.completeRows <- x[-miss.RowIdx, ]
    x.missingRows <- x[miss.RowIdx, ]
    miss.list <- order(rowSums(is.na(x.missingRows)))
    for (i in seq(nrow(x.missingRows))) {
        target <- x.missingRows[miss.list[i], ]
        dist.list <- similarityCal(target, x.completeRows, sim.method)
        neighborsIdx <- order(abs(dist.list), decreasing = T)[1:k]
        missColIdx <- which(is.na(target))
        estimation <- sapply(missColIdx, function(j) {
            weight <- dist.list[neighborsIdx]
            weightedAvg <- weight \%*\% x.completeRows[neighborsIdx, 
                j]/sum(weight)
            return(weightedAvg)
        })
        target[missColIdx] <- estimation
        x.missingRows[miss.list[i], ] <- target
        x.completeRows <- rbind(x.completeRows, target)
    }
    x[miss.RowIdx, ] <- x.missingRows
    return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
