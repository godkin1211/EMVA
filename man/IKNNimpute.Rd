\name{IKNNimpute}
\alias{IKNNimpute}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
IKNN
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
IKNNimpute(x, k = 10, sim.method = "euclidean", iter = 3, e = 0.001)
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
  \item{iter}{
%%     ~~Describe \code{iter} here~~
}
  \item{e}{
%%     ~~Describe \code{e} here~~
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
function (x, k = 10, sim.method = "euclidean", iter = 3, e = 0.001) 
{
    cat(sprintf("k = \%g, iter= \%g, e = \%g\n", k, iter, e))
    missIdx <- is.na(x)
    rowNum <- nrow(x)
    miss.RowIdx <- which(rowSums(missIdx) != 0)
    x.ravged <- RAVGimpute(x)
    cat("Row average imputation completed!\n")
    x.miss <- (cbind(1:rowNum, x))[miss.RowIdx, ]
    x.complete <- x.ravged
    cat("The size of x.complete:", dim(x.complete), "\n")
    err <- 99
    for (r in 1:iter) {
        x.old <- x.complete
        err.old <- err
        cat(sprintf("Start the \%g cycle of iknn imputation\n", 
            r))
        x.imputed <- t(apply(x.miss, 1, function(j) {
            rowIdx <- j[1]
            j.origin <- j[-1]
            neighbor.pool <- x.complete[-rowIdx, ]
            target <- x.complete[rowIdx, ]
            dist.list <- similarityCal(target, neighbor.pool, 
                method = sim.method)
            neighborsIdx <- order(dist.list, decreasing = T)[1:k]
            missColIdx <- which(is.na(j.origin))
            estimation <- sapply(missColIdx, function(h) {
                weight <- dist.list[neighborsIdx]
                weightedAvg <- weight \%*\% neighbor.pool[neighborsIdx, 
                  h]/sum(weight)
                return(weightedAvg)
            })
            j.origin[missColIdx] <- estimation
            return(j.origin)
        }))
        cat(sprintf("The \%g cycle imputation has been completed\n", 
            r))
        x.complete[miss.RowIdx, ] <- x.imputed
        err <- sum((x.old[missIdx] - x.complete[missIdx])^2)
        cat("err:", err, "\n")
        if (r > 1) 
            if ((err.old/err < 4) | (err < e)) 
                break
    }
    x <- x.complete
    return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
