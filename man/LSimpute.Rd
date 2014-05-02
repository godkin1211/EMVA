\name{LSimpute}
\alias{LSimpute}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
LS
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LSimpute(x, k = 15, sim.method = "pearson", e = 1e-06)
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
function (x, k = 15, sim.method = "pearson", e = 1e-06) 
{
    rowNum <- nrow(x)
    colNum <- ncol(x)
    missIdx <- is.na(x)
    miss.RowIdx <- which(rowSums(missIdx) != 0)
    x.completeRows <- x[-miss.RowIdx, ]
    x.missingRows <- x[miss.RowIdx, ]
    x.imputed <- t(apply(x.missingRows, 1, function(i) {
        missColIdx <- which(is.na(i))
        dist.list <- similarityCal(i[-missColIdx], x.completeRows[, 
            -missColIdx], method = sim.method)
        neighborsIdx <- order(dist.list, decreasing = T)[1:k]
        neighbor.sim <- dist.list[neighborsIdx]
        weight <- ((neighbor.sim^2)/(1 - (neighbor.sim)^2 + e))^2
        fit.list <- lapply(neighborsIdx, function(m) {
            fit <- lm(i[-missColIdx] ~ x.completeRows[m, -missColIdx])
            reg <- c(fit$coefficients[[1]], fit$coefficients[[2]])
            return(reg)
        })
        estimation <- sapply(missColIdx, function(j) {
            regression <- rep(0, k)
            for (n in seq(k)) {
                a0 <- fit.list[n][[1]][1]
                b0 <- fit.list[n][[1]][2]
                regression[n] <- a0 + b0 * x.completeRows[neighborsIdx[n], 
                  j]
            }
            weightedAvg <- (weight \%*\% regression)/sum(weight)
            return(weightedAvg)
        })
        i[missColIdx] <- estimation
        return(i)
    }))
    x[miss.RowIdx, ] <- x.imputed
    return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
