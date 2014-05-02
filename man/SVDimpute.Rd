\name{SVDimpute}
\alias{SVDimpute}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
SVD
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SVDimpute(x, k = 15, iters = 10, sim.method = "euclidean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{iters}{
%%     ~~Describe \code{iters} here~~
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
function (x, k = 15, iters = 10, sim.method = "euclidean") 
{
    missing.matrix <- is.na(x)
    numMissing <- sum(missing.matrix)
    print(paste("imputing on", numMissing, "missing values with matrix size", 
        nrow(x) * ncol(x), sep = " "))
    if (numMissing == 0) 
        return(x)
    miss.colIdx <- which(apply(missing.matrix, 2, function(i) any(i)))
    x.missing <- (rbind(1:ncol(x), x))[, miss.colIdx]
    x.imputed <- apply(x.missing, 2, function(j) {
        colIdx <- j[1]
        j.original <- j[-1]
        missing.rows <- which(missing.matrix[, colIdx])
        if (length(missing.rows) == nrow(x)) 
            warning(paste("Column", colIdx, "is completely missing", 
                sep = " "))
        j.original[missing.rows] <- mean(j.original[-missing.rows])
        j.original
    })
    x[, miss.colIdx] <- x.imputed
    missing.matrix2 <- is.na(x)
    x[missing.matrix2] <- 0
    for (i in 1:iters) {
        print(paste("Running iteration", i, sep = " "))
        x.svd <- svd(x, nu = k, nv = k)
        x.svd <- x.svd$u \%*\% diag(x.svd$d[1:k], nrow = k, ncol = k) \%*\% 
            t(x.svd$v)
        x[missing.matrix] <- x.svd[missing.matrix]
    }
    return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
