\name{USRimpute}
\alias{USRimpute}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
User-defined imputation function.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
USRimpute(x, lambda, verbose = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (x, lambda, verbose = F) 
{
    missing.matrix <- is.na(x)
    numMissing <- sum(missing.matrix)
    print(paste("imputing on", numMissing, "missing values with matrix size", 
        nrow(x) * ncol(x), sep = " "))
    if (numMissing == 0) 
        return(x)
    missing.cols.indices <- which(apply(missing.matrix, 2, function(i) any(i)))
    x.missing <- (rbind(1:ncol(x), x))[, missing.cols.indices]
    x.missing.imputed <- apply(x.missing, 2, function(j) {
        colIndex = j[1]
        j.original = j[-1]
        missing.rows = which(missing.matrix[, colIndex])
        if (length(missing.rows) == nrow(x)) 
            warning(paste("Column", colIndex, "is completely missing", 
                sep = " "))
        j.original[missing.rows] = mean(j.original[-missing.rows])
        j.original
    })
    x[, missing.cols.indices] <- x.missing.imputed
    missing.matrix2 <- is.na(x)
    x[missing.matrix2] <- 0
    x.svd <- svd(x)
    lambda.indices <- which(x.svd$d < lambda)
    if (length(lambda.indices) > 0) 
        d.augmented <- c(x.svd$d[-lambda.indices] - lambda, rep(0, 
            length(lambda.indices)))
    else d.augmented <- x.svd$d - lambda
    x.imputed <- x.svd$u \%*\% diag(d.augmented) \%*\% x.svd$v
    return(list(x = x.imputed, missing.matrix = missing.matrix))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
