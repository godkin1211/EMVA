\name{KNNimpute}
\alias{KNNimpute}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
KNN
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
KNNimpute(x, k = 15, sim.method = "euclidean")
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
function (x, k = 15, sim.method = "euclidean") 
{
    miss.idx <- is.na(x)
    x.missRowIdx <- which(rowSums(miss.idx) != 0)
    x.miss <- x[x.missRowIdx, ]
    x.comp <- x[-x.missRowIdx, ]
    x.imp <- t(apply(x.miss, 1, function(i) {
        missidx <- which(is.na(i))
        vec.usable <- i[-missidx]
        neighbor.pool <- x.comp[, -missidx]
        similarities <- similarityCal(vec.usable, neighbor.pool, 
            method = sim.method)
        neighborhood <- order(similarities, decreasing = TRUE)[1:k]
        imp.values <- sapply(missidx, function(j) weight.avg <- similarities[neighborhood] \%*\% 
            x.comp[neighborhood, j]/sum(similarities[neighborhood]))
        i[missidx] <- imp.values
        return(i)
    }))
    x[x.missRowIdx, ] <- x.imp
    return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
