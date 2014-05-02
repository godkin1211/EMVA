\name{CPPcal}
\alias{CPPcal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
This function is for CPP calculation.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CPPcal(mat.com, mat.imp, k, iter.max = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat.com}{
%%     ~~Describe \code{mat.com} here~~
}
  \item{mat.imp}{
%%     ~~Describe \code{mat.imp} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{iter.max}{
%%     ~~Describe \code{iter.max} here~~
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
function (mat.com, mat.imp, k, iter.max = 10) 
{
    ans.clustered <- kmeans(mat.com, k, iter.max)
    ans.clust <- ans.clustered$cluster
    imp.clustered <- kmeans(mat.imp, k, iter.max)
    imp.clust <- imp.clustered$cluster
    cpptable <- matrix(nr = k, nc = k)
    rownames(cpptable) <- sprintf("C.ans.\%i", 1:k)
    colnames(cpptable) <- sprintf("C.imp.\%i", 1:k)
    for (i in 1:k) {
        ans.genes.inthiscluster <- which(ans.clust == i)
        for (j in 1:k) {
            imp.genes.inthiscluster <- which(imp.clust == j)
            common.member <- length(intersect(ans.genes.inthiscluster, 
                imp.genes.inthiscluster))
            cpptable[i, j] <- common.member
        }
    }
    cpp <- sum(apply(cpptable, 1, max))/sum(cpptable)
    return(cpp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
