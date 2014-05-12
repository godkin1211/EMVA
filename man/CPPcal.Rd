\name{CPPcal}
\alias{CPPcal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates CPP (Conserved Pair Proportion).
}
\description{
An implementation of CPP. See A. G. de Brevern et al., 2004.
}
\usage{
CPPcal(mat.com, mat.imp, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat.com}{Complete data matrix.}
  \item{mat.imp}{Imputed data matrix.}
  \item{k}{Number of groups, as in \code{\link{cutree}.}}
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
Alexandre G de Brevern, Serge Hazout and Alain Malpertuy. Influence of microarray experiments missing values on the stability of gene groups by hierarchical clustering. \emph{BMC Bioinformatics}, 5:114.
}
\author{Chiachun Chiu and Weisheng Wu.}
\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (mat.com, mat.imp, k) 
{
    #ans.clustered <- kmeans(mat.com, k, iter.max)
    ans.clustered<-hclust(dist(mat.com),method="complete")
    #ans.clust <- ans.clustered$cluster
    ans.clust<-cutree(ans.clustered,k=k)
    #imp.clustered <- kmeans(mat.imp, k, iter.max)
    imp.clustered<-hclust(dist(mat.imp),method="complete")
    #imp.clust <- imp.clustered$cluster
    imp.clust<-cutree(imp.clustered,k=k)
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
