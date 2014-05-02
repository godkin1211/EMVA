\name{BLCIcal}
\alias{BLCIcal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This function is for BLCI calculation.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BLCIcal(mat.com, mat.imp, resp.type = "Pattern discovery", nperms = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat.com}{
%%     ~~Describe \code{mat.com} here~~
}
  \item{mat.imp}{
%%     ~~Describe \code{mat.imp} here~~
}
  \item{resp.type}{
%%     ~~Describe \code{resp.type} here~~
}
  \item{nperms}{
%%     ~~Describe \code{nperms} here~~
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
function (mat.com, mat.imp, resp.type = "Pattern discovery", 
    nperms = 100) 
{
    genenames <- paste("GENE", as.character(1:nrow(mat.com)), 
        sep = "")
    geneid <- as.character(1:nrow(mat.com))
    ans.data <- list(x = mat.com, eigengene.number = 10, geneid = geneid, 
        genenames = genenames)
    imp.data <- list(x = mat.imp, eigengene.number = 10, geneid = geneid, 
        genenames = genenames)
    ans.samr.obj <- samr(ans.data, resp.type = resp.type, nperms = nperms)
    imp.samr.obj <- samr(imp.data, resp.type = resp.type, nperms = nperms)
    ans.delta <- samr.compute.delta.table(ans.samr.obj)
    imp.delta <- samr.compute.delta.table(imp.samr.obj)
    ans.siggene.table <- samr.compute.siggenes.table(ans.samr.obj, 
        del = 0, ans.data, ans.delta, all.genes = T)
    imp.siggene.table <- samr.compute.siggenes.table(imp.samr.obj, 
        del = 0, imp.data, imp.delta, all.genes = T)
    ans.totalgenes <- rbind(ans.siggene.table$genes.lo, ans.siggene.table$genes.up)
    imp.totalgenes <- rbind(imp.siggene.table$genes.lo, imp.siggene.table$genes.up)
    ans.siggene.judge <- as.numeric(ans.totalgenes[, 7])
    imp.siggene.judge <- as.numeric(imp.totalgenes[, 7])
    ans.siggenes <- ans.totalgenes[ans.siggene.judge < 10, , 
        drop = F]
    imp.siggenes <- imp.totalgenes[imp.siggene.judge < 10, , 
        drop = F]
    ans.siggene.list <- ans.siggenes[, 2]
    imp.siggene.list <- imp.siggenes[, 2]
    ans.nonsiggenes <- ans.totalgenes[ans.siggene.judge >= 10, 
        , drop = F]
    imp.nonsiggenes <- imp.totalgenes[imp.siggene.judge >= 10, 
        , drop = F]
    ans.nonsiggene.list <- ans.nonsiggenes[, 2]
    imp.nonsiggene.list <- imp.nonsiggenes[, 2]
    ans.imp.siggenes.int <- intersect(ans.siggene.list, imp.siggene.list)
    ans.imp.nonsiggenes.int <- intersect(ans.nonsiggene.list, 
        imp.nonsiggene.list)
    blci <- length(ans.imp.siggenes.int)/length(ans.siggene.list) + 
        length(ans.imp.nonsiggenes.int)/length(ans.nonsiggene.list) - 
        1
    return(blci)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
