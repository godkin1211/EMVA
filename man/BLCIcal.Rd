\name{BLCIcal}
\alias{BLCIcal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculates biomarker list concordance index for differentially expressed gene detection.}
\description{
Implementation of BLCI computation. See S. Oh et al., 2011.
}
\usage{
BLCIcal(mat.com, mat.imp, resp.type = "Pattern discovery", nperms = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat.com}{Complete data matrix.}
  \item{mat.imp}{Imputed data matrix.}
  \item{resp.type}{ Setting for problem type, as in \code{\link{samr}.}
}
  \item{nperms}{ Number of permutations used to estimate discovery rates, as in \code{\link{samr}.}
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
Sunghee Oh, Dongwan D. Kang, Guy N. Brock and George C. Tseng. Biological impact of missing-value imputation on downstream analyses of gene expression profiles. \emph{Bioinformatics}, 27(1): 78-86, 2011.
}
\author{
Chiachun Chiu and Weisheng Wu.
}
\note{Requires \code{\link{samr}}.}

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
\keyword{ differentially expressed genes }
\keyword{ missing value imputation }% __ONLY ONE__ keyword per line
