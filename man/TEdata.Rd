\name{TEdata}
\alias{TEdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
The constructor of TEData object.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TEdata(MVdata.obj, miss.rate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MVdata.obj}{
%%     ~~Describe \code{MVdata.obj} here~~
}
  \item{miss.rate}{
%%     ~~Describe \code{miss.rate} here~~
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
function (MVdata.obj, miss.rate) 
{
    headname <- MVdata.obj@headname
    missRow.idx <- which(rowSums(MVdata.obj@miss.idx) != 0)
    ori.data <- MVdata.obj@ori.data[-missRow.idx, ]
    miss.data <- ori.data
    num.fixxed.row <- 300
    fixxed.row.idx <- sample(seq(nrow(miss.data)), num.fixxed.row)
    fixxed.row.idx <- fixxed.row.idx[order(fixxed.row.idx)]
    temp.data <- miss.data[-fixxed.row.idx, ]
    miss.data.size <- nrow(miss.data) * ncol(miss.data)
    MV.num <- round(miss.data.size * miss.rate)
    temp.data[sample(seq(nrow(temp.data) * ncol(temp.data)), 
        MV.num)] <- NA
    miss.data[-fixxed.row.idx, ] <- temp.data
    miss.idx <- is.na(miss.data)
    new("TEData", headname = headname, ori.data = ori.data, miss.idx = miss.idx, 
        miss.rate = miss.rate, miss.data = miss.data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
