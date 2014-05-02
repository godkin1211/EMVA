\name{MVdata}
\alias{MVdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
The constructor of MVData object.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MVdata(inputfile = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputfile}{
%%     ~~Describe \code{inputfile} here~~
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
function (inputfile = NULL) 
{
    if (inputfile != "" || inputfile != NULL) {
        headname <- readLines(inputfile)[1]
        ori.data <- as.matrix(read.table(inputfile, header = T, 
            row.names = 1, sep = "\t", skip = 1))
        miss.idx <- is.na(ori.data)
        new("MVData", headname = headname, ori.data = ori.data, 
            miss.idx = miss.idx)
    }
    else {
        data(lymphoma)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
