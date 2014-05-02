\name{MVData-class}
\Rdversion{1.1}
\docType{class}
\alias{MVData-class}
\alias{getMissRate,MVData-method}
\alias{getSize,MVData-method}
\alias{show,MVData-method}

\title{Class \code{"MVData"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MVData", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{headname}:}{Object of class \code{"character"} ~~ }
    \item{\code{ori.data}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{miss.idx}:}{Object of class \code{"matrix"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{getMissRate}{\code{signature(x = "MVData")}: ... }
    \item{getSize}{\code{signature(x = "MVData")}: ... }
    \item{show}{\code{signature(object = "MVData")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("MVData")
}
\keyword{classes}
