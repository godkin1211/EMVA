\name{IMPData-class}
\Rdversion{1.1}
\docType{class}
\alias{IMPData-class}
\alias{evaluator,IMPData-method}

\title{Class \code{"IMPData"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("IMPData", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{ori.data}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{imputed.data}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{miss.idx}:}{Object of class \code{"matrix"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{evaluator}{\code{signature(x = "IMPData")}: ... }
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
showClass("IMPData")
}
\keyword{classes}
