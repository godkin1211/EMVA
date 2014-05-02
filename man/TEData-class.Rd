\name{TEData-class}
\Rdversion{1.1}
\docType{class}
\alias{TEData-class}
\alias{impute,TEData-method}
\alias{show,TEData-method}

\title{Class \code{"TEData"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TEData", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{miss.data}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{miss.rate}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{headname}:}{Object of class \code{"character"} ~~ }
    \item{\code{ori.data}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{miss.idx}:}{Object of class \code{"matrix"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{MVData}"}, directly.
}
\section{Methods}{
  \describe{
    \item{impute}{\code{signature(x = "TEData")}: ... }
    \item{show}{\code{signature(object = "TEData")}: ... }
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
showClass("TEData")
}
\keyword{classes}
