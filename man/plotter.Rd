\name{plotter}
\alias{plotter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
A function used to plot the result of IMDE.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotter(x, index = "NRMSE", fileformat = "png", xlab)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{index}{
%%     ~~Describe \code{index} here~~
}
  \item{fileformat}{
%%     ~~Describe \code{fileformat} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
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
function (x, index = "NRMSE", fileformat = "png", xlab) 
{
    x.transform <- cbind(melt(x, variable.name = "Method", value.name = index), 
        xlab)
    if (index == "NRMSE") {
        p1 <- ggplot(x.transform, aes(x = xlab, y = NRMSE, colour = Method)) + 
            geom_line(aes(group = Method), size = 1.5) + geom_point(aes(shape = Method), 
            size = 7) + xlab("Missing Rate (\%)") + ylab("NRMSE")
        if (fileformat == "png") {
            png(filename = "nrmse.png", width = 1600, height = 600, 
                bg = "white", res = NA)
            print(p1)
            dev.off()
        }
        else if (fileformat == "tiff") {
            tiff(filename = "nrmse.tiff", width = 1600, height = 600, 
                bg = "white", res = NA)
            print(p1)
            dev.off()
        }
        else if (fileformat == "jpg") {
            jpeg(filename = "nrmse.jpg", width = 1600, height = 600, 
                bg = "white", res = NA)
            print(p1)
            dev.off()
        }
        else {
            bmp(fileformat == "nrmse.bmp", width = 1600, height = 600, 
                bg = "white", res = NA)
            print(p1)
            dev.off()
        }
    }
    else if (index == "CPP") {
        p1 <- ggplot(x.transform, aes(x = xlab, y = CPP, colour = Method)) + 
            geom_line(aes(group = Method), size = 1.5) + geom_point(aes(shape = Method), 
            size = 7) + scale_x_continuous(breaks = number_ticks(20)) + 
            xlab("Missing Rate (\%)") + ylab("CPP")
        if (fileformat == "png") {
            png(filename = "cpp.png", width = 1600, height = 600, 
                bg = "white", res = NA)
            print(p1)
            dev.off()
        }
        else if (fileformat == "tiff") {
            tiff(filename = "cpp.tiff", width = 1600, height = 600, 
                bg = "white", res = NA)
            print(p1)
            dev.off()
        }
        else if (fileformat == "jpg") {
            jpeg(filename = "cpp.jpg", width = 1600, height = 600, 
                bg = "white", res = NA)
            print(p1)
            dev.off()
        }
        else {
            bmp(fileformat == "cpp.bmp", width = 1600, height = 600, 
                bg = "white", res = NA)
            print(p1)
            dev.off()
        }
    }
    else {
        p1 <- ggplot(x.transform, aes(x = xlab, y = BLCI, colour = Method)) + 
            geom_line(aes(group = Method), size = 1.5) + geom_point(aes(shape = Method), 
            size = 7) + scale_x_continuous(breaks = number_ticks(20)) + 
            xlab("Missing Rate (\%)") + ylab("BLCI")
        if (fileformat == "png") {
            png(filename = "blci.png", width = 1600, height = 600, 
                bg = "white", res = NA)
            print(p1)
            dev.off()
        }
        else if (fileformat == "tiff") {
            tiff(filename = "blci.tiff", width = 1600, height = 600, 
                bg = "white", res = NA)
            print(p1)
            dev.off()
        }
        else if (fileformat == "jpg") {
            jpeg(filename = "blci.jpg", width = 1600, height = 600, 
                bg = "white", res = NA)
            print(p1)
            dev.off()
        }
        else {
            bmp(fileformat == "blci.bmp", width = 1600, height = 600, 
                bg = "white", res = NA)
            print(p1)
            dev.off()
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
