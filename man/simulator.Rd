\name{simulator}
\alias{simulator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
This function is used to perform a complete performance-comparison simulation.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~

}
\usage{
simulator(inputData, selected.md = c("zero", "ravg", "knn", "svd", "lls"), runs = 5, ksets = list(knn.k = 10, iknn.k = 10, sknn.k = 10, svd.k = 10, ls.k = 10, lls.k = 50), performance.idx = "nrmse", missRates = c(1, 5, 10, 15, 20), output.format = "png")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputData}{
%%     ~~Describe \code{inputData} here~~
}
  \item{selected.md}{
%%     ~~Describe \code{selected.md} here~~
}
  \item{runs}{
%%     ~~Describe \code{runs} here~~

}
  \item{ksets}{
%%     ~~Describe \code{ksets} here~~
}
  \item{performance.idx}{
%%     ~~Describe \code{performance.idx} here~~
}
  \item{missRates}{
%%     ~~Describe \code{missRates} here~~
}
  \item{output.format}{
%%     ~~Describe \code{output.format} here~~
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
function (inputData, selected.md = c("zero", "ravg", "knn", "svd", 
    "lls"), runs = 5, ksets = list(k1 = 10, k2 = 10, k3 = 10, 
    k4 = 10, k5 = 10, k6 = 50), performance.idx = "nrmse", missRates = c(1, 
    5, 10, 15, 20), output.format = "png") 
{
    if (!is.object(inputData)) 
        stop("Data type error!")
    if (sum(selected.md \%in\% c("zero", "ravg", "knn", "iknn", 
        "sknn", "svd", "ls", "lls", "usr")) != length(selected.md)) 
        stop("Please choose a correct method we provide.")
    if (sum(performance.idx \%in\% c("nrmse", "cpp", "blci")) != 
        length(performance.idx)) 
        stop("Please choose a correct method we provide")
    available.md <- c("zero", "ravg", "knn", "iknn", "sknn", 
        "svd", "ls", "lls", "usr")
    table.construct <- strmacro(index = "NRMSEs", expr = function(x, 
        y) {
        index.table <- as.data.frame(matrix(nr = length(x), nc = length(y)))
        colnames(index.table) <- y
        rownames(index.table) <- x
        return(index.table)
    })
    table.construct.nrmse <- table.construct(index = "NRMSEs")
    table.construct.cpp <- table.construct(index = "CPPs")
    table.construct.blic <- table.construct(index = "BLCIs")
    if ("nrmse" \%in\% performance.idx) {
        NRMSEs.table <- table.construct.nrmse(missRates, selected.md)
        for (i in selected.md) {
            switch(match.arg(i, available.md), zero = {
                zero.nrmse <- rep(0, runs)
            }, ravg = {
                ravg.nrmse <- rep(0, runs)
            }, knn = {
                knn.nrmse <- rep(0, runs)
            }, iknn = {
                iknn.nrmse <- rep(0, runs)
            }, sknn = {
                sknn.nrmse <- rep(0, runs)
            }, svd = {
                svd.nrmse <- rep(0, runs)
            }, ls = {
                ls.nrmse <- rep(0, runs)
            }, lls = {
                lls.nrmse <- rep(0, runs)
            }, usr = {
                usr.nrmse <- rep(0, runs)
            })
        }
    }
    if ("cpp" \%in\% performance.idx) {
        CPPs.table <- table.construct.cpp(missRates, selected.md)
        for (i in selected.md) {
            switch(match.arg(i, available.md), zero = {
                zero.cpp <- rep(0, runs)
            }, ravg = {
                ravg.cpp <- rep(0, runs)
            }, knn = {
                knn.cpp <- rep(0, runs)
            }, iknn = {
                iknn.cpp <- rep(0, runs)
            }, sknn = {
                sknn.cpp <- rep(0, runs)
            }, svd = {
                svd.cpp <- rep(0, runs)
            }, ls = {
                ls.cpp <- rep(0, runs)
            }, lls = {
                lls.cpp <- rep(0, runs)
            }, usr = {
                usr.cpp <- rep(0, runs)
            })
        }
    }
    if ("blci" \%in\% performance.idx) {
        BLCIs.table <- table.construct.blci(missRates, selected.md)
        for (i in selected.md) {
            switch(match.arg(i, available.md), zero = {
                zero.blci <- rep(0, runs)
            }, ravg = {
                ravg.blci <- rep(0, runs)
            }, knn = {
                knn.blci <- rep(0, runs)
            }, iknn = {
                iknn.blci <- rep(0, runs)
            }, sknn = {
                sknn.blci <- rep(0, runs)
            }, svd = {
                svd.blci <- rep(0, runs)
            }, ls = {
                ls.blci <- rep(0, runs)
            }, lls = {
                lls.blci <- rep(0, runs)
            }, usr = {
                usr.blci <- rep(0, runs)
            })
        }
    }
    for (i in 1:length(missRates)) {
        missRate <- missRates[i]/100
        for (j in seq(runs)) {
            testData <- TEdata(inputData, missRate)
            for (m in selected.md) {
                if (m == "zero") {
                  imputedData <- impute(testData, method = "zero")
                  for (p in performance.idx) {
                    performance <- evaluator(imputedData, method = p)
                    if (p == "nrmse") 
                      zero.nrmse[j] <- performance
                    else if (p == "cpp") 
                      zero.cpp[j] <- performance
                    else zero.blci[j] <- performance
                  }
                }
                else if (m == "ravg") {
                  imputedData <- impute(testData, method = "ravg")
                  for (p in performance.idx) {
                    performance <- evaluator(imputedData, method = p)
                    if (p == "nrmse") 
                      ravg.nrmse[j] <- performance
                    else if (p == "cpp") 
                      ravg.cpp[j] <- performance
                    else ravg.blci[j] <- performance
                  }
                }
                else if (m == "knn") {
                  imputedData <- impute(testData, method = "knn", 
                    k = ksets$k1)
                  for (p in performance.idx) {
                    performance <- evaluator(imputedData, method = p)
                    if (p == "nrmse") 
                      knn.nrmse[j] <- performance
                    else if (p == "cpp") 
                      knn.cpp[j] <- performance
                    else knn.blci[j] <- performance
                  }
                }
                else if (m == "sknn") {
                  imputedData <- impute(testData, method = "sknn", 
                    k = ksets$k2)
                  for (p in performance.idx) {
                    performance <- evaluator(imputedData, method = p)
                    if (p == "nrmse") 
                      sknn.nrmse[j] <- performance
                    else if (p == "cpp") 
                      sknn.cpp[j] <- performance
                    else sknn.blci[j] <- performance
                  }
                }
                else if (m == "iknn") {
                  imputedData <- impute(testData, method = "iknn", 
                    k = ksets$k3)
                  for (p in performance.idx) {
                    performance <- evaluator(imputedData, method = p)
                    if (p == "nrmse") 
                      iknn.nrmse[j] <- performance
                    else if (p == "cpp") 
                      iknn.cpp[j] <- performance
                    else iknn.blci[j] <- performance
                  }
                }
                else if (m == "svd") {
                  imputedData <- impute(testData, method = "svd", 
                    k = ksets$k4)
                  for (p in performance.idx) {
                    performance <- evaluator(imputedData, method = p)
                    if (p == "nrmse") 
                      svd.nrmse[j] <- performance
                    else if (p == "cpp") 
                      svd.cpp[j] <- performance
                    else svd.blci[j] <- performance
                  }
                }
                else if (m == "ls") {
                  imputedData <- impute(testData, method = "ls", 
                    k = ksets$k5)
                  for (p in performance.idx) {
                    performance <- evaluator(imputedData, method = p)
                    if (p == "nrmse") 
                      ls.nrmse[j] <- performance
                    else if (p == "cpp") 
                      ls.cpp[j] <- performance
                    else ls.blci[j] <- performance
                  }
                }
                else if (m == "lls") {
                  imputedData <- impute(testData, method = "lls", 
                    k = ksets$k6)
                  for (p in performance.idx) {
                    performance <- evaluator(imputedData, method = p)
                    if (p == "nrmse") 
                      lls.nrmse[j] <- performance
                    else if (p == "cpp") 
                      lls.cpp[j] <- performance
                    else lls.blci[j] <- performance
                  }
                }
                else if (m == "usr") {
                  imputedData <- impute(testData, method = "usr")
                  for (p in performance.idx) {
                    performance <- evaluator(imputedData, method = p)
                    if (p == "nrmse") 
                      usr.nrmse[j] <- performance
                    else if (p == "cpp") 
                      usr.cpp[j] <- performance
                    else usr.blci[j] <- performance
                  }
                }
                else {
                }
            }
        }
        for (p in performance.idx) {
            for (m in selected.md) {
                if (p == "nrmse") {
                  eval(parse(text = paste("NRMSEs.table$", m, 
                    "[", i, "]<-mean(", m, ".", p, ")", sep = "")))
                }
                else if (p == "cpp") {
                  eval(parse(text = paste("CPPs.table$", m, "[", 
                    i, "]<-mean(", m, ".", p, ")", sep = "")))
                }
                else {
                  eval(parse(text = paste("BLCIs.table$", m, 
                    "[", i, "]<-mean(", m, ".", p, ")", sep = "")))
                }
            }
        }
    }
    for (p in performance.idx) {
        if (p == "nrmse") 
            plotter(NRMSEs.table, index = "NRMSE", fileformat = output.format, 
                xlab = missRates)
        else if (p == "cpp") 
            plotter(CPPs.table, index = "CPP", fileformat = output.format, 
                xlab = missRates)
        else plotter(BLCIs.table, index = "BLCI", fileformat = output.format, 
            xlab = missRates)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
