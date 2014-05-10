simulator<-function(inputData,
                    selected.md=c("zero","ravg","knn","svd","lls"),
                    runs=5,
                    ksets=list(knn.k=10,iknn.k=10,sknn.k=10,svd.k=10,ls.k=10,lls.k=50),
                    performance.idx="nrmse",
                    missRates=c(1,5,10,15,20),
                    output.format="png",
					email="",
					passwd=""
					)
{

  # Check each argument
    if (!is.object(inputData)) stop("Data type error!")
    if (sum(selected.md %in% c("zero","ravg","knn","iknn","sknn","svd","ls","lls","usr")) != length(selected.md))
        stop("Please choose a correct method we provide.")
    if (sum(performance.idx %in% c("nrmse","cpp","blci")) != length(performance.idx))
        stop("Please choose a correct method we provide")
  # Construct an adequate data-structure to store the results from each performance evaluator.
    available.md <- c("zero","ravg","knn","iknn","sknn","svd","ls","lls","usr")
    table.construct <- strmacro(index="NRMSEs",
                                expr= function(x,y){
                                    index.table <- as.data.frame(matrix(nr=length(x),nc=length(y)))
                                    colnames(index.table) <- y
                                    rownames(index.table) <- x
                                    return(index.table)
                                    })

    table.construct.nrmse <- table.construct(index="NRMSEs")
    table.construct.cpp<-table.construct(index="CPPs")
    table.construct.blci<-table.construct(index="BLCIs")

    if ("nrmse" %in% performance.idx) {
        NRMSEs.table <- table.construct.nrmse(missRates, selected.md)
        for (i in selected.md) {
            switch(match.arg(i, available.md),
                zero = {zero.nrmse <- rep(0, runs)},
                ravg = {ravg.nrmse <- rep(0, runs)},
                knn = {knn.nrmse <- rep(0, runs)},
                iknn = {iknn.nrmse <- rep(0, runs)},
                sknn = {sknn.nrmse <- rep(0, runs)},
                svd = {svd.nrmse <- rep(0, runs)},
                ls = {ls.nrmse <- rep(0, runs)},
                lls = {lls.nrmse <- rep(0, runs)},
                usr = {usr.nrmse <- rep(0, runs)})
        }
    }

    if ("cpp" %in% performance.idx) {
        CPPs.table <- table.construct.cpp(missRates, selected.md)
        for (i in selected.md) {
            switch(match.arg(i, available.md),
                zero = {zero.cpp <- rep(0, runs)},
                ravg = {ravg.cpp <- rep(0, runs)},
                knn = {knn.cpp <- rep(0, runs)},
                iknn = {iknn.cpp <- rep(0, runs)},
                sknn = {sknn.cpp <- rep(0, runs)},
                svd = {svd.cpp <- rep(0, runs)},
                ls = {ls.cpp <- rep(0, runs)},
                lls = {lls.cpp <- rep(0, runs)},
                usr = {usr.cpp <- rep(0, runs)})
        }
    }

    if ("blci" %in% performance.idx) {
        BLCIs.table <- table.construct.blci(missRates, selected.md)
        for (i in selected.md) {
            switch(match.arg(i, available.md),
                zero = {zero.blci <- rep(0, runs)},
                ravg = {ravg.blci <- rep(0, runs)},
                knn = {knn.blci <- rep(0, runs)},
                iknn = {iknn.blci <- rep(0, runs)},
                sknn = {sknn.blci <- rep(0, runs)},
                svd = {svd.blci <- rep(0, runs)},
                ls = {ls.blci <- rep(0, runs)},
                lls = {lls.blci <- rep(0, runs)},
                usr = {usr.blci <- rep(0, runs)})
        }
    }
    
    ## Simulation procedure
    for (i in 1:length(missRates)) {
        
        missRate <- missRates[i]/100
        
        for (j in seq(runs)){
            # Generating testing data
			cat("\n==============\nMissing Rate:",missRates[i],"%\n==============\n\n")
			cat("\n########\n# Run:",j,"\n########\n")
            testData <- TEdata(inputData, missRate)
            for (m in selected.md) {
                if (m == "zero") {
                    cat("***************\n* ZEROimpute.......\n***************\n")
                    imputedData <- impute(testData, method="zero")
                    for (p in performance.idx) {
                        performance <-evaluator(imputedData,method=p)
                        if (is.na(performance) | is.nan(performance)) stop("Performance is NA!")
                        if (p == "nrmse") zero.nrmse[j] <- performance
                        else if (p == "cpp") zero.cpp[j] <- performance
                        else zero.blci[j] <- performance
                    }
                } else if (m == "ravg") {
                    cat("***************\n* RAVGimpute.......\n***************\n")
                    imputedData <- impute(testData, method="ravg")
                    for (p in performance.idx) {
                        performance <-evaluator(imputedData,method=p)
                        if (is.na(performance) | is.nan(performance)) stop("Performance is NA!")
                        if (p == "nrmse") ravg.nrmse[j] <- performance
                        else if (p == "cpp") ravg.cpp[j] <- performance
                        else ravg.blci[j] <- performance
                    }
                } else if (m == "knn") {
                    cat("***************\n* KNNimpute.......\n***************\n")
                    imputedData <- impute(testData, method="knn", k=ksets$knn.k)
                    for (p in performance.idx) {
                        performance <-evaluator(imputedData,method=p)
                        if (is.na(performance) | is.nan(performance)) stop("Performance is NA!")
                        if (p == "nrmse") knn.nrmse[j] <- performance
                        else if (p == "cpp") knn.cpp[j] <- performance
                        else knn.blci[j] <- performance
                    }
                } else if (m == "sknn") {
                    cat("***************\n* SKNNimpute.......\n***************\n")
                    imputedData <- impute(testData, method="sknn", k=ksets$sknn.k)
                    for (p in performance.idx) {
                        performance <-evaluator(imputedData,method=p)
                        if (is.na(performance) | is.nan(performance)) stop("Performance is NA!")
                        if (p == "nrmse") sknn.nrmse[j] <- performance
                        else if (p == "cpp") sknn.cpp[j] <- performance
                        else sknn.blci[j] <- performance
                    }
                } else if (m == "iknn") {
                    cat("***************\n* IKNNimpute.......\n***************\n")
                    imputedData <- impute(testData, method="iknn", k=ksets$iknn.k)
                    for (p in performance.idx) {
                        performance <-evaluator(imputedData,method=p)
                        if (is.na(performance) | is.nan(performance)) stop("Performance is NA!")
                        if (p == "nrmse") iknn.nrmse[j] <- performance
                        else if (p == "cpp") iknn.cpp[j] <- performance
                        else iknn.blci[j] <- performance
                    }
                } else if (m == "svd") {
                    cat("***************\n* SVDimpute.......\n***************\n")
                    imputedData <- impute(testData, method="svd", k=ksets$svd.k)
                    for (p in performance.idx) {
                        performance <-evaluator(imputedData,method=p)
                        if (is.na(performance) | is.nan(performance)) stop("Performance is NA!")
                        if (p == "nrmse") svd.nrmse[j] <- performance
                        else if (p == "cpp") svd.cpp[j] <- performance
                        else svd.blci[j] <- performance
                    }
                } else if (m == "ls") {
                    cat("***************\n* LSimpute.......\n***************\n")
                    imputedData <- impute(testData, method="ls", k=ksets$ls.k)
                    for (p in performance.idx) {
                        performance <-evaluator(imputedData,method=p)
                        if (is.na(performance) | is.nan(performance)) stop("Performance is NA!")
                        if (p == "nrmse") ls.nrmse[j] <- performance
                        else if (p == "cpp") ls.cpp[j] <- performance
                        else ls.blci[j] <- performance
                    }
                } else if (m == "lls") {
                    cat("***************\n* LLSimpute.......\n***************\n")
                    imputedData <- impute(testData, method="lls", k=ksets$lls.k)
                    for (p in performance.idx) {
                        performance <-evaluator(imputedData,method=p)
                        if (is.na(performance) | is.nan(performance)) stop("Performance is NA!")
                        if (p == "nrmse") lls.nrmse[j] <- performance
                        else if (p == "cpp") lls.cpp[j] <- performance
                        else lls.blci[j] <- performance
                    }
                } else {
                    cat("***************\n* USRimpute.......\n***************\n")
                    imputedData <- impute(testData, method="usr")
                    for (p in performance.idx) {
                        performance <-evaluator(imputedData,method=p)
                        if (is.na(performance) | is.nan(performance)) stop("Performance is NA!")
                        if (p == "nrmse") usr.nrmse[j] <- performance
                        else if (p == "cpp") usr.cpp[j] <- performance
                        else usr.blci[j] <- performance
                    }
                }
            }
        }

        for (p in performance.idx) {
            for (m in selected.md) {
                if (p == "nrmse") {
                    if (eval(parse(text=paste("mean(",m,".",p,")",sep=""))) == NA | eval(parse(text=paste("mean(",m,".",p,")",sep=""))) == NaN) stop("mean(m,p) is missing!")
                    eval(parse(text=paste("NRMSEs.table$",m,"[",i,"]<-mean(",m,".",p,")",sep="")))
                    cat("Results:\n")
                    print(NRMSEs.table)
                } else if (p == "cpp") {
                    eval(parse(text=paste("CPPs.table$",m,"[",i,"]<-mean(",m,".",p,")",sep="")))
                    cat("Results:\n")
                    print(CPPs.table)
                } else {
                    eval(parse(text=paste("BLCIs.table$",m,"[",i,"]<-mean(",m,".",p,")",sep="")))
                    cat("Results:\n")
                    print(BLCIs.table)
                }
            }
        }
    }
    # Plotting
    for (p in performance.idx) {
        if (p == "nrmse")
            plotter(NRMSEs.table, index="NRMSE", fileformat=output.format, xlab=missRates)
        else if (p == "cpp")
            plotter(CPPs.table, index="CPP", fileformat=output.format, xlab=missRates)
        else
            plotter(BLCIs.table, index="BLCI", fileformat=output.format, xlab=missRates)
    }

	# Job completion notification
	if (email != "" && passwd != "" ) send.mail(from = email, to = email,
											subject = "IMDE: job-completion notification",
											body = "Your job has been completed!",
											smtp = list(host.name="smtp.gmail.com", port=465, user.name = email, passwd=passwd, ssl=TRUE),
										    authenticate = TRUE, send =TRUE)
}
