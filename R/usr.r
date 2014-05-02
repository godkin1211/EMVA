USRimpute <- function(x, lambda, verbose=F) {
    missing.matrix <- is.na(x)
    numMissing <- sum(missing.matrix)
    print(paste("imputing on", numMissing, "missing values with matrix size", nrow(x) * ncol(x), sep=" "))
    if(numMissing == 0) return (x) 
    missing.cols.indices <- which(apply(missing.matrix, 2, function(i) any(i) )) 
    x.missing <- (rbind(1:ncol(x), x))[,missing.cols.indices]
    x.missing.imputed <- apply(x.missing, 2, function(j) {
        colIndex = j[1]
        j.original = j[-1]
        missing.rows = which(missing.matrix[,colIndex])
        if(length(missing.rows) == nrow(x)) warning( paste("Column",colIndex,"is completely missing",sep=" ") )
        j.original[missing.rows] = mean(j.original[-missing.rows])
        j.original
    })
    x[,missing.cols.indices] <- x.missing.imputed
    missing.matrix2 <- is.na(x)
    x[missing.matrix2] <- 0
    
    x.svd <- svd(x)
    lambda.indices <- which(x.svd$d < lambda)
    if(length(lambda.indices) > 0) 
        d.augmented <- c(x.svd$d[-lambda.indices] - lambda, rep(0, length(lambda.indices)))
    else 
        d.augmented <- x.svd$d - lambda 
    
    x.imputed <- x.svd$u %*% diag(d.augmented) %*% x.svd$v
    return(list(x = x.imputed, missing.matrix = missing.matrix))
}    