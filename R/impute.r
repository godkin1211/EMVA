###########################
# Imputation Algorithms   #
###########################
# Zero
ZEROimpute<-function(x) {
    x[is.na(x)]<-0
    return(x)
}

#Row average
RAVGimpute<-function(x){
    x.imp<-t(apply(x,1,function(i){
        i[is.na(i)]<-mean(i,na.rm=TRUE); i
    }))
    return(x.imp)
}

# KNN
KNNimpute<-function(x, k=15, sim.method="euclidean"){
    miss.idx<-is.na(x)
    x.missRowIdx<-which(rowSums(miss.idx) != 0)
    x.miss<-x[x.missRowIdx,]
    x.comp<-x[-x.missRowIdx,]
    x.imp<-t(apply(x.miss,1,function(i){
        missidx<-which(is.na(i))
        vec.usable<-i[-missidx]
        neighbor.pool<-x.comp[,-missidx]
        similarities<-similarityCal(vec.usable,neighbor.pool,method=sim.method)
        neighborhood<-order(similarities,decreasing=TRUE)[1:k]
        imp.values<-sapply(missidx,function(j) weight.avg<-similarities[neighborhood] %*% x.comp[neighborhood,j] / sum(similarities[neighborhood]))
        i[missidx]<-imp.values
        return(i)
    }))
    x[x.missRowIdx,]<-x.imp
    return(x)
}

# IKNN
IKNNimpute<-function(x, k=10, sim.method="euclidean", iter=3, e=1e-3) {
    cat(sprintf("k = %g, iter= %g, e = %g\n",k,iter,e))
    missIdx<-is.na(x)
    rowNum<-nrow(x)
    miss.RowIdx<-which(rowSums(missIdx) != 0)
    x.ravged<-RAVGimpute(x)
    cat("Row average imputation completed!\n")
    x.miss<-(cbind(1:rowNum, x))[miss.RowIdx,]
    x.complete<-x.ravged
    cat("The size of x.complete:",dim(x.complete),"\n")
    err<-99
    for (r in 1:iter) {
        x.old<-x.complete
        err.old<-err
        cat(sprintf("Start the %g cycle of iknn imputation\n", r))
        x.imputed<-t(apply(x.miss, 1, function(j) {
            rowIdx<-j[1]
            j.origin<-j[-1]
            neighbor.pool<-x.complete[-rowIdx,]
            target<-x.complete[rowIdx,]
            dist.list<-similarityCal(target, neighbor.pool, method=sim.method)
            neighborsIdx<-order(dist.list,decreasing=T)[1:k]
            missColIdx<-which(is.na(j.origin))
            estimation<-sapply(missColIdx, function(h){
                weight<-dist.list[neighborsIdx]
                weightedAvg<-weight %*% neighbor.pool[neighborsIdx,h]/sum(weight)
                return(weightedAvg)
            })
            j.origin[missColIdx]<-estimation
            return(j.origin)
        }))
        cat(sprintf("The %g cycle imputation has been completed\n",r))
        x.complete[miss.RowIdx,]<-x.imputed
        err<-sum((x.old[missIdx]-x.complete[missIdx])^2)
        cat("err:",err,"\n")
        if (r>1) if ((err.old/err < 4) | (err < e)) break 
    }
    x<-x.complete
    return(x)
}

# SKNN
SKNNimpute<-function(x, k=10, sim.method="euclidean"){
    rowNum<-nrow(x)
    colNum<-ncol(x)
    missIdx<-is.na(x)
    miss.RowIdx<-which(rowSums(missIdx)!=0)
    x.completeRows<-x[-miss.RowIdx,]
    x.missingRows<-x[miss.RowIdx,]
    miss.list<-order(rowSums(is.na(x.missingRows)))
    for (i in seq(nrow(x.missingRows))) {
        target<-x.missingRows[miss.list[i],]
        dist.list<-similarityCal(target, x.completeRows, sim.method)
        neighborsIdx<-order(abs(dist.list),decreasing=T)[1:k]
        missColIdx<-which(is.na(target))
        estimation<-sapply(missColIdx, function(j){
            weight<-dist.list[neighborsIdx]
            weightedAvg<-weight %*% x.completeRows[neighborsIdx,j] / sum(weight)
            return(weightedAvg)
        })
        target[missColIdx]<-estimation
        x.missingRows[miss.list[i],]<-target
        x.completeRows<-rbind(x.completeRows, target)
    }
    x[miss.RowIdx,]<-x.missingRows
    return(x)
}

# SVD
SVDimpute<-function(x, k=15, iters=10, sim.method="euclidean"){
    missing.matrix <- is.na(x)
    numMissing <- sum(missing.matrix)
    print(paste("imputing on", numMissing, "missing values with matrix size", nrow(x)*ncol(x), sep=" "))
    if(numMissing == 0) return (x)
    miss.colIdx <- which(apply(missing.matrix, 2, function(i) any(i) )) 
    x.missing <- (rbind(1:ncol(x), x))[,miss.colIdx]
    x.imputed <- apply(x.missing, 2, function(j) {
        colIdx <- j[1]
        j.original <- j[-1]
        missing.rows <- which(missing.matrix[,colIdx])
        if(length(missing.rows) == nrow(x))
            warning( paste("Column",colIdx,"is completely missing",sep=" ") )
        j.original[missing.rows] <- mean(j.original[-missing.rows])
        j.original
    })
    x[,miss.colIdx] <- x.imputed
    missing.matrix2 <- is.na(x)
    x[missing.matrix2] <- 0
    for(i in 1:iters) {
        print(paste("Running iteration", i, sep=" "))
        x.svd <- svd(x, nu=k, nv=k)
        x.svd <- x.svd$u %*% diag(x.svd$d[1:k],nrow=k,ncol=k) %*% t(x.svd$v)
        x[missing.matrix] <- x.svd[missing.matrix]
    }
    return(x)
}

# LSimpute
LSimpute<-function(x, k=15, sim.method="pearson", e=1e-6){
    rowNum<-nrow(x)
    colNum<-ncol(x)
    missIdx<-is.na(x)
    miss.RowIdx<-which(rowSums(missIdx)!=0)
    x.completeRows<-x[-miss.RowIdx,]
    x.missingRows<-x[miss.RowIdx,]
    x.imputed<-t(apply(x.missingRows,1,function(i){
        missColIdx<-which(is.na(i))
        dist.list<-similarityCal(i[-missColIdx], x.completeRows[,-missColIdx], method=sim.method)
        neighborsIdx<-order(dist.list,decreasing=T)[1:k]
        neighbor.sim<-dist.list[neighborsIdx]
        weight<-((neighbor.sim^2)/(1-(neighbor.sim)^2+e))^2
        fit.list<-lapply(neighborsIdx, function(m){
            fit<-lm(i[-missColIdx] ~ x.completeRows[m,-missColIdx])
            reg<-c(fit$coefficients[[1]],fit$coefficients[[2]])
            return(reg)
        })
        estimation<-sapply(missColIdx,function(j){
            regression<-rep(0,k)
            for (n in seq(k)) {
                a0<-fit.list[n][[1]][1]
                b0<-fit.list[n][[1]][2]
                regression[n]<-a0+b0*x.completeRows[neighborsIdx[n],j]
            }
            weightedAvg<-(weight %*% regression)/sum(weight)
            return(weightedAvg)
        })
        i[missColIdx]<-estimation
        return(i)
    }))
    x[miss.RowIdx,]<-x.imputed
    return(x)
}

# LLS
LLSimpute<-function(x, k=50, sim.method="euclidean"){
    missidx<-is.na(x)
    missRowIdx<-which(rowSums(missidx) != 0)
    x.completePart<-x[-missRowIdx,]
    x.missingPart<-x[missRowIdx,]
    imputed<-t(apply(x.missingPart,1,function(i){
        missColIdx<-which(is.na(i))
        dist.list<-similarityCal(i[-missColIdx],x.completePart[,-missColIdx],sim.method)
        neighborIdx<-order(dist.list,decreasing=T)[1:k]
        A<-x.completePart[neighborIdx,-missColIdx,drop=FALSE]
        b<-x.completePart[neighborIdx,missColIdx,drop=FALSE]
        Cp<-ginv(t(A))
        w<-i[-missColIdx,drop=FALSE]
        X<-Cp %*% w
        ans<-t(b) %*% X
        i[missColIdx]<-ans
        return(i)
    }))
    x[missRowIdx,]<-imputed
    cat("LLS imputation completed!\n")
    return(x)
}

# User-defined
#USRimpute<-function(usr.FUN,x, ...) { usr.FUN(x, ...) }