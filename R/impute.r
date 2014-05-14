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
    # Impute the row that are all-missing or almost-all-missing with column-average.
    if (length(which(rowSums(is.na(x))==ncol(x)))) for (i in which(rowSums(is.na(x))==ncol(x))) x[i,] <- colMeans(x,na.rm=T)
    
    if (length(which(rowSums(is.na(x))==ncol(x)-1))) {
        for(i in which(rowSums(is.na(x))==ncol(x)-1)) {
            missCol <- which(is.na(x[i,]))
            x[i,missCol] <- colMeans(x,na.rm=T)[missCol]
        }
    }
    
    miss.idx <- is.na(x)
    x.missRowIdx <- which(rowSums(miss.idx) != 0)
    x.miss <- x[x.missRowIdx,]
    x.comp <- x[-x.missRowIdx,]
    x.imp <- t(apply(x.miss,1,function(i){
        missidx <- which(is.na(i))
        vec.usable <- i[-missidx]
        neighbor.pool <- x.comp[,-missidx,drop=F]
        similarities <- similarityCal(vec.usable,neighbor.pool,method=sim.method)
        neighborhood <- order(similarities,decreasing=TRUE)[1:k]
        imp.values <- sapply(missidx,function(j) {
            weight.avg<-similarities[neighborhood] %*% x.comp[neighborhood,j,drop=F] / sum(similarities[neighborhood])
            weight.avg
        })
        i[missidx] <- imp.values
        return(i)
    }))
    x[x.missRowIdx,] <- x.imp
    return(x)
}

# IKNN
IKNNimpute<-function(x, k=10, sim.method="euclidean", iter=3, e=1e-3) {
    missIdx <- is.na(x)
    rowNum <- nrow(x)
    colNum <- ncol(x)
    # Impute the row that are all-missing or almost-all-missing with column-average.
    if (length(which(rowSums(is.na(x))==colNum))) for (i in which(rowSums(is.na(x))==colNum)) x[i,] <- colMeans(x,na.rm=T)
    
    if (length(which(rowSums(is.na(x))==colNum-1))) {
        for(i in which(rowSums(is.na(x))==colNum-1)) {
            missCol <- which(is.na(x[i,]))
            x[i,missCol] <- colMeans(x,na.rm=T)[missCol]
        }
    }
    
    missIdx<-is.na(x)
    miss.RowIdx<-which(rowSums(missIdx) != 0)
    x.ravged<-RAVGimpute(x)
    x.miss<-(cbind(1:rowNum, x))[miss.RowIdx,]
    x.complete<-x.ravged
    err<-99
    for (r in 1:iter) {
        x.old<-x.complete
        err.old<-err
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
        x.complete[miss.RowIdx,]<-x.imputed
        err<-sum((x.old[missIdx]-x.complete[missIdx])^2)
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
    # Impute the row that are all-missing or almost-all-missing with column-average.
    if (length(which(rowSums(missIdx)==colNum))) for (i in which(rowSums(missIdx)==colNum)) x[i,] <- colMeans(x,na.rm=T)
    
    if (length(which(rowSums(missIdx)==colNum-1))) {
        for(i in which(rowSums(missIdx)==colNum-1)) {
            missCol <- which(is.na(x[i,]))
            x[i,missCol] <- colMeans(x,na.rm=T)[missCol]
        }
    }
    missIdx<-is.na(x)
    miss.RowIdx<-which(rowSums(missIdx)!=0)
    x.completeRows<-x[-miss.RowIdx,]
    x.missingRows<-x[miss.RowIdx,]
    miss.list<-order(rowSums(is.na(x.missingRows)))
    for (i in seq(nrow(x.missingRows))) {
        target<-x.missingRows[miss.list[i],]
        missColIdx<-which(is.na(target))
        dist.list<-similarityCal(target[-missColIdx], x.completeRows[,-missColIdx], sim.method)
        neighborsIdx<-order(abs(dist.list),decreasing=T)[1:k]
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
SVDimpute<-function(x, k=15, iters=10){
    if (ncol(x) < k) k<-ncol(x)
    missIdx <- is.na(x)
	# First, replace MVs with row-average/col-average
	x.imputed<-RAVGimpute(x)
    x <- x.imputed
	# Second, using singular value decompostion to re-estimation original MVs.
    for(i in 1:iters) {
        x.svd <- svd(x, nu=k, nv=k)
        x.svd <- x.svd$u %*% diag(x.svd$d[1:k],nrow=k,ncol=k) %*% t(x.svd$v)
        x[missIdx] <- x.svd[missIdx]
    }
    return(x)
}

# LSimpute
LSimpute<-function(x, k=15, sim.method="euclidean", e=1e-6){
    rowNum<-nrow(x)
    colNum<-ncol(x)
    missIdx<-is.na(x)
    # Impute the row that are all-missing or almost-all-missing with column-average.
    if (length(which(rowSums(missIdx)==colNum))) for (i in which(rowSums(missIdx)==colNum)) x[i,] <- colMeans(x,na.rm=T)
    
    if (length(which(rowSums(missIdx)==colNum-1))) {
        for(i in which(rowSums(missIdx)==colNum-1)) {
            missCol <- which(is.na(x[i,]))
            x[i,missCol] <- colMeans(x,na.rm=T)[missCol]
        }
    }

    miss.RowIdx<-which(rowSums(is.na(x))!=0)
    x.completeRows<-x[-miss.RowIdx,]
    x.missingRows<-x[miss.RowIdx,]
    x.imputed<-t(apply(x.missingRows,1,function(i){
        missColIdx<-which(is.na(i))
        dist.list<-similarityCal(i[-missColIdx], x.completeRows[,-missColIdx,drop=F], method=sim.method)
        neighborsIdx<-order(dist.list,decreasing=T)[1:k]
        neighbor.sim<-dist.list[neighborsIdx]
        weight<-((neighbor.sim^2)/(1-(neighbor.sim)^2+e))^2
        fit.list<-lapply(neighborsIdx, function(m){
            fit<-lm(i[-missColIdx] ~ x.completeRows[m,-missColIdx])
            if (is.na(fit$coefficients[[2]])) fit$coefficients[[2]] <- 0
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
    colNum <- ncol(x)
    missidx <- is.na(x)
    # Impute the row that are all-missing or almost-all-missing with column-average.
    if (length(which(rowSums(missidx)==colNum))) for (i in which(rowSums(missidx)==colNum)) x[i,] <- colMeans(x,na.rm=T)
    
    if (length(which(rowSums(missidx)==colNum-1))) {
        for(i in which(rowSums(missidx)==colNum-1)) {
            missCol <- which(is.na(x[i,]))
            x[i,missCol] <- colMeans(x,na.rm=T)[missCol]
        }
    }
    missidx <- is.na(x)
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
    return(x)
}

# User-defined
#USRimpute<-function(x, ...) {  }
