##########################
# Performance evaluator
##########################
# NRMSE
NRMSEcal<-function(mat.com, mat.imp, miss.idx) {
    # Check the content of each argument
    if (sum(is.na(mat.com))!=0) stop("There exists NAs in ans!\n")
    if (sum(is.na(mat.imp))!=0) stop("There exists NAs in imputed data!\n")
    if (!is.logical(miss.idx)) stop("Problems in miss.idx!\n")
    if (length(mat.com[miss.idx])!=length(mat.imp[miss.idx])) stop("Length inequivalent!\n")
    if (!checkObj(mat.com)) {
        if (checkObj(mat.com))
            mat.com<-as.matrix(mat.com)
        else
            stop("Invalid object type on mat.com!")
    }
    if (!checkObj(mat.imp)) {
        if (checkObj)
            mat.com<-as.matrix(mat.com)
        else
            stop("Invalid object type on mat.com!")
    }
    nrmse<-sqrt(mean((mat.com[miss.idx]-mat.imp[miss.idx])^2)/var(mat.com[miss.idx]))
    return(nrmse)
}

# CPP
CPPcal<-function(mat.com,mat.imp){
    if (sum(is.na(mat.com))!=0) stop("There exists NAs in ans!\n")
    if (sum(is.na(mat.imp))!=0) stop("There exists NAs in imputed data!\n")
    if (!checkObj(mat.com)) {
        if (checkObj(mat.com))
            mat.com<-as.matrix(mat.com)
        else
            stop("Invalid object type on mat.com!")
    }
    if (!checkObj(mat.imp)) {
        if (checkObj)
            mat.com<-as.matrix(mat.com)
        else
            stop("Invalid object type on mat.com!")
    }
    # Using Affinity propagation (AP) clustering to determine the selection of k in kmeans
    if (! "apcluster" %in% sessionInfo()$otherPkgs) library(apcluster)
    cluster.result <- apcluster(negDistMat(r=2), mat.com)
    k <- length(cluster.result@clusters)
    ans.clustered<-kmeans(mat.com,k)
    ans.clust<-ans.clustered$cluster
    imp.clustered<-kmeans(mat.imp,k)
    imp.clust<-imp.clustered$cluster
    cpptable<-matrix(nr=k,nc=k)
    rownames(cpptable)<-sprintf("C.ans.%i",1:k)
    colnames(cpptable)<-sprintf("C.imp.%i",1:k)
    for (i in 1:k) {
        ans.genes.inthiscluster<-which(ans.clust==i)
        for (j in 1:k) {
            imp.genes.inthiscluster<-which(imp.clust==j)
            common.member<-length(intersect(ans.genes.inthiscluster,imp.genes.inthiscluster))
            cpptable[i,j]<-common.member
        }
    }
    cpp<-sum(apply(cpptable,1,max))/sum(cpptable)
    return(cpp)
}

# BLCI
BLCIcal<-function(mat.com, mat.imp, design.matrix=NULL, contrasts=NULL){
    if (! "limma" %in% sessionInfo()$otherPkgs) library(limma)
    if (sum(is.na(mat.imp))) stop("There exists NAs in imputed matrix!")
    if (sum(is.na(mat.imp))!=0) stop("There exists NAs in imputed data!\n")
    if (!checkObj(mat.com)) {
        if (checkObj(mat.com))
            mat.com<-as.matrix(mat.com)
        else
            stop("Invalid object type on mat.com!")
    }
    if (!checkObj(mat.imp)) {
        if (checkObj)
            mat.com<-as.matrix(mat.com)
        else
            stop("Invalid object type on mat.com!")
    }
    DE.table.com <- DE.table.imp <-NULL
    if (is.null(design.matrix)) {
        cat("\nBecause design.matrix does not be set, 
            EMVA will simply assume that this is two-group experiment.\n")
        colNum <- ncol(mat.com)
        if (colNum %% 2 != 0) stop("The number of experiments/conditions is an odd number!")
        Group <- factor(rep(c("WT","Mu"), each=colNum/2), levels=c("WT","Mu"))
        design.matrix <- model.matrix(~0+Group)
        colnames(design.matrix) <- levels(Group)
        fit.com <- lmFit(mat.com, design.matrix)
        fit.imp <- lmFit(mat.imp, design.matrix)
        cont.matrix <- makeContrasts(Mu-WT, levels=design.matrix)
    } else {
        fit.com <- lmFit(mat.com, design.matrix)
        fit.imp <- lmFit(mat.imp, design.matrix)
        cont.matrix <- contrasts
    }
    
    fit2.com <- contrasts.fit(fit.com, cont.matrix)
    fit2.com <- eBayes(fit2.com, 0.01)
    DE.table.com <- topTable(fit2.com, adjust="fdr", sort.by="B", number=nrow(mat.com))
    
    fit2.imp <- contrasts.fit(fit.imp, cont.matrix)
    fit2.imp <- eBayes(fit2.imp, 0.01)
    DE.table.imp <- topTable(fit2.imp, adjust="fdr", sort.by="B", number=nrow(mat.imp))

    num.siggenes.com <- sum(DE.table.com$adj.P.Val <= 0.05)
    num.nonsiggenes.com <- nrow(mat.com) - num.siggenes.com
    num.siggenes.imp <- sum(DE.table.imp$adj.P.Val <= 0.05)
    num.nonsiggenes.imp <- nrow(mat.imp) - num.siggenes.imp
    
    blci <- num.siggenes.com/num.siggenes.imp + num.nonsiggenes.com/num.nonsiggenes.imp - 1
    
    #if (length(ans.imp.siggenes.int) != 0 & length(ans.siggene.list) != 0)
    #    blci<-length(ans.imp.siggenes.int)/length(ans.siggene.list) + length(ans.imp.nonsiggenes.int)/length(ans.nonsiggene.list) - 1
    #else
    #    blci<-length(ans.imp.nonsiggenes.int)/length(ans.nonsiggene.list) - 1

    return(blci)
}
