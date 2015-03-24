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
    if (! "apcluster" %in% sessionInfo()$otherPkgs) library(aplcuster)
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
BLCIcal<-function(mat.com, mat.imp, resp.type="Pattern discovery", eg=10, nperms=100){
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
    genenames<-paste("GENE",as.character(1:nrow(mat.com)),sep="")
    geneid<-as.character(1:nrow(mat.com))
    ans.data<-list(x=mat.com, eigengene.number=eg, geneid=geneid, genenames=genenames)
    imp.data<-list(x=mat.imp, eigengene.number=eg, geneid=geneid, genenames=genenames)
    ans.samr.obj<-samr(ans.data, resp.type=resp.type, nperms=nperms)
    imp.samr.obj<-samr(imp.data, resp.type=resp.type, nperms=nperms)
    ans.delta<-samr.compute.delta.table(ans.samr.obj)
    imp.delta<-samr.compute.delta.table(imp.samr.obj)
    ans.siggene.table<-samr.compute.siggenes.table(ans.samr.obj, del=0, ans.data, ans.delta, all.genes=T)
    imp.siggene.table<-samr.compute.siggenes.table(imp.samr.obj, del=0, imp.data, imp.delta, all.genes=T)
    ans.totalgenes<-rbind(ans.siggene.table$genes.lo, ans.siggene.table$genes.up)
    imp.totalgenes<-rbind(imp.siggene.table$genes.lo, imp.siggene.table$genes.up)
    ans.siggene.judge<-as.numeric(ans.totalgenes[,7])
    imp.siggene.judge<-as.numeric(imp.totalgenes[,7])
    ans.siggenes<-ans.totalgenes[ans.siggene.judge<10,,drop=F]
    imp.siggenes<-imp.totalgenes[imp.siggene.judge<10,,drop=F]
    ans.siggene.list<-ans.siggenes[,2]
    imp.siggene.list<-imp.siggenes[,2]
    ans.nonsiggenes<-ans.totalgenes[ans.siggene.judge >= 10,,drop=F]
    imp.nonsiggenes<-imp.totalgenes[imp.siggene.judge >= 10,,drop=F]
    ans.nonsiggene.list<-ans.nonsiggenes[,2]
    imp.nonsiggene.list<-imp.nonsiggenes[,2]
    ans.imp.siggenes.int<-intersect(ans.siggene.list,imp.siggene.list)
    ans.imp.nonsiggenes.int<-intersect(ans.nonsiggene.list,imp.nonsiggene.list)
    if (length(ans.imp.siggenes.int) != 0 & length(ans.siggene.list) != 0)
        blci<-length(ans.imp.siggenes.int)/length(ans.siggene.list) + length(ans.imp.nonsiggenes.int)/length(ans.nonsiggene.list) - 1
    else
        blci<-length(ans.imp.nonsiggenes.int)/length(ans.nonsiggene.list) - 1

    return(blci)
}
