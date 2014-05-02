## Setting class
# Set Class for missing Data
setClass("MVData", representation(
            headname="character",
            ori.data="matrix",
            miss.idx="matrix"
        )
)

# Set Class for testing data
setClass("TEData", representation(
            miss.data="matrix",
            miss.rate="numeric"
        ),
        contains="MVData"
)

# Set Class for imputed data
setClass("IMPData",representation(
            ori.data="matrix",
            imputed.data="matrix",
            miss.idx="matrix"
        )
)


## Setting Generic Function
setGeneric("impute",function(x,method,...) standardGeneric("impute"))
setGeneric("getSize", function(x) standardGeneric("getSize"))
setGeneric("getMissRate", function(x) standardGeneric("getMissRate"))
setGeneric("evaluator", function(x,method,...) standardGeneric("evaluator"))

## Setting constructors
# Constructor for MVdata
MVdata<-function(inputfile=NULL) {
    if (inputfile != "" || inputfile != NULL) {
	    headname<-readLines(inputfile)[1]
	    ori.data<-as.matrix(read.table(inputfile,header=T,row.names=1,sep="\t",skip=1))
	    miss.idx<-is.na(ori.data)
	    new("MVData",headname=headname,ori.data=ori.data,miss.idx=miss.idx)
    } else {
        data(lymphoma)
    }
}

# Constructor for TEdata
TEdata<-function(MVdata.obj, miss.rate) {
    headname<-MVdata.obj@headname
    missRow.idx<-which(rowSums(MVdata.obj@miss.idx)!=0)
    ori.data<-MVdata.obj@ori.data[-missRow.idx,]
    miss.data<-ori.data
    num.fixxed.row<-300
    fixxed.row.idx<-sample(seq(nrow(miss.data)),num.fixxed.row)
    fixxed.row.idx<-fixxed.row.idx[order(fixxed.row.idx)]
    temp.data<-miss.data[-fixxed.row.idx,]
    miss.data.size<-nrow(miss.data)*ncol(miss.data)
    MV.num<-round(miss.data.size*miss.rate)
    temp.data[sample(seq(nrow(temp.data)*ncol(temp.data)),MV.num)]<-NA
    miss.data[-fixxed.row.idx,]<-temp.data
    miss.idx<-is.na(miss.data)
    new("TEData",headname=headname, ori.data=ori.data, miss.idx=miss.idx, miss.rate=miss.rate, miss.data=miss.data)
}

IMPdata<-function(x.ori,x.imputed,miss.idx){
    new("IMPData",ori.data=x.ori,imputed.data=x.imputed,miss.idx=miss.idx)
}

## Setting Methods
# Show Methods
setMethod("show","MVData",function(object) {
	cat("The class of the object:",class(object),"\n")
    cat("The size of this data:", getSize(object),"\n")
    cat("The missing percentage of this data:",getMissRate(object),"%\n")
})
setMethod("show","TEData",function(object) {
    cat("This is testing data generated with a MV data\n")
    cat("The class of the object:",class(object),"\n")
    cat("The size of this data:", getSize(object),"\n")
    cat("The missing percentage of this data:",getMissRate(object),"%\n")
})

# Getters
setMethod("getSize","MVData",function(x) dim(x@ori.data))
setMethod("getMissRate","MVData", function(x) sum(x@miss.idx)*100/(nrow(x@ori.data)*ncol(x@ori.data)))

# Setters

# SetMethod for impute
setMethod("impute","TEData", function(x, method, ...) {
    imputedData<-switch(method,
            zero = ZEROimpute(x@miss.data),
            ravg = RAVGimpute(x@miss.data),
            knn = KNNimpute(x@miss.data, ...),
            sknn = SKNNimpute(x@miss.data, ...),
            iknn = IKNNimpute(x@miss.data, ...),
            svd = SVDimpute(x@miss.data, ...),
            ls = LSimpute(x@miss.data, ...),
            lls = LLSimpute(x@miss.data, ...),
            usr = USRimpute(x@miss.data, ...)
    )
    output<-IMPdata(x@ori.data,imputedData,x@miss.idx)
})

# SetMethod for evaluator
setMethod(f = "evaluator", signature = "IMPData", definition = function(x, method, ...) {
    switch(method,
            nrmse = NRMSEcal(x@ori.data, x@imputed.data, x@miss.idx),
            blci = BLCIcal(x@ori.data, x@imputed.data,...),
            cpp = CPPcal(x@ori.data, x@imputed.data, k=10,...)
    )
})
