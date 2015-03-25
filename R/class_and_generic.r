## Setting class
# Set Class for missing Data
setClass(Class="MVData", 
         representation(headname="character", ori.data="matrix", miss.idx="matrix"),
         validity=function(object) {
             #cat("~~~ MVData: inspector ~~~ \n")
             if (!is.matrix(object@ori.data)) stop("The object type of input is not a matrix.")
             if (!is.numeric(object@ori.data)) stop("The data type of input is not numeric.")
             return(TRUE)
         }
)

# Set Class for testing data
setClass("TEData", representation(
            miss.data="matrix",
            miss.rate="numeric"
        ),
        contains="MVData",
        validity=function(object) {
            #cat("~~~ TEData: inspector ~~~ \n")
            if (!is.matrix(object@ori.data) || !is.numeric(object@ori.data)) stop("Complete data is not numeric or not matrix-type!")
            if (!is.matrix(object@miss.data) || !is.numeric(object@miss.data)) stop("Missing data is not numeric or not matrix-type!")
            if (nrow(object@ori.data)!=nrow(object@miss.data) || ncol(object@ori.data)!=ncol(object@miss.data))
                stop("The size of complete matrix is not equal to the size of missing matrix!")
            return(TRUE)            
        }
        # In TEData class, ori.data means the data without any missing values
)

# Set Class for imputed data
setClass(Class="IMPData", 
        representation(comp.data="matrix", imputed.data="matrix", miss.idx="matrix"),
        validity=function(object) {
            #cat("~~~ IMPData: inspector ~~~\n")
            if ( nrow(object@comp.data)!=nrow(object@imputed.data) || ncol(object@comp.data)!=ncol(object@imputed.data) ) {
                stop("The size of the complete matrix is not equal to the size of the imputed matrix!")
            } else if (sum(is.na(object@imputed.data)) != 0) {
                stop("The imputed matrix still contains missing values!")
            } else {}
            return(TRUE)
        }
)


## Setting Generic Function
setGeneric("impute",function(object,method,...) standardGeneric("impute"))
setGeneric("getSize", function(object) standardGeneric("getSize"))
setGeneric("getMissRate", function(object) standardGeneric("getMissRate"))
setGeneric("evaluator", function(object,method,...) standardGeneric("evaluator"))
setGeneric("checkObj", function(object) standardGeneric("checkObj"))

## Setting constructors
# Constructor for MVdata
MVdata<-function(inputfile=NULL) {
    if (inputfile != "" || inputfile != NULL) {
        if (file.exists(inputfile)) {
	        headname<-readLines(inputfile)[1]
	        ori.data<-as.matrix(read.table(inputfile,header=T,row.names=1,sep="\t",skip=1))
            if (!is.matrix(ori.data) && !is.numeric(ori.data)) stop("Input data is not matrix-type or not numeric!")
	        miss.idx<-is.na(ori.data)
			if (sum(miss.idx) == 0) stop("There're no missing entries in this matrix!")
	        new("MVData",headname=headname,ori.data=ori.data,miss.idx=miss.idx)
        } else {
            stop("No such file!")
        }
    } else {
        cat("No input file given, so loading a sample dataset inside EMVA.\n")
        data(lymphoma)
    }
}

# Constructor for TEdata
TEdata<-function(MVdata.obj, miss.rate) {
    headname<-MVdata.obj@headname
    missRow.idx<-which(rowSums(MVdata.obj@miss.idx)!=0)
    # Here, ori.data is a complete matrix without any missing entries
    ori.data<-MVdata.obj@ori.data[-missRow.idx,]
    miss.data<-ori.data
    num.fixxed.row<-300
    if ( num.fixxed.row > (nrow(ori.data)/2) ) {
        warning("The number of rows of this data matrix is less than 300!")
        num.fixxed.row <- nrow(ori.data)/2
    }
    # Randomly select a plenty of rows and keep them unchanged.
    fixxed.row.idx<-sample(seq(nrow(miss.data)),num.fixxed.row)
    fixxed.row.idx<-fixxed.row.idx[order(fixxed.row.idx)]
    temp.data<-miss.data[-fixxed.row.idx,]
    if ( (nrow(temp.data)+length(fixxed.row.idx)) != nrow(ori.data) ) stop("There're some rows missing!\n")
    # Computing the number of entries should be made missing
    miss.data.size<-nrow(miss.data)*ncol(miss.data)
    MV.num<-round(miss.data.size*miss.rate)
    temp.data[sample(seq(nrow(temp.data)*ncol(temp.data)),MV.num)]<-NA
    miss.data[-fixxed.row.idx,]<-temp.data
    miss.idx<-is.na(miss.data)
    new("TEData", headname=headname, ori.data=ori.data, 
        miss.idx=miss.idx, miss.rate=miss.rate, miss.data=miss.data)
}

# Constructor for IMPdata

IMPdata<-function(x.comp,x.imputed,miss.idx){
    new("IMPData", comp.data=x.comp, imputed.data=x.imputed, miss.idx=miss.idx)
}

## Setting Methods
# Show Methods
setMethod("show","MVData",function(object) {
	cat("The class of the object:", class(object), "\n")
    cat("The size of this data:", getSize(object), "\n")
    cat("The missing percentage of this data:", getMissRate(object), "%\n")
})

setMethod("show","TEData",function(object) {
    cat("This is testing data which is generated from a MV data\n")
    cat("The class of the object:", class(object), "\n")
    cat("The size of this data:", getSize(object), "\n")
    cat("The missing percentage of this data:", getMissRate(object), "%\n")
})

setMethod("show", "IMPData", function(object) {
	cat("This is imputed data after imputing a testing data or a MV data\n")
	cat("The class of the object:", class(object), "\n")
	cat("The size of this data:", getSize(object), "\n")
	cat("The missing percentage of this data:", getMissRate(object), "%\n")
})

# Getters
setMethod("getSize","MVData",function(object) dim(object@ori.data))
setMethod("getSize","IMPData",function(object) dim(object@imputed.data))
setMethod("getMissRate","MVData", function(object) sum(object@miss.idx)*100/(nrow(object@ori.data)*ncol(object@ori.data)))
setMethod("getMissRate","IMPData", function(object) sum(is.na(object@imputed.data))*100/(nrow(object@imputed.data)*ncol(object@imputed.data)))


# Check function
setMethod("checkObj","numeric", function(object) {
    if (!is.vector(object)) {
        #cat("This obejct is not a vector!")
        return(FALSE)
    } else {
        #cat("This object is a vector!")
        return(TRUE)
    }
})
setMethod("checkObj","matrix", function(object) {
    if (!is.matrix(object)) {
        #cat("This object is not a matrix!")
        return(FALSE)
    } else {
        #cat("This object is a matrix!")
        return(TRUE)
    }
})
setMethod("checkObj","data.frame", function(object) {
    if (!is.data.frame(object)) {
        #cat("This object is not a data frame!")
        return(FALSE)
    } else {
        #cat("This object is a data frame!")
        return(TRUE)
    }
})

# SetMethod for impute
setMethod("impute","TEData", function(object, method, ...) {
    imputedData<-switch(method,
            zero = ZEROimpute(object@miss.data),
            ravg = RAVGimpute(object@miss.data),
            knn = KNNimpute(object@miss.data, ...),
            sknn = SKNNimpute(object@miss.data, ...),
            iknn = IKNNimpute(object@miss.data, ...),
            svd = SVDimpute(object@miss.data, ...),
            ls = LSimpute(object@miss.data, ...),
            lls = LLSimpute(object@miss.data, ...),
            usr = USRimpute(object@miss.data, ...)
    )
    output<-IMPdata(object@ori.data,imputedData,object@miss.idx)
})

# SetMethod for evaluator
setMethod(f = "evaluator", signature = "IMPData", definition = function(object, method, ...) {
    switch(method,
            nrmse = NRMSEcal(object@comp.data, object@imputed.data, object@miss.idx),
            blci = BLCIcal(object@comp.data, object@imputed.data, ...),
            cpp = CPPcal(object@comp.data, object@imputed.data, ...)
    )
})
