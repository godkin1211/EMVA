###########################
# Similarity Calculator
###########################
similarityCal<-function(vec, mat, method="euclidean"){
    if (!checkObj(vec)) {
        if (checkObj(vec)) 
            vec <- as.vector(vec)
        else
            stop("The object type of vec is not valid!")
    }
    
    if (!checkObj(mat)) {
        if (checkObj(mat)) 
            mat <- matrix(mat, nr=1)
        else
            stop("The object type of mat is not valid!")
    }
    
    epsilon <- 1e-8
    methods <- c("euclidean","pearson","cosine")
    switch(match.arg(method,methods),
            euclidean=1/sqrt(rowSums((mat-matrix(vec,nc=length(vec),nr=nrow(mat),byrow=TRUE))^2)+epsilon),
            pearson=apply(mat,1,function(i) abs(cor(i,vec,use="everything",method="pearson"))),
            cosine=apply(mat,1,function(i) abs(sum(i * vec)/sqrt(sum(i^2)*sum(vec^2)))),
            {stop("Pleas choose one method from \"euclidean\", \"pearson\" or \"cosine\"")}
    )
}