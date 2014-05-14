###########################
# Similarity Calculator
###########################
similarityCal<-function(vec, mat, method="euclidean"){
    epsilon <- 1e-8
    methods <- c("euclidean","pearson","cosine")
    switch(match.arg(method,methods),
           euclidean=1/sqrt(rowSums((mat-matrix(vec,nc=length(vec),nr=nrow(mat),byrow=TRUE))^2)+epsilon),
           pearson=apply(mat,1,function(i) abs(cor(i,vec,use="everything",method="pearson"))),
           cosine=apply(mat,1,function(i) abs(sum(i * vec)/sqrt(sum(i^2)*sum(vec^2))))
    )
}