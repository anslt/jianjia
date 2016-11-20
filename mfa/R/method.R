summary_eigen<-function(x){
  if(class(x)!="mfa"){
    stop("\n'object' need to be mfa")
  }
  eigens<-x$eigen
  singular<-sqrt(eigens)
  cumulative<-cumsum(eigens)
  interia<-eigens/cumulative[length(cumulative)]
  i_cunmulative<-cumsum(interia)
  result<-data.frame("singular_values"=singular,"eigenvalue"=eigens,
                     "cumulative"=cumulative,"interia"=interia,
                     "cumulative"=i_cunmulative)
  result
}

#Three types-observation,variable and table
contribution<-function(x,type){
  if(class(x)!="mfa"){
    stop("\n'object' need to be mfa")
  }
  if(type=="observation"){
    summation=colSums(x$common_factor_scores^2  %*% diag(x$mass))
    percent=apply(x$common_factor_scores,1,function(y) y^2*(x$mass)/summation)
    colnames(percent)<-paste0(rep("dim",length(percent[1,])),seq(length(percent[1,])))
    rownames(percent)<-paste0(rep("ob",length(percent[,1])),seq(length(percent[,1])))
    round(percent,3)
  }
  else if(type=="variable"){
    percent=diag(rep(x$weight,diff(x$divide))) %*% x$loadings^2  
    colnames(percent)<-paste0(rep("dim",length(percent[1,])),seq(length(percent[1,])))
    rownames(percent)<-paste0(rep("vr",length(percent[,1])),seq(length(percent[,1])))
    percent
  }
  else if(type=="table"){
    tmp=diag(rep(x$weight,diff(x$divide))) %*% x$loadings^2
    percent=matrix(ncol=dim(x$loadings)[2])
    for(num in seq(length(x$divide)-1)){
      table<-tmp[x$divide[num]:(x$divide[num+1]-1),]
      table<-colSums(table)
      if(any(is.na(percent))){
        percent<-table
      }
      else{
        percent<-rbind(percent,table)
      }
    }
    colnames(percent)<-paste0(rep("dim",length(percent[1,])),seq(length(percent[1,])))
    rownames(percent)<-paste0(rep("tb",length(percent[,1])),seq(length(percent[,1])))
    round(percent,3)
  }
  else{
    print("The type is Error")
  }
}
#To find the sum of the diagonal entried of the matrix
tr<-function(matrix){
  if(dim(matrix)[1]!=dim(matrix)[2]){
    stop("\nThe matrix need to be a squared matrix")
  }
  sum(diag(matrix))
}

RV<-function(table1,table2){
  numerator=tr((table1 %*% t(table1)) %*% (table2 %*% t(table2)))
  denominator=sqrt(tr(table1 %*% t(table1) %*% table1 %*% t(table1))*
                     tr(table2 %*% t(table2) %*% table2 %*% t(table2)))
  numerator/denominator
}

RV_table<-function(dataset,sets){
  N<-length(sets)
  result<-matrix(ncol=N,nrow=N)
  store<-list()
  for(x in seq(N)){
    store[[x]]=as.matrix(data[,sets[[x]]])
  }
  for(i in seq(N)){
    for(j in seq(N)){
      result[i,j]<-RV(store[[i]],store[[j]])
    }
  }
  colnames(result)<-paste0(rep("tb",length(N)),seq(N))
  rownames(result)<-paste0(rep("tb",length(N)),seq(N))
  result
}

#We give the original normailized data and the number of block 
#ro construct the table
LG<-function(table1,table2){
  a1<-svd(table1)$d[1]
  a2<-svd(table2)$d[1]
  tr((table1 %*% t(table1)) %*% (table2 %*% t(table2)))*(1/a1^2)*(1/a2^2)
}

