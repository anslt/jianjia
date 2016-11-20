
# data<-read.table("wines.csv",sep=",") 

# Function to compute MFA
mfa1<-function(data,sets,ncomps=NULL,center=TRUE,scale=TRUE){

  # attention: the first line of data should be the name of variables
  nrow=nrow(data)
  ncol=0
  for(ele in sets){
    ncol<-ncol+length(ele)
  }
  # if SETS is a list of numeric vectors with the position of the active variables in the data table
  # if sets is a list of character vectors with the names of the active variables
  if (is.numeric(unlist(sets)) | is.character(unlist(sets))){
    selected_data<-matrix(numeric(ncol*nrow),ncol=ncol,nrow=nrow)
    # N is the vector to store the ncol of each sub-matrix
    N=numeric(length(sets))
    # write selected data
    k<-1
    i<-1
    for (v in sets){
      for(j in v){
        selected_data[,k]<-data[,j]
        N[i]<-N[i]+1
        k<-k+1
      }
      i<-i+1
    }
  }
    else{
      stop("\nsets value is not correct!")
    }

  # ncomps: indicating how many number of components (i.e. factors) are to be extracted
  if(length(ncomps) != 0){
    selected_data<-selected_data[1:ncomps,]
    nrow<-ncomps
  }
  # scale the data
  if(is.numeric(scale) & length(scale)==1){
    y<-scale(selected_data,center = TRUE,FALSE)
    X<-NULL
    # normalizing each column such that the sum of the square values of all its elements is equal to 1
    for (i in 1:length(y[1,]))
    {
      sum1<-sum(y[,i]^2)
      X<-cbind(X,y[,i]*scale/sqrt(sum1))
    }
  }
  else{
    X<-scale(selected_data,center,scale)
  }
  # step 1 PCA of Each Data Table
    # centering such that its mean=0 
    summation<-numeric(length(N))
    alpha<-numeric(length(N))
    raw<-NULL
    for (i in 1:length(N))
    {
      if(i==1){
        summation[i]=0
      }
      else {
        summation[i]<-sum(N[1:i-1])
      }
    # SVD of each table
      K<-X[,(summation[i]+1):(summation[i]+N[i])]
    # the weight of a table is obtained from the first singular value of its PCA 
      alpha[i]<-1/(svd(K)$d[1]^2)
      raw<-c(raw,rep(alpha[i],N[i]))
    }
  # Step 2
  A<-diag(raw)
  M<-diag(1/nrow,nrow,nrow)
  X_new <- sqrt(M) %*% X %*% sqrt(A)
  P <- svd(X_new)$u / sqrt(1/nrow)
  Y <- 1/sqrt(raw)
  Y <- diag(Y)
  Q <- Y %*% svd(X_new)$v  
  
  # 1) vector containing the eigenvalues
  d <- t(P) %*% M %*% X %*% A %*% Q
  eigen <- diag(d)[diag(d)>1E-05]^2
  # 2) matrix of common factor scores
  Factor_scores <- P %*% d
  # 3) partial factor scores
  P_F<-array()
  dimension<-as.character(1:length(N))
  for (i in 1:length(N)){
    P_Fi <- length(N)*alpha[i] * X[,(summation[i]+1):(summation[i]+N[i])] %*% Q[(summation[i]+1):(summation[i]+N[i]),]
    attr(P_F,dimension[i])=P_Fi
  }
  # 4) matrix of loadings (factor loadings) = Q
  result<-list(
    eigen=eigen,
    common_factor_scores=Factor_scores,
    partial_factor_scores=P_F,
    loadings=Q,
    divide=cumsum(c(c(0),N))+1,
    weight=alpha,
    mass=rep(1/nrow,nrow),
    scale_x=X
  )
  class(result) <- "mfa"
  return(result)
}



