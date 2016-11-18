#' @title mfa
#' @description a brief description
#' @param data: data that need to be analysed
#' @export
#' @return result:...
#' @examples \dontrun{
#' an example}
# data<-read.table("wines.csv",sep=",")

# Function to compute MFA
mfa<-function(data,sets,ncomps=NULL,center=TRUE,scale=TRUE){
  # attention: the first line of data should be the name of variables
  nrow<-nrow(data)-1
  # if SETS is a list of numeric vectors with the position of the active variables in the data table
  if (is.numeric(sets)){
    ncol<-length(sets)
    selected_data<-matrix(numeric(ncol*nrow),ncol=ncol,nrow=nrow)
    # N is the vector to store the ncol of each sub-matrix
    N=numeric(10)
    # write selected data
    j<-1
    for (i in sets){
      selected_data[,j]<-as.numeric(as.character(data[2:(nrow(data)),i]))
      if(i<=7){N[1]=N[1]+1}
      else if(i<=13){N[2]=N[2]+1}
        else if(i<=19){N[3]=N[3]+1}
          else if(i<=24){N[4]=N[4]+1}
            else if(i<=30){N[5]=N[5]+1}
              else if(i<=35){N[6]=N[6]+1}
                else if(i<=39){N[7]=N[7]+1}
                  else if(i<=45){N[8]=N[8]+1}
                    else if(i<=50){N[9]=N[9]+1}
                     else if(i<=58){N[10]=N[10]+1}
      j<-j+1
    }
  }

    # if sets is a list of character vectors with the names of the active variables
    else if(is.character(sets)){
      raw<-NULL
      N<-numeric(10)
      j<-0
      for(i in 2:ncol(data)){
        if(any(sets==data[1,i])){
        raw<-cbind(raw,as.numeric(as.character(data[2:nrow(data),i])))
        j<-j+1
        if(i<=7){N[1]=N[1]+1}
          else if(i<=13){N[2]=N[2]+1}
            else if(i<=19){N[3]=N[3]+1}
              else if(i<=24){N[4]=N[4]+1}
                else if(i<=30){N[5]=N[5]+1}
                  else if(i<=35){N[6]=N[6]+1}
                    else if(i<=39){N[7]=N[7]+1}
                      else if(i<=45){N[8]=N[8]+1}
                        else if(i<=50){N[9]=N[9]+1}
                          else if(i<=58){N[10]=N[10]+1}
        }
      }
      selected_data<-raw
    }
    else{
      stop("sets value is not correct!")
    }
  # ncomps: indicating how many number of components (i.e. factors) are to be extracted
  if(length(ncomps) != 0){
    selected_data<-selected_data[1:ncomps,]
  }
  # scale the data
  selected_data<-scale(selected_data,center,scale)

  # step 1 PCA of Each Data Table
    # centering such that its mean=0
    y<-scale(selected_data,center = TRUE,FALSE)
    X<-NULL
    # normalizing each column such that the sum of the square values of all its elements is equal to 1
    for (i in 1:length(y[1,]))
    {
      sum<-sum(y[,i]^2)
      X<-cbind(X,y[,i]/sqrt(sum))
    }
    N <- N[N>0]

    alpha<-numeric(length(N))
    raw<-NULL
    for (i in 1:length(N))
    {
      if(i==1){
        sum[i]=0
      }
      else {
        sum[i]<-sum(N[1:i-1])
      }
    # SVD of each table
      K<-X[,(sum[i]+1):(sum[i]+N[i])]
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
  c<-as.character(1:length(N))
  for (i in 1:length(N)){
    P_Fi <- length(N)*alpha[i] * X[,(sum[i]+1):(sum[i]+N[i])] %*% Q[(sum[i]+1):(sum[i]+N[i]),]
    attr(P_F,c[i])=P_Fi
  }
  # 4) matrix of loadings (factor loadings) = Q
  result<-list(
    eigen=eigen,
    common_factor_scores=Factor_scores,
    partial_factor_scores=P_F,
    loadings=Q
  )
  class(result) <- "mfa"
  return(result)
}


