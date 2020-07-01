#------------------------------------------#
# Subcodes (general)   by Jun Song
#------------------------------------------#
# : Generally used R codes in my research
##################################################################
# 1. gam : estimating gamma of the gaussian kernel
# 2. center : center the multivariate x
# 3. matpower : Moore-Penrose inverse ; threshold for eigenvalue = 10^(-15)
# 4. matpower2 : Moore-Penrose inverse with Tychonoff regularization; threshold for eigenvalue = 10^(-15)
# 5. x.ip2 : give Gram matrix with given kernel
# 6. sym : make a symmetric matrix
# 7. spearman : spearman correlation
# 8. normmat : give n x n  norm matrix (ij)= ||X_i - X_j ||
# 9. makematrix : convert list to matrix (this is for functional data with balanced timepoints)
# 10. makelist : convert matrix to list (reverse of makematrix)
# 11. makedataframe : convert list(functionla data) into dataframe
# 12. optionout : just print the option (just for self-test)
# 13. choose.tset : Let A , B two sets of time points. For each time point in A, find the element in B that is closest to that time point.
# 14. cov.ftn : Matern covariance function
# 15. matern.ftn : Simulate Matern process
# 16. bdiag : make block diagonal matrix (for estimation with vector valued functions)
# 17. random.beta : generate random function with specified RKHS
# 18. diag.mat: block diagonal matrix of same matri
# 19. plotpred : plot with labels
##################################################################


##################################################################
# gam(normX) estimating width-parameter of kernel function -----
##################################################################
# input : normX : n x n matrix such that 
#       (normX)_ij = ||x_i - x_j||^2
# value : estimated gamma 
##################################################################
gam <- function(normx){
  n <- nrow(normx)
  avg <- 0 
  avg <- sum(sqrt(normx))/2 + sum(diag(sqrt(normx)))/2
  avg <- avg/(n*(n-1)/2)
  gamma <- 1/(avg)^2
  return(gamma)}

##############################################################
#      function 9: center X (n by p matrix)------        
##############################################################
center <- function(x){
  return(t(t(x)-apply(x,2,mean)))}

##################################################################
# matpower(x,a) : power of matrix (moore penrose type)------
##################################################################
# input : x : matrix
# value : x^a
##################################################################
matpower <- function(a,alpha){ # only for symmetric matrix
  tmp <- eigen(a,sym=TRUE)
  eval <- tmp$values
  evec <- tmp$vectors
  # remove vectors if the eigenvalue <0 or smaller then 10^(-15)
  m <- sum(eval > 10^(-15))
  return(evec[,1:m]%*%diag((eval[1:m])^alpha)%*%
           t(evec[,1:m]))}

##################################################################
# matpower2(x,a,epsilon) : power of matrix(regularized moore-penrose)------
##################################################################
# input : x : matrix
# value : x^a
##################################################################
matpower2 <- function(a,alpha,epsilon){ # only for symmetric matrix
  n <- nrow(a)
  scale <- eigen(a,sym=TRUE)$val[1]
  #   scale = 1
  tmp <- eigen(a+scale*epsilon*diag(n),sym=TRUE)
  eval <- tmp$values
  evec <- tmp$vectors
  # remove vectors if the eigenvalue <0 or smaller then 10^(-15)
  m <- sum(eval > 10^(-15))
  return(evec[,1:m]%*%diag((eval[1:m])^alpha)%*%
           t(evec[,1:m]))}
##################################################################
# x.ip2 :  inner product matrix when x is vector-valued function------
##################################################################
# Inner product matrix of X
# output : n x n matrix,  (x.ip)ij= <xi,xj>_{HT}
# x : list of x
# KT : gram matrix of T  which is used in the inner product
###############################################################
x.ip2 <- function(x, tt, xdata=NULL, t.unique, KT, et=0.02, tequal=FALSE, dim=2) { 
  n <- length(tt)
  xxip <- matrix(0, n, n)
  if(tequal){ 
    pp <- length(t.unique)
    V.inv <- matpower2(KT,-1,et)
    if(dim==1) {
      xmat <- makematrix(x, pp)
      xxip <- xmat %*% V.inv %*% KT %*% V.inv %*% t(xmat)
    }
    if(dim>1) {
      for(d in 1:dim) {
        xmat <- makematrix(xdata[,(2+d)],pp)
        xxip <- xxip + xmat %*% V.inv %*% KT %*% V.inv %*% t(xmat)
      }
    }
  }
  if(!tequal){
    for(i in 1:(n-1)){ 
      ni <- length(tt[[i]])
      index.i <- rep(0,ni)
      for(k in 1:ni){
        index.i[k] <- which(tt[[i]][k]==t.unique)
      }
      Vii <- KT[index.i,index.i]
      Vii.inv <- matpower2(Vii,-1,et)
      if(dim==1) xxip[i,i] <- t(x[[i]]) %*% Vii.inv %*% Vii %*% Vii.inv %*% x[[i]]
      if(dim>1) {
        for(d in 1:dim)  xxip[i,i] <-  xxip[i,i] + t(x[[i]][,d]) %*% Vii.inv %*% Vii %*% Vii.inv %*% x[[i]][,d]
      }
      
      for(j in (i+1):n){
        nj <- length(tt[[j]])
        index.j <- rep(0,nj)
        for(k in 1:nj){
          index.j[k] <- which(tt[[j]][k]==t.unique)
        }
        Vjj <- KT[index.j,index.j]
        Vij <- KT[index.i,index.j]
        Vjj.inv <- matpower2(Vjj,-1,et)
        if(dim==1) xxip[i,j] <- t(x[[i]]) %*% Vii.inv %*% Vij %*% Vjj.inv %*% x[[j]]
        if(dim>1){
          for(d in 1:dim) xxip[i,j] <- xxip[i,j] + t(x[[i]][,d]) %*%  Vii.inv %*% Vij %*% Vjj.inv %*% x[[j]][,d]
        }
      }
    }
    
    ## x.ip [n , n]
    ni <- length(tt[[n]])
    index.n <- rep(0,ni)
    for(k in 1:ni){
      index.n[k] <- which(tt[[n]][k]==t.unique)
    }
    Vnn <- KT[index.n,index.n]
    Vnn.inv <- matpower2(Vnn,-1,et)
    for(d in 1:dim) xxip[n,n] <- xxip[n,n] + t(x[[n]][,d]) %*% Vnn.inv %*% Vnn %*% Vnn.inv %*% x[[n]][,d]
    
    xxip <- xxip + t(xxip) - diag(diag(xxip))  
    
  }
  
  return(xxip)
}


##################################################################
# sym(a) : make a symmertric------
##################################################################
sym <- function(a){
  return((a + t(a))/2)}

##################################################################
# spearman(x,y) : spearman correlation------
##################################################################
spearman = function(x,y){
  return(cor(rank(x),rank(y)))
}
##################################################################
# normmat(x, ftn=FALSE, xip) ------
##################################################################
# Input - x is matrix (n x p)
#       - x is a function
#            - xip : n x n matrix xip_ij = <x_i, x_j>
# Output - n x n matrix  out_ij = \| x_i - x_j \|
##################################################################
normmat <- function(x, ftn=FALSE, xip) {
  if(!ftn){
    xx <- x %*% t(x)
    dxx <- diag(xx)
    normmat.x <- t(-2 * xx + dxx) + dxx  
    normmat.x[abs(normmat.x)<(10^(-15))]<-0
  }
  if(ftn){
    normmat.x <- t(-2 * xip + diag(xip)) + diag(xip)
  }
  
  return(normmat.x)
}

##################################################################
# makelist(x) : convert a matrix x  into list------
##################################################################
# Input - x is matrix (n x p)
# Output - list of length n, for each i,  out[[i]] = x[i,]
##################################################################
makelist <- function(x){
  lapply(seq_len(nrow(x)), function(i) x[i,])  
}

##################################################################
# makematrix(x) : convert a dataframe/list x  into matrix------
##################################################################
# Input - x is a list (n x  p)  // dataframe with $x as vector
# Output - matrix of n x p
##################################################################
makematrix <- function(x, p){
  if(is.data.frame(x)) matrix(x$x, ncol = p, byrow = TRUE)
  else if(is.list(x)) matrix(unlist(x), ncol = p, byrow = TRUE)
  else if(is.vector(x)) matrix(x, ncol = p, byrow = TRUE)
}

##################################################################
# make.dataframe  : matrix -> x
##################################################################
##################################################################
makedataframe <- function(x,tt){
  n <- nrow(x);nt <- ncol(x)
  xx <- NULL
  xx$id <- rep(1:n,each=nt)
  xx$tt <- as.vector(t(tt))
  xx$x <- as.vector(t(x))
  xx <- data.frame(xx)
  return(xx)
}
##################################################################
# optionout
##################################################################
##################################################################
optionout <- function(x){
  cat(gsub(',',';',x))
}
##################################################################
# choose.tset(A,B)  # choose #(A) of B that is closest to A
##################################################################
# Input - A : vector, B : vector  (or a set )
# Output - #(A) of B that is closest to A
##################################################################
choose.tset <- function(A,B){
  ni <- length(A)
  index <- rep(0,ni)
  for(i in 1:ni){
    index[i] <- which.min(abs(A[i] - B))
  }
  return(list(index=index, val=B[index]))
}

##################################################################
# Simulate Matern process
##################################################################
# matern.ftn(J, N, nu, ssq=1)
# Input : J, N, nu, ssq ; 
#  - J : number of time points
#  - N : number of curves
#  - nu : smoothness parameter
#  - ssq : variance parameter
###################################################################
# Output : J x N matrix X[,i] : observed i-th curve on 1:J/J \in [0,1]
###################################################################
cov.ft <- function(d, nu, ssq=1){ # when nu > 0 
  td <- sqrt(2*nu)*d
  ssq * 1/(gamma(nu)*2^(nu-1)) * (td)^nu * besselK(td,nu)
}

mppower <- function(matrix,power,ignore){ # symmetric matrix
  eig <- eigen(matrix, symm=TRUE)
  eval <- eig$values
  evec <- eig$vectors
  m <- length(eval[eval>ignore])
  tmp <- evec[,1:m]%*%diag(eval[1:m]^power)%*% t(evec[,1:m])
  return(tmp)
}

matern.ftn <- function(J,N, nu, ssq=1){   # J : time, N : number of curves
  cov.mat <- matrix(0, J, J)
  tt <- 1:J/J
  cov.vec <- cov.ft(d=tt[2:J], nu=nu, ssq=ssq) 
  for(j in 1:J-1){
    cov.mat[j,((j+1):J)] <- cov.vec[1:(J-j)]
  }
  cov.mat <- cov.mat + t(cov.mat) + diag(1, J)
  cov.mat.sq <- mppower(cov.mat,1/2,10^(-14))
  #cov.array.sq <- mppower(cov.array, 1/2, 10^(-16))
  return((cov.mat.sq %*% matrix(rnorm(N*J), J, N)))
  #   plot.ts(Matern.array, 
  #           plot.type="single", 
  #           ylab = "Matern process",
  #           col = 1:N,  
  #           main=paste("nu = ",nu[i]))
}
#plot.ts(matern.ftn(J=30,N=100,nu=1, ssq=1), plot.type="single",col=1:100)


##################################################################
# bdiag2 : Construct Block Diagonal
##################################################################
# Input - M : matrix, d : the number of matrices 
# Input - Or M N are matrix 
# Output - diag(M,M,M,..., M) or diag(M,N) 
##################################################################
bdiag2 <- function(M, d, N=NULL){
  nc <- ncol(M); nr <- nrow(M)
  if(is.null(N)){
    out <- matrix(0, nr*d, nc*d)
    for(i in 1:d){
      out[(((i-1)*nr + 1) : (i*nr)), (((i-1)*nc + 1) : (i*nc))] <- M
    }
  }
  if(!is.null(N)){
    if(!is.null(M)){
      nc2 <- ncol(N); nr2 <- nrow(N)
      out <- matrix(0, nr+nr2, nc+nc2)
      out[(1:nr),(1:nc)] <- M
      out[((nr+1):(nr+nr2)), ((nc+1):(nc+nc2))] <- N  
    }
    if(is.null(M)) out <- N
  }
  return(out)
}

########
# Random function
#######
random.beta <- function(N, gamma, sigmasq,  fdrobj){
  tt <- fdrobj$tt
  KT <- fdrobj$KT
  KT.inv <- fdrobj$KT.inv
  Sxx.inv <- fdrobj$Sxx.inv
  a <- rnorm(N,0,sigmasq)
  ttt <- runif(N)
  x <- rep(0,length(tt))
  for(k in 1:N){
    x <- x + a[k] * exp(-gamma*(ttt[k]-tt)^2)
  }
  return(list(fhat=Sxx.inv %*% KT.inv %*% x, tt=tt, KT=KT, KT.inv=KT.inv, et=fdrobj$et) )
}


#-----------------------------------------------------------#
# diag.mat(A, q)  :  diag(A, q)
#-----------------------------------------------------------#
diag.mat = function(A, q){
  n=nrow(A);p=ncol(A)
  out = matrix(0, n*q, p*q)
  for(i in 1:q){
    out[ (n*(i-1)+1) : (n*i), (p*(i-1)+1):(p*i)] = A
  }
  return(out)
}


#-----------------------------------------------------------#
# plotpred----
#-----------------------------------------------------------#


plotpred <- function(x, n.tuning, vec, xlab="x", ylab="y", main="", ind.info = NULL){
  n <- nrow(x)
  if(!is.null(ind.info)){
    ind.info <- as.numeric(factor(ind.info))  # change indicator number
    nclass <- length(unique(ind.info))
    plot(x[which(ind.info==1),vec[1]], x[which(ind.info==1),vec[2]],
         xlim=c(min(x[, vec[1]]),max(x[, vec[1]])), 
         ylim=c(min(x[, vec[2]]),max(x[, vec[2]])),
         xlab=xlab, ylab=ylab, main=main, pch=16)
    for(i in 2:(nclass)){
      points(x[which(ind.info==i),vec[1]], x[which(ind.info==i),vec[2]], col=i,pch=16)
    }
  }else {
    if(n %% n.tuning != 0 ) stop('n is not multiple of n.tuning')
    nclass <- n/n.tuning
    plot(x[1:n.tuning, vec[1]],x[1:n.tuning, vec[2]],
         xlim=c(min(x[, vec[1]]),max(x[, vec[1]])), 
         ylim=c(min(x[, vec[2]]),max(x[, vec[2]])),
         xlab=xlab, ylab=ylab, main=main,pch=paste(1))
    for(i in 1:(nclass-1)){
      points(x[(i*n.tuning+1):((i+1)*n.tuning), vec[1]],x[(i*n.tuning+1):((i+1)*n.tuning), vec[2]], col=i+1, pch=paste(i+1))
    }
  }
}



#-----------------------------------------------------------#
# plot function with labels----
#-----------------------------------------------------------#
# p = 1 or 2
fdplot = function(x, y,tt, n.point=NULL){
  n=dim(x)[1]; nt=dim(x)[2];p = dim(x)[3]
  if(is.null(n.point))n.point = n
  if(!is.null(dim(x))){ # is.null is required if data are a list-class
    p = dim(x)[3]
    if(is.na(p)) p = 1  
  }
  
  if(is.list(x)){
    n = lenght(x)
    i=which(y==1)[1]
    plot(x[[1]][[i]], x[[2]][[i]], xlim=c(min(x[[1]]),max(x[[1]])),ylim=c(min(x[[2]]),max(x[[2]])),pch=16,
         col = y[i])
    for(i in 1:n.point){
      points(x[[1]][[i]], x[[2]][[i]], pch=16, col=y[i])
    }
  }else if(p==2){
    i=which(y==1)[1]
    plot(x[i,,1], x[i,,2], col=y[1], xlim=c(min(x[,,1]),max(x[,,1])),ylim=c(min(x[,,2]),max(x[,,2])), pch=16)
    for(i in 1:n.point){
      points(x[i,,1], x[i,,2], col=y[i],  pch=16)
    }
    for(i in 1:n.point) lines(x[i,,1],x[i,,2], col=y[i])
  }else if(p==1){
    plot(tt,x[1,], type='l', ylim=c(min(x),max(x)), xlab="t", ylab="X(t)", main = "Model I", col=y)
    for(i in 1:n.point)lines(tt,x[i,], col=(y[i]))
  }
  
}