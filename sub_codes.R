# update April 8, 2016
# - separate c code into sub_codes.cpp / bspline.cpp  
library(RcppArmadillo)
library(fda)
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}
if(get_os()=="osx"){
  Rcpp::sourceCpp('~/CloudStation/Codes/Functional Dimension Reduction/cpp_code/bspline.cpp')
  Rcpp::sourceCpp('~/CloudStation/Codes/Functional Dimension Reduction/cpp_code/sub_codes.cpp')
}
if(get_os()=="windows"){
  Rcpp::sourceCpp('C:/Users/Jignic/CloudStation/Codes/Functional Dimension Reduction/cpp_code/bspline.cpp')
  Rcpp::sourceCpp('C:/Users/Jignic/CloudStation/Codes/Functional Dimension Reduction/cpp_code/sub_codes.cpp')
}

##############################################################
# data2fdr   4/8/2016
##############################################################
# make data to list (variables for fdr)
# it can be applied linear FDR, nonlinear FDR, Kernal FPCA
##############################################################
# data frame  v1 : id, v2 : tt, v3 : x_1, v4: x_2, v5 : x_3 ....
#
# output
# 1. [x] 
# 2. <x_i, x_j>
# 3. normX = \|X_i - X_j \|^2
data2fdr <- function(x=NULL, tx=NULL, xftn=FALSE,  xtequal=FALSE, gcvall=FALSE, lambdascale=2,
                     ftnstruct="fda", sht=2, maxnbasis=50, basisname="bspline", pca.basis=FALSE,
                     dimx=1,ncv1=10, et=0.05){
  require(fda)
  fpca=NULL
  if(is.data.frame(x)){ # change x and tt as a list
    # data frame  v1 : id, v2 : tt, v3 : x_1, v4: x_2, v5 : x_3 ....
    xdata <- x
    x <- lapply(split(xdata,xdata$id), function(xdata) matrix(as.numeric(as.matrix(xdata[,3:(2+dimx)])),ncol=dimx))
    tx <-  lapply(split(xdata,xdata$id), function(xdata) as.numeric(xdata[,2]))
  }
  if(!xftn){  # x : n x p matrix
    x <- as.matrix(x)
    #if(center==TRUE) x <- center(x)
    n <- nrow(as.matrix(x))
    normX <- normmat(x)
    basis=NULL;xcoef=NULL  ## those are not needed in vectors
    ntx=NULL; xip <- tcrossprod(x)
    basis <- NULL ; xcoef=NULL;
    binnerprod=NULL;lambdaout=NULL
  }
  if(xftn){  # x : list  tx : list
    n <- length(x)
    xt.unique <- as.vector(sort(unique(unlist(tx))))
    ntx <- length(xt.unique) 
    if(ftnstruct=="rkhs"){
      normmat.t <- normmat(as.matrix(xt.unique))
      gt <- gam(normmat.t)/sht
      KTX <- exp(-gt*normmat.t)  # gram matrix for T  
      xip <- xiplist(n, x, tx, xt.unique, KTX, et, as.numeric(xtequal), dimx)
      basis <- NULL ; xcoef=NULL;
      binnerprod=KTX;lambdaout=NULL
      stop('need more coding')
    }
    if(ftnstruct=="fda"){
      tmin <- min(xt.unique) ; tmax<-max(xt.unique)
      nbasis <- ntx
      if(ntx>maxnbasis) nbasis <- maxnbasis
      if(nbasis%%2==0) nbasis <- nbasis - 1 ## odd number of basis
      if(!pca.basis){  # use the basis 
        fpca <- NULL
        basis <- eval(as.name(paste("create.",basisname,".basis",sep='')))(c(tmin,tmax), nbasis=nbasis) # maxnbasis
        binnerprod <- eval(as.name(paste(basisname,"pen",sep='')))(basis,0)
        RR<-eval(as.name(paste(basisname,"pen",sep='')))(basis,2)
        lambda_all<- exp(seq(log(10^(-8)),log(10^(-lambdascale)),length=ncv1))
        if(dimx==1){
          #X <- eval.basis(xt.unique, basis)
          #xmat <- makematrix(x,ntx)
          #### Using same smoothing parameters for all curves
          #       gcvf.mat <- bsplinegcv_mat(xmat=X, rr=RR, ystar=t(as.matrix(xmat)), lambdagrid=lambda_all)
          #       xcoef <- gcvf.mat$mucoef
          #### Using different smoothing parameters for all curves
          temps <-bspline_unbal(x,tx,basis,RR,lambdagrid=lambda_all, gcvall=as.numeric(gcvall))
          xcoef <- temps$mucoef; lambdaout <- temps$lambdaout
          #if(center==TRUE) xcoef <- t(center(t(xcoef)))
          # gcvf$mucoef : a matrix of nbasis X number of functions   [,1] : coefficients of X_1
          xip <- t(xcoef) %*% binnerprod %*% xcoef
        }
        if(dimx>1){
          # do component-wise
          xcoef <- array(0,c(dimx,nbasis,n) )
          xip <- matrix(0,n,n)
          for(d in 1:dimx){
            xd <- lapply(x,function(x) x[,d])
            temps <-bspline_unbal(xd,tx,basis,RR,lambdagrid=lambda_all, gcvall=as.numeric(gcvall))
            xcoef[d,,] <- temps$mucoef; lambdaout <- temps$lambdaout
            #if(center==TRUE) xcoef[d,,] <- t(center(t(xcoef[d,,])))
            # gcvf$mucoef : a matrix of nbasis X number of functions   [,1] : coefficients of X_1
            xip <- xip +  t(xcoef[d,,]) %*% binnerprod %*% xcoef[d,,]
          }
        }
      }
      if(pca.basis){ # use PCA basis
        nbasis <- 50
        basis <- eval(as.name(paste("create.",basisname,".basis",sep='')))(c(tmin,tmax), nbasis=nbasis) # maxnbasis
        RR<-eval(as.name(paste(basisname,"pen",sep='')))(basis,2)
        lambda_all<- exp(seq(log(10^(-8)),log(10^(-lambdascale)),length=ncv1))
        if(dimx==1){
          # get the functions as bspline basis first
          temps <-bspline_unbal(x,tx,basis,RR,lambdagrid=lambda_all, gcvall=as.numeric(gcvall))
          xcoef <- temps$mucoef; lambdaout <- temps$lambdaout
          
          # do PCA
          x.fd <- fd(xcoef,basis)
          fpca <- pca.fd(x.fd, nharm=maxnbasis)
          binnerprod <-  inprod(fpca$harmonics, fpca$harmonics)
          RR <- as.matrix(inprod(fpca$harmonics, fpca$harmonics, Lfdobj1 = 2, Lfdobj2 = 2))
          
          temps <-bspline_pca(x,tx, basis, fpca$harmonics,RR,lambdagrid=lambda_all, gcvall=as.numeric(gcvall))
          xcoef <- temps$mucoef; lambdaout <- temps$lambdaout
          xip <- t(xcoef) %*% binnerprod %*% xcoef
        }
      }
      
    }
    
    normX = t(-2 * xip + diag(xip)) + diag(xip)
  }
  if(!exists("lambdaout")) lambdaout=NULL
  return(list(n=n,ntx=ntx,xcoef=xcoef, xip=xip,  normX=normX, basis=basis, binnerprod=binnerprod,lambdaout=lambdaout, fpca=fpca))
}



##############################################################
#      function 9: center X (n by p matrix)------        
##############################################################
center <- function(x){
  return(t(t(x)-apply(x,2,mean)))}

dataframetolist <- function(xdata, dim=2){
  idlist = unique(xdata$id)
  n=length(idlist)
  x=NULL;tt=NULL;
  for(i in 1:n){
    index = which(xdata$id==idlist[i])
    tt[[i]] = as.vector(xdata[index,2])
    x[[i]] = as.matrix((xdata[index,3:(2+dim)]))
  }
  return(list(tt=tt,x=x))
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
# makelist(x) : convert a matrix x  into list------
##################################################################
# Input - x is matrix (n x p)
# Output - list of length n, for each i,  out[[i]] = x[i,]
# Missing : put 999.999 on missing part
##################################################################
makelist <- function(x, missing=FALSE){
  if(!missing){
    return(lapply(seq_len(nrow(x)), function(i) x[i,]))
  }
  if(missing){
    return(lapply(seq_len(nrow(x)), function(i) x[i,which(x[i,]!=999.999)]))
  }
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
# matpower Moore-Penrose Inverse ---- 
##################################################################
# Input - A : matrix , alpha : power, epsilon : regularization parameter
# Output - (A+epsilon * \lambda_1 I)^{+\alpha}
#    \lambda_1 is the first eigenvalue of A
##################################################################
matpower <- function(A, alpha, epsilon=0){
  if(epsilon==0) return(mppower(A,alpha))
  if(epsilon>0) return(matpower2(A, alpha, epsilon))
  if(epsilon<0) stop('epsilon should be a nonnegative real number')
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
# bdiag : Construct Block Diagonal
##################################################################
# Input - M : matrix, d : the number of matrices 
# Input - Or M N are matrix 
# Output - diag(M,M,M,..., M) or diag(M,N) 
##################################################################
bdiag <- function(M, d, N=NULL){
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



##################################################################
# Multivarate sign and rank-----
##################################################################
# input : vector or scalar
# output 
#   if scalar  : -1, 0, 1 depending on sign 
mult.sign <- function(x){  # x is vector or scalar, or n x p - matrix 
  if(is.vector(x)){
    xnorm <- sqrt(sum(x^2))
    if(xnorm == 0) return(0)
    if(xnorm !=0) return(x/xnorm)  
  }
  if(is.matrix(x)){
    out <- matrix(0,nrow(x),ncol(x))
    xnorm <- apply(x,1,function(x) sqrt(sum(x^2)))
    out[xnorm == 0] = 0
    out[xnorm !=0 ] = (x/xnorm)[xnorm!=0]
    return(out)
  }
}
mult.rank <- function(x){  # x is n-vector or n x p - matrix 
  n <- nrow(x)
  if(is.null(n)) {
    n <- length(x)
    out <- rep(0,n)
    for(i in 1:n){
      r <- rep(0,n)
      for(j in 1:n){
        r[j] <- mult.sign(x[i]-x[j])
      }
      out[i] <- mean(r)
    }
    n <- NULL
  }
  if(!is.null(n)){
    p <- ncol(x)
    out <- matrix(0,n,p)
    for(i in 1:n){
      r <- matrix(0,n,p)
      for(j in 1:n){
        r[j,] <- mult.sign(x[i,]-x[j,])
      }
      out[i,] <- apply(r,2,mean)
    }
  }
  return(out)
}

mult.cor <- function(u, v){ # v should be matrix
  if(is.matrix(u)) nu = nrow(u)
  if(is.vector(u)) nu = length(u)
  if(is.matrix(v)) nv = nrow(v)
  if(is.vector(v)) nv = length(v)
  if(nu != nv) cat("error. observed number is different")
  n = nu
  suu = var(u); svv = var(v)
  suv = cov(u,v); 
  if(is.matrix(u)) suu.inv.half = matpower(suu,-1/2)
  if(is.vector(u)) suu.inv.half = 1/sqrt(suu)
  if(is.matrix(v)) svv.inv = matpower(svv,-1)
  if(is.vector(v)) svv.inv = 1/svv
  return(sum(diag(suu.inv.half %*% suv %*% svv.inv %*% t(suv) %*% suu.inv.half)))
}



mult.cor2 <- function(model, fdrobj){ # v should be matrix
  u <- model$true.beta; vv <- fdrobj$pred; d <- model$d
  result <-rep(0,3)
  for(i in 1:3){
    v <- vv[[i]][,1:d] 
    if(is.matrix(u)) nu = nrow(u)
    if(is.vector(u)) nu = length(u)
    if(is.matrix(v)) nv = nrow(v)
    if(is.vector(v)) nv = length(v)
    if(nu != nv) cat("error. observed number is different")
    n = nu
    suu = var(u); svv = var(v)
    suv = cov(u,v); 
    if(is.matrix(u)) suu.inv.half = matpower(suu,-1/2)
    if(is.vector(u)) suu.inv.half = 1/sqrt(suu)
    if(is.matrix(v)) svv.inv = matpower(svv,-1)
    if(is.vector(v)) svv.inv = 1/svv
    result[i] <- (sum(diag(suu.inv.half %*% suv %*% svv.inv %*% t(suv) %*% suu.inv.half)))
  }
  return(result)
}

plotpred <- function(x, n.tuning, vec, xlab="x", ylab="y", main="", ind.info = NULL){
  n <- nrow(x)
  if(!is.null(ind.info)){
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

######################################
# misrate : mis classification rate of total
######################################
misrate <- function(pred, class){
  n=length(class)
  return(sum(pred!=class)/n)
}
######################################
# misrate2 : mis-classification rate by each class
# output[1,i] : False negative  // Recall
#               # of mis-classification i-th class / # of objects in i-th class (true)
# output[2,i] : False positive  // Precision
#               # of mis-classification i-th class / # of objects in i-th class (predicted)
######################################
misrate2 <- function(pred,class){
  n <- length(class)
  class.info <- levels(class)
  result <- matrix(0, 2, nclass)
  for(i in 1:length(class.info)){
    ni = sum(class==class.info[i])
    ni2 = sum(pred==class.info[i])
    result[1,i] <- sum((pred==class) & (pred==class.info[i]))/ni
    result[2,i] <- sum((pred==class) & (pred==class.info[i]))/ni2
  }
  return(result)
}


ktftn <- function(coef,ttt,tseq,gt){
  nt=length(ttt); nseq = length(tseq); n = dim(coef)[2]
  evalmat <- matrix(0, n, nseq)
  for(i in 1:n){
    for(k in 1:nt){
      evalmat[i,] <- evalmat[i,] + exp(-gt * (tseq - ttt[k])^2 ) * coef[k,i]
    }  
  }
  return(evalmat)
}


image_points = function(image, x, y, cex = 1, pos = NULL) {
  if (length(x) != length(y)) {
    stop("length(x)!=length(y): check your data")
  }
  dim.x = dim(image)[2]  #image width
  dim.y = dim(image)[1]  #image height
  if (dim.x == dim.y) {
    # obtian the ratio of width to height or height to width
    ratio.x = ratio.y = 1
  } else if (dim.x < dim.y) {
    ratio.x = dim.x/dim.y
    ratio.y = 1
  } else {
    ratio.x = 1
    ratio.y = dim.y/dim.x
  }
  cex = cex/10  #how large the image should be, divided by 10 so that it matches more closely to plotting points
  pin = par()$pin  #pin provides the width and height of the _active graphic device_
  pin.ratio = pin/max(pin)  #take the ratio
  usr = par()$usr  #usr provides the lower.x, lower.y, upper.x, upper.y values of the plotable region
  
  # combine the active device dimensions, the image dimensions, and the
  # desired output size
  image.size.y = (usr[4] - usr[3]) * pin.ratio[1] * cex
  image.size.x = (usr[2] - usr[1]) * pin.ratio[2] * cex
  for (i in 1:length(x)) {
    # plot each point pos can be NULL (default) or 1, 2, 3, or 4, corresponding
    # to centered (defualt), bottom, left, top, right, respectively.
    if (is.null(pos)) {
      # centered at (x,y), define the bottom/top and left/right boundaries of the
      # image
      x.pos = c(x[i] - (image.size.x * ratio.x)/2, x[i] + (image.size.x * 
                                                             ratio.x)/2)
      y.pos = c(y[i] - (image.size.y * ratio.y)/2, y[i] + (image.size.y * 
                                                             ratio.y)/2)
      
      rasterImage(image, x.pos[1], y.pos[1], x.pos[2], y.pos[2])
    } else if (pos == 1) {
      x.pos = c(x[i] - (image.size.x * ratio.x)/2, x[i] + (image.size.x * 
                                                             ratio.x)/2)
      y.pos = c(y[i] - (image.size.y * ratio.y), y[i])
    } else if (pos == 2) {
      x.pos = c(x[i] - (image.size.x * ratio.x), x[i])
      y.pos = c(y[i] - (image.size.y * ratio.y)/2, y[i] + (image.size.y * 
                                                             ratio.y)/2)
    } else if (pos == 3) {
      x.pos = c(x[i] - (image.size.x * ratio.x)/2, x[i] + (image.size.x * 
                                                             ratio.x)/2)
      y.pos = c(y[i], y[i] + (image.size.y * ratio.y))
    } else if (pos == 4) {
      x.pos = c(x[i], x[i] + (image.size.x * ratio.x))
      y.pos = c(y[i] - (image.size.y * ratio.y)/2, y[i] + (image.size.y * 
                                                             ratio.y)/2)
    }
    
    rasterImage(image, x.pos[1], y.pos[1], x.pos[2], y.pos[2])  #plot image
  }
}


#########################################
# detect outlier from multivariate
#########################################
# example
# outliers <- which(detect_outlier(ppp)>10)
detect_outlier <- function(x, nsd=1){
  # x is n X p matrix
  # output : detecting the outlier point (distance based)
  normx <- normmat(x)   # n X n matrix : (i,j) element = ||X_i - X_j||^2
  alldist <- normx[lower.tri(normx)]  # a vector of all distnace
  avg.dist <- mean(alldist); sd.dist <- sd(alldist) # get the mean and sd
  return(apply((normx > (avg.dist+nsd * sd.dist)), 2, sum)) # a vector : i-th element is the number of j's such that ||X_i-X_j||>avg.dist+nsd * sd.dist
}
