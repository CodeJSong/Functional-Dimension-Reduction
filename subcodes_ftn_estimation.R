#---------------------------------------------#
# functional estimation basis expansion (subcodes) by Jun Song
#---------------------------------------------#
require(fda)

#---------------------------------------------#
## 1. Function estimation ####
#---------------------------------------------#
# Function estimation-1) using fda (basis expansion)
# Function estimation-2) using rkhs (with rkhs inner product)
## 1-1. (fda) Function estimation----
# input: matrix/array -> get.fd
#        list -> get.fd.sparse
#        a single function -> get.ftn
# output: coefficient, basis
#-----------------------------------------------------------#
## 1-1-1. function: get.fd  - get.fd from matrix  / array ----
#-----------------------------------------------------------#
# input:  x: nt X n   matrix  or nt X n X p array
#         tt:  a vector of time points
# output:   coef: nbasis X ncurve  / nbasis X ncurve X n_covariates
#-----------------------------------------------------------#
get.fd = function(xraw,tt,basisname=NULL,nbasis=15,ncv=10,basis=NULL){
  nt=dim(xraw)[1];n=dim(xraw)[2];p=dim(xraw)[3]
  if(is.na(p)) p = 1
  tmin=min(tt);tmax=max(tt)
  if(is.null(basis)){
    if(is.null(basisname)) stop('error: either basis or basisname required.')
    if(basisname=='bspline'){
      basis=create.bspline.basis(c(tmin,tmax),nbasis=nbasis)
      kt = bsplinepen(basis,0)
    }else if(basisname=='fourier'){
      if(nbasis %% 2 ==0) nbasis=nbasis +1   # odd number of basis functions
      basis=create.fourier.basis(c(tmin,tmax), nbasis=nbasis)
      kt=diag(nbasis)
    }
  }
  nbasis=basis$nbasis
  
  if(p==1){
    lambda_all<- exp(seq(log(5*10^(-15)), log(5*10^(0)), len=ncv))
    nlam<-length(lambda_all)
    gcv_all<-matrix(0,n, nlam)
    ms_all_m<-array(0, c(nbasis, n,nlam))
    for(i in 1:nlam){
      lambda=lambda_all[i]
      fd_par<-fdPar(fdobj=basis,Lfdobj=2,lambda=lambda)
      tmp<-smooth.basis(tt,xraw,fdParobj=fd_par) 
      gcv_all[,i]<-tmp$gcv
      ms_all_m[,,i]<-tmp$fd$coef
    }
    lam.ind = apply(gcv_all,1,which.min)
    lambdas<-lambda_all[lam.ind]
    xcoo = matrix(0,nbasis,n)
    for(i in 1:n){
      xcoo[,i]= ms_all_m[,i,lam.ind[i]]
    }
  }else if(p>1){
    lambda_all<- exp(seq(log(5*10^(-15)), log(5*10^(0)), len=ncv))
    nlam<-length(lambda_all)
    gcv_all<-array(0,c(n,p, nlam))
    ms_all_m<-array(0, c(nbasis, n,p,nlam))
    for(i in 1:nlam){
      lambda=lambda_all[i]
      fd_par<-fdPar(fdobj=basis,Lfdobj=2,lambda=lambda)
      tmp<-smooth.basis(tt,xraw,fdParobj=fd_par) 
      gcv_all[,,i]<-tmp$gcv
      ms_all_m[,,,i]<-tmp$fd$coef
    }
    lam.ind = apply(gcv_all,c(1,2),which.min)
    
    lambdas<-matrix(0,n,p)
    xcoo = array(0,c(nbasis,n,p))
    for(i in 1:n){
      lambdas[i,] = lambda_all[lam.ind[i,]]
      xcoo[,i,]= ms_all_m[,i,,lam.ind[i]]
    }
  }
  
  return(list(coef=xcoo,lambdas=lambdas, basis=basis))
}

#-----------------------------------------------------------#
##  1-1-2. function: get.fd.sparse -- x-t as list. ----
#-----------------------------------------------------------#
# x^1 = x[[1]] 
# x^2 = x[[2]]
# ...
# x[[j]][[i]] = x^j_i(t)   i=1,...,n,, j=1,...,p
# x[[1]][[1]] = x^1_1(t)
# x[[1]][[2]] = x^1_2(t) ...   
#-----------------------------------------------------------#
get.fd.sparse = function(x, tt,basisname=NULL,nbasis=15,ncv=10){
  p = length(x)
  n = length(x[[1]])
  
  tmin=min(unlist(tt));tmax=max(unlist(tt))
  
  if(basisname=='bspline'){
    basis=create.bspline.basis(c(tmin,tmax),nbasis=nbasis)
    kt = bsplinepen(basis,0)
  }else if(basisname=='fourier'){
    if(nbasis %% 2 ==0) nbasis=nbasis +1   # odd number of basis functions
    basis=create.fourier.basis(c(tmin,tmax), nbasis=nbasis)
    kt=diag(nbasis)
  }
  
  xcoo = array(0,c(nbasis, n, p))
  for(j in 1:p){
    coef_j = foreach(i=1:n,  .combine='cbind', .export = c("get.fd.one"), .packages = c("fda")) %dopar%
      get.fd.one(x[[j]][[i]], tt[[j]][[i]], basis, ncv)
    xcoo[,,j]= coef_j
  }
  return(list(coef=xcoo, basis=basis))
}
## 1-1-3. other related functions----
#-----------------------------------------------------------#
#       function: get.ftn   // for one function
#-----------------------------------------------------------#
get.ftn = function(x,tt,basis=NULL,basisname=NULL,nbasis=NULL,ncv=10){
  require(fda)
  nt = length(tt)
  # Get basis functoins and gram matrix
  if(is.null(basis)){
    tmin=min(tt);tmax=max(tt)
    
    if(basisname=='bspline'){
      basis=create.bspline.basis(c(tmin,tmax),nbasis=nbasis)
      kt = bsplinepen(basis,0)
    }else if(basisname=='fourier'){
      if(nbasis %% 2 ==0) nbasis=nbasis +1   # odd number of basis functions
      basis=create.fourier.basis(c(tmin,tmax), nbasis=nbasis)
      kt=diag(nbasis)
    }
  }
  nbasis=basis$nbasis
  lambda_all<- exp(seq(log(5*10^(-15)), log(5*10^(0)), len=ncv))
  nlam<-length(lambda_all)
  gcv_all<-rep(0, nlam)
  ms_all_m<-matrix(0, nbasis,nlam)
  for(i in 1:nlam){
    lambda=lambda_all[i]
    fd_par<-fdPar(fdobj=basis,Lfdobj=2,lambda=lambda)
    tmp<-smooth.basis(tt,x,fdParobj=fd_par) 
    gcv_all[i]<-tmp$gcv
    ms_all_m[,i]<-tmp$fd$coef
  }
  
  lam.ind = which.min(gcv_all)
  lambdas<-lambda_all[lam.ind]
  xcoo = ms_all_m[,lam.ind]
  
  return(list(coef=xcoo,lambdas=lambdas, basis=basis))
}




#-----------------------------------------------------------#
#       function: get.fd.missing  - get.fd with a few NA from matrix  / array
#-----------------------------------------------------------#
get.fd.missing = function(x,tt,basisname,nbasis,ncv=10){
  require(fda)
  nt=dim(x)[1];n=dim(x)[2];p=dim(x)[3]
  if(is.na(p)) p = 1
  tmin=min(tt);tmax=max(tt)
  if(basisname=='bspline'){
    basis=create.bspline.basis(c(tmin,tmax),nbasis=nbasis)
    kt = bsplinepen(basis,0)
  }else if(basisname=='fourier'){
    if(nbasis %% 2 ==0) nbasis=nbasis +1   # odd number of basis functions
    basis=create.fourier.basis(c(tmin,tmax), nbasis=nbasis)
    kt=diag(nbasis)
  }
  
  if(p==1){
    if(sum(is.na(x))>0){
      xcoo = get.fd.p1(x,tt,basis,ncv)
    }else{
      xcoo = get.fd(x,tt,basisname,nbasis,ncv)$coef
    }
  }else if(p>1){
    xcoo = array(0,c(nbasis,n,p))
    missing.ind = unique(which(is.na(x),arr.ind=T)[,3])
    xcoo[,,-missing.ind] = get.fd(x[,,-missing.ind],tt,basisname,nbasis,ncv)$coef
    nna = length(missing.ind)
    xna = array(0,c(nbasis,n,nna))
    
    for(j in 1:nna){
      xna[,,j] = get.fd.p1(x[,,missing.ind[j]],tt,basis,ncv)
    }
    xcoo[,,missing.ind] = xna
  }
  return(list(coef=xcoo,basis=basis))
}
get.fd.p1 = function(x,tt,basis,ncv){
  nt=dim(x)[1];n=dim(x)[2]
  nbasis=basis$nbasis
  na.ind = which(is.na(x),arr.ind=T)[,2]
  xok = x[,-na.ind]
  xcoo.ok=get.fd(xok,tt, ncv=ncv,basis=basis)$coef
  xna = x[,na.ind]
  nna = ncol(xna)
  xcoo.na = matrix(0,nbasis,nna)
  for(kk in 1:nna){
    xx = xna[,kk]
    xcoo.na[,kk] =get.ftn(xx[!is.na(xx)], tt[!is.na(xx)],basis=basis,ncv=ncv)$coef
  }
  xcoo = matrix(0,nbasis,n)
  xcoo[,-na.ind] = xcoo.ok
  xcoo[,na.ind] = xcoo.na
  return(xcoo)
}

#-----------------------------------------------------------#
#       function: get.fd.one -- one function
#-----------------------------------------------------------#
# input :  x : vector , tt: vector
#-----------------------------------------------------------#
# x^1 = x[[1]] 
get.fd.one = function(x,tt,basis,ncv){
  lambda_all<- exp(seq(log(5*10^(-15)), log(5*10^(0)), len=ncv))
  nlam<-length(lambda_all)
  gcv_all<-rep(0, nlam)
  nbasis=basis$nbasis
  ms_all_m<-matrix(0, nbasis, nlam)
  for(i in 1:nlam){
    lambda=lambda_all[i]
    fd_par<-fdPar(fdobj=basis,Lfdobj=2,lambda=lambda)
    tmp<-smooth.basis(tt,x,fdParobj=fd_par) 
    gcv_all[i]<-tmp$gcv
    ms_all_m[,i]<-tmp$fd$coef
  }
  lam.ind = which.min(gcv_all)
  return(ms_all_m[,lam.ind])
}


#-----------------------------------------------------------#
#  fpca: functional PCA
#-----------------------------------------------------------#
# input :  ftn: get.fd object    (requires "coef" and "basis")
#-----------------------------------------------------------#
fpca = function(ftn){
  temp = ftn$coef
  n = dim(temp)[2]; p = dim(temp)[3]
  
  if(is.na(p)) p = 1
  xcoef = temp
  
  basis =ftn$basis
  basisname=basis$type
  if(basisname=="bspline") GB = bsplinepen(basis,0)
  if(basisname=="fourier") GB = fourierpen(basis,0)
  
  
  one = matrix(1, n, 1)
  Q = diag(n) - one %*% t(one) / n
  B.half = matpower(GB,0.5)
  
  if(p==1){
    Sigma = B.half %*% xcoef %*% Q %*% t(xcoef) %*% B.half / n
    egn = eigen(Sigma, sym=TRUE)
    B.inv.half = matpower(GB,-0.5)
    pred = Q %*% t(xcoef) %*% GB %*% B.inv.half %*% egn$vec
    out = list(pred=pred, eval=egn$val, GB=GB, mat =  B.inv.half %*% egn$vec)
  }else if(p>1){ # BX is the same for now...
    M.half = B.half %*% xcoef[,,1] %*% Q
    B.inv.half = matpower(GB,-0.5)
    D.inv.half = B.inv.half
    BB = GB %*% xcoef[,,1] %*% Q
    for(j in 2:p){
      D.inv.half = as.matrix(bdiag(D.inv.half, B.inv.half))
      M.half = rbind(M.half, B.half %*% xcoef[,,j] %*% Q)
      BB = rbind(BB, GB %*% xcoef[,,j] %*% Q)
    }
    DD = M.half
    M.half = M.half %*% t(M.half)/n
    egn = eigen(M.half, sym=TRUE)
    pred = t(BB) %*% D.inv.half %*% egn$vec
    out = list(pred=pred, eval=egn$val,GB=GB,  mat=D.inv.half %*% egn$vec)
  }
  return(out)
}



#-----------------------------------------------------------#
# get.fd.fpc.uni   convert the basis to fpc -- univariate functional data
#-----------------------------------------------------------#
# input :  ftn: get.fd object    (requires "coef" and "basis")
# output : coef: M X n matrix coefficient
#          fpc.coef: fpca result
#          fpcmat: conversion matrix
#-----------------------------------------------------------#
get.fd.fpc.uni = function(ftn){
  # basic
  temp = ftn$coef
  n = dim(temp)[2]; p = dim(temp)[3]
  if(p>1) stop('use a function for multivariate functional data')
  # fpc info
  fpc.out = fpca(ftn)
  GB = fpc.out$GB
  mat = fpc.out$mat
  
  # converter
  fpcmat = t(mat) %*% GB
  
  return(list(coef = fpcmat%*% temp, fpc.coef = fpc.out$mat, fpcmat = fpcmat))
}

#-----------------------------------------------------------#
# eval.fpc  evaluate the function based on fpc-based coordinate representation
#-----------------------------------------------------------#
# input :  arg: evaluate on where  t1,t2,,..,t10
#          coef.fpc --   M X 1 vector/matrix, or M X n matrix
#          fd.fpc: get.fd.fpc object,
#          ftn: get.fd object    (requires "coef" and "basis")
# output : x(t1),..(xtn)
#-----------------------------------------------------------#

eval.fpc = function(arg, coef.fpc, fd.fpc, ftn){
  phimat = fd.fpc$fpc.coef
  out = t(coef.fpc) %*% t(phimat) %*% t(eval.basis(arg, basis))
  return(out)
}