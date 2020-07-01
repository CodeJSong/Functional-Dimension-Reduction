#############################################
# functional dimension reduction (subcodes)
#############################################
require(fda)
#####################################
# Nonlinear regression
#####################################
################################
# kernmat:  k(X,Y)  = exp(-\|X-Y\|^2/h)
################################
# xtype=scalar -> multivariate x
# is.fd(x) -> fda package object
# xtype='function' -> matrix x, but functional object (different inner product)
kernmat = function(x, y=NULL, h=.1, xtype='scalar', kernt=NULL, gamma=1){
  normxy = normxy(x, y, h, xtype, kernt)
  return(exp(-gamma*normxy/h))
}
normxy = function(x, y=NULL, h=.1, xtype='scalar', kernt=NULL){
  if(xtype=='function' & (!is.fd(x))){
    n1=dim(x)[1]
    nt=dim(x)[2]
    if(is.null(y)) y=x
    n2=dim(y)[1]
    tto=seq(from=0,to=1,length=nt)
    kt=gramt(tto,tto,kernt)
    xx = x%*% kt %*% t(x)
    yy = y%*% kt %*% t(y)
    xy = x %*% kt %*% t(y)
  }
  if(xtype=='scalar'){
    n1=dim(x)[1]
    if(is.null(y)) y=x
    n2=dim(y)[1]
    if(is.null(n1)){
      n1 = length(x); n2=length(y)
    }
    xx = x%*%t(x)
    yy = y%*%t(y)
    xy = x %*% t(y)
    
  }
  if(is.fd(x)){
    if(is.null(y)) y = x
    x.tes = y
    n=dim(x$coef)[2]
    px=dim(x$coef)[3]
    if(is.na(px)) px=1
    basis.name = x$basis$type
    basis = x$basis
    nbasis = x$basis$nbasis
    xcoo = (x$coef)
    if(basis.name == "bspline")kt = bsplinepen(basis,0)
    if(basis.name == "fourier") kt = diag(nbasis)
    x.test.coo = (x.tes$coef)
    n.test=dim(x.test.coo)[2]
    
    if(px>1){
      stop('multiple function is not coded')
    }else if(px==1) {
      xcoo = t(xcoo)
      x.test.coo = t(x.test.coo)
    }
    x = xcoo
    y = x.test.coo
    n1=n;n2=n.test
    xx = x%*% kt %*% t(x)
    yy = y%*% kt %*% t(y)
    xy = x %*% kt %*% t(y)
    
  }
  xxrow = matrix(diag(xx),n1,n2)
  yycol = matrix(diag(yy),n1,n2,byrow=T)
  
  normxy = xxrow + yycol - 2*xy
  normxy[normxy<0]=0
  return(normxy)
}
gammax = function(norm){
  n = (nrow(norm)+ncol(norm))/2
  sigma=sum(sqrt(norm))/(2*choose(n,2))
  return(1/(2*sigma^2))
}
###########################################################################
#       function: nonlinear.est  (nonlinear function-on-function regression)
###########################################################################
# y : response in training data
# x : predictor in training
# x.test : predictor in test
# output : y.hat - estimated y(training), y.test.hat = predicted y (test set)
nonlinear.est = function(y,x,x.test=NULL, ytype='function', xtype='scalar', kern=kern, h=0.1,cvh=TRUE, kfolds=5,ncv=20){
  if(ytype=='function')n = ncol(y)
  if(ytype=='scalar') n = length(y)
  if(is.null(x.test)) x.test=x
  
  # h adjustment
  norm.all = normxy(x,x.test, xtype=xtype, kernt = kern)
  gamma = gammax(norm.all)
  ############
  # k-fold cross validation for h over training set
  ############
  if(cvh){
    folds <- rep_len(1:kfolds, n)
    folds <- sample(folds, n)
    
    hs=exp(seq(log(10^(-6)),log(10),len=ncv))
    cvout = matrix(0,kfolds,ncv)
    
    for(k in 1:kfolds){
      
      fold <- which(folds == k)
      if(ytype=='function'){
        y.tr = y[,-fold]
        y.tes = y[,fold]
      }else if(ytype=='scalar'){
        y.tr = y[-fold]
        y.tes = y[fold]}
      
      # x
      if(!is.null(dim(x))){
        x.tr = x[-fold,]
        x.tes = x[fold,]    
      }
      
      
      if(is.null(dim(x))){
        x.tr = x[-fold]
        x.tes = x[fold]
      }
      
      for(icv in 1:ncv){
        hcv = hs[icv]
        kx = kernmat(x.tr,x.tes, h=hcv, xtype=xtype, kernt = kern, gamma)
        y.test.hat = nonlinear.comp(y.tr,kx, ytype=ytype)
        cvout[k,icv]= sum((y.test.hat - y.tes)^2)
      }
    }
    cv = apply(cvout,2,mean)
    h = hs[which.min(cv)]
  }else{
    cvout=0
  }
  kx = kernmat(x,x, h=h, xtype=xtype, kernt = kern, gamma)
  kx.test = kernmat(x,x.test, h=h, xtype=xtype, kernt = kern, gamma)
  ###
  #y hat 
  ### 
  y.hat = nonlinear.comp(y,kx,ytype=ytype)
  ### 
  # y test hat
  ### 
  y.test.hat = nonlinear.comp(y,kx.test,ytype=ytype)
  return(list(y.hat=y.hat,y.test.hat=y.test.hat, cv=cvout, h=h))
}
nonlinear.comp = function(y,kx,ytype='function'){
  kx.sum = apply(kx,2,sum)  # kx.sum[j] = sum(kx[,j])
  if(ytype=='function'){
    nt = nrow(y); n.test=ncol(kx)
    y.hat = matrix(0,nt,n.test)
    for(j in 1:n.test){
      nu = apply(t(y) * kx[,j],2,sum)
      y.hat[,j] = nu/kx.sum[j]
    }  
  }else if(ytype=='scalar'){
    n.test=ncol(kx)
    y.hat = rep(0,n.test)
    for(j in 1:n.test){
      nu = sum(y * kx[,j])
      y.hat[j] = nu/kx.sum[j]
    }
  }
  
  return(y.hat)
}

###########################################################################
#      gram matrix in the t domain, using rkhs
###########################################################################
gramt=function(tte,tto,kern){
  complexity=1 
  ltte=length(tte);ltto=length(tto)
  if (kern=="gauss"){
    a1=matrix(tte^2,ltte,ltto);a2=tte%*%t(tto);a3=t(matrix(tto^2,ltto,ltte))
    a=a1-2*a2+a3
    b1=matrix(tto^2,ltto,ltto);b2=tto%*%t(tto);b3=t(matrix(tto^2,ltto,ltto))
    b=b1-2*b2+b3
    sigma=sum(sqrt(b))/(2*choose(ltto,2));gamma=complexity/(2*sigma^2)
    ktmat=exp(-gamma*(a))}
  if(kern=="brown"){
    arr=array(0,c(ltte,ltto,2))
    arr[,,1]=matrix(tte,ltte,ltto);arr[,,2]=t(matrix(tto,ltto,ltte))
    ktmat=apply(arr,c(1,2),min)}
  return(t(ktmat))
}

gramy=function(y){
  n=length(y)
  k2=y%*%t(y);k1=t(matrix(diag(k2),n,n));k3=t(k1);k=k1-2*k2+k3
  sigma=sum(sqrt(k))/(2*choose(n,2));gamma=1/(2*sigma^2)
  return(exp(-gamma*(k1-2*k2+k3)))
}

###########################################################################
#      gram matrix in the x domain, using rkhs gauss kernel 
###########################################################################
gramx=function(x,kt, p=1,gamma=NULL){
  n=dim(x)[1]
  p=dim(x)[3]
  if(is.na(p)) p=1
  if(p==1){
    xip = x %*% kt %*% t(x)
  }else{
    xip = matrix(0, n, n)
    for(j in 1:p){
      xip=xip + x[,,j]%*%kt%*%t(x[,,j])
    }  
  }
  k = t(-2 * xip + diag(xip)) + diag(xip)
  k[k<0]=0
  if(is.null(gamma)){
    sigma=sum(sqrt(k))/(2*choose(n,2));
    gamma=1/(2*sigma^2)
  }
  return(list(gram=exp(-gamma*(k)),gam=gamma, norm=k))
}
gam = function(k){
  k[k<0]=0
  sigma=sum(sqrt(k))/(2*choose(n,2));
  gamma=1/(2*sigma^2)
  return(gamma)
}


###########################################################################
#       function: discrete kernel
###########################################################################
gram.dis=function(y){
  n=length(y);yy=matrix(y,n,n);diff=yy-t(yy);vecker=rep(0,n^2)
  vecker[c(diff)==0]=1;vecker[c(diff)!=0]=0 
  return(matrix(vecker,n,n))}


###########################################################################
#       function: Moor-Panrose power
###########################################################################
mppower=function(matrix,power,ignore=10^(-15)){
  eig = eigen(matrix)
  eval = eig$values
  evec = eig$vectors
  m = length(eval[abs(eval)>ignore])
  tmp = evec[,1:m]%*%diag(eval[1:m]^power)%*%t(evec[,1:m])
  return(tmp)
}

mppower2=function(matrix,power,ignore=10^(-14)){
  tmp = svd(matrix)
  eval = tmp$d
  evec = tmp$u
  eval[abs(eval)<ignore]=0
  tmp = evec%*%diag(eval^power)%*%t(evec)
  return(tmp)
}

###########################################################################
#       function: matpower
###########################################################################

matpower=function(a,alpha){
  a=(a+t(a))/2
  tmp=eigen(a,sym=T)
  return(tmp$vectors%*%diag((tmp$values)^alpha)%*%t(tmp$vectors))}
matpower2 = function(a, alpha, eps)
  matpower(a+eps*onorm(a)*diag(nrow(a)),alpha)
###########################################################################
#       function: operator norm of a matrix
###########################################################################
onorm=function(a) return(eigen(sym(a), sym=T)$values[1])
###########################################################################
#       function: symmetrize a matrix 
###########################################################################
sym=function(a) return(round((a+t(a))/2,8))



###########################################################################
#       function: get.fd  - get.fd from matrix  / array
###########################################################################
get.fd = function(x,tt,kern=NULL,nbasis=NULL,ncv=10,basis=NULL){
  require(fda)
  nt=dim(x)[1];n=dim(x)[2];p=dim(x)[3]
  if(is.na(p)) p = 1
  tmin=min(tt);tmax=max(tt)
  if(!is.null(kern)){
    if(kern=='bspline'){
      basis=create.bspline.basis(c(tmin,tmax),nbasis=nbasis)
      kt = bsplinepen(basis,0)
    }else if(kern=='fourier'){
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
      tmp<-smooth.basis(tt,x,fdParobj=fd_par) 
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
      tmp<-smooth.basis(tt,x,fdParobj=fd_par) 
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
  
  return(list(coef=xcoo,lambdas=lambdas, kt=kt, basis=basis))
}



###########################################################################
#       function: get.fd.missing  - get.fd with a few NA from matrix  / array
###########################################################################
get.fd.missing = function(x,tt,kern,nbasis,ncv=10){
  require(fda)
  nt=dim(x)[1];n=dim(x)[2];p=dim(x)[3]
  if(is.na(p)) p = 1
  tmin=min(tt);tmax=max(tt)
  if(kern=='bspline'){
    basis=create.bspline.basis(c(tmin,tmax),nbasis=nbasis)
    kt = bsplinepen(basis,0)
  }else if(kern=='fourier'){
    if(nbasis %% 2 ==0) nbasis=nbasis +1   # odd number of basis functions
    basis=create.fourier.basis(c(tmin,tmax), nbasis=nbasis)
    kt=diag(nbasis)
  }
  
  if(p==1){
    if(sum(is.na(x))>0){
      xcoo = get.fd.p1(x,tt,basis,ncv)
    }else{
      xcoo = get.fd(x,tt,kern,nbasis,ncv)$coef
    }
  }else if(p>1){
    xcoo = array(0,c(nbasis,n,p))
    missing.ind = unique(which(is.na(x),arr.ind=T)[,3])
    xcoo[,,-missing.ind] = get.fd(x[,,-missing.ind],tt,kern,nbasis,ncv)$coef
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
###############################################################
#       get.fd: get functional object
###############################################################
get.fd.list = function(x,tt, basis,ncv=20){
  p=length(x)
  n=length(x[[1]])
  nbasis=basis$nbasis
  coefs = array(0,c(nbasis,n,p))
  for(j in 1:p){
    for(i in 1:n){
      coefs[,i,j]= get.fd1(x[[j]][[i]], tt[[j]][[i]], basis=basis,ncv=ncv)$coef
    }
  }
  return(list(coef=coefs, basis=basis))
}
get.fd1 = function(x, tt, basis=NULL, kern='bspline', nbasis=10, ncv=20){
  # 1 function
  # x: vector, tt: vector
  require(fda)
  if(is.null(basis)){
    tmin=min(tt);tmax=max(tt)
    if(kern=='bspline'){
      basis=create.bspline.basis(c(tmin,tmax),nbasis=nbasis)
      kt = bsplinepen(basis,0)
    }else if(kern=='fourier'){
      if(nbasis %% 2 ==0) nbasis=nbasis +1   # odd number of basis functions
      basis=create.fourier.basis(c(tmin,tmax), nbasis=nbasis)
      kt=diag(nbasis)
    }  
  }else{
    nbasis=basis$nbasis
    if(basis$type=='bspline') kt = bsplinepen(basis,0)
    if(basis$type=='fourier')kt=diag(nbasis)
  }
  lambda_all<- exp(seq(log(5*10^(-15)), log(5*10^(0)), len=ncv))
  nlam<-length(lambda_all)
  gcv_all<-matrix(0,nlam)
  ms_all_m<-array(0, c(nbasis, nlam))
  for(i in 1:nlam){
    lambda=lambda_all[i]
    fd_par<-fdPar(fdobj=basis,Lfdobj=2,lambda=lambda)
    tmp<-smooth.basis(tt,x,fdParobj=fd_par) 
    gcv_all[i]<-tmp$gcv
    ms_all_m[,i]<-tmp$fd$coef
  }
  lam.ind = which.min(gcv_all)
  lambda<-lambda_all[lam.ind]
  xcoo = ms_all_m[,lam.ind]
  
  
  return(list(coef=xcoo,lambdas=lambda, kt=kt, basis=basis))
}
###########################################################################
#       function: get.ftn   // for one function
###########################################################################

get.ftn = function(x,tt,basis=NULL,kern=NULL,nbasis=NULL,ncv=10){
  require(fda)
  nt = length(tt)
  # Get basis functoins and gram matrix
  if(is.null(basis)){
    tmin=min(tt);tmax=max(tt)
    
    if(kern=='bspline'){
      basis=create.bspline.basis(c(tmin,tmax),nbasis=nbasis)
      kt = bsplinepen(basis,0)
    }else if(kern=='fourier'){
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



###############################################################
# function: Q = I - J/n
###############################################################
qmat = function(n) return(diag(n)-rep(1,n)%*%t(rep(1,n))/n)



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

###########################################################################
#      compute the coordinate of x or y for one function
###########################################################################
xcoovec=function(f,tte,tto,ridge,kern){
  kt=gramt(tto,tto,kern)
  scale=eigen(kt)$values[1]
  ktinv=matpower(kt+scale*ridge*diag(nrow(kt)),-1)
  kt1=t(gramt(tte,tto,kern))
  coo=c(ktinv%*%f);evalu=c(kt1%*%ktinv%*%f)
  return(list(coo=coo,evalu=evalu))}
###########################################################################
#       compute the coordinate of x or y for n functions
###########################################################################
xcoomat=function(ff,tte,ttt,ridge,kern){
  n=dim(ff)[1];pf=dim(ff)[3];nb=dim(ttt)[2]
  if(is.na(pf)) pf=1
  if(pf==1){
    ffcoo=numeric()
    for(i in 1:n) ffcoo=rbind(ffcoo,xcoovec(ff[i,],tte,ttt[i,],ridge,kern)$coo)  
  }else if(pf>1){
    ffcoo = array(0,c(n,nb,pf))
    for(d in 1:pf){
      for(i in 1:n){
        ffcoo[i,,d]=xcoovec(ff[i,,d],tte,ttt[i,],ridge,kern)$coo
      }
    }
  }
  
  return(ffcoo)
}

###############################################################
# diag.mat(A, q)  :  diag(A, q)
##############################################################
diag.mat = function(A, q){
  n=nrow(A);p=ncol(A)
  out = matrix(0, n*q, p*q)
  for(i in 1:q){
    out[ (n*(i-1)+1) : (n*i), (p*(i-1)+1):(p*i)] = A
  }
  return(out)
}
###########################################################################
#      gcv tuning for ex and ey when y is scalar or function
###########################################################################
gcvxy=function(x,y,eps,which,kernt,et,ytype){
  nt=dim(x)[2];n=dim(x)[1];Q=diag(n)-rep(1,n)%*%t(rep(1,n))/n
  tto=seq(from=0,to=1,length=nt);tte=seq(from=0,to=1,
                                         length=2*(nt-1)+1);ttt=t(matrix(tto,nt,n))
  kt=gramt(tto,tto,kernt)
  if(which=="ey"&ytype=="scalar") 
  {xcoo=xcoomat(x,tte,ttt,et,kernt);G1=gramx(xcoo,kt)$gram;G2=gramy(y)}
  if(which=="ey"&ytype=="function") 
  {xcoo=xcoomat(x,tte,ttt,et,kernt);G1=gramx(xcoo,kt)$gram
  ycoo=xcoomat(y,tte,ttt,et,kernt);G2=gramx(ycoo,kt)$gram}
  if(which=="ex"&ytype=="scalar")  
  {xcoo=xcoomat(x,tte,ttt,et,kernt);G2=gramx(xcoo,kt)$gram;G1=gramy(y)}
  if(which=="ex"&ytype=="function") 
  {xcoo=xcoomat(x,tte,ttt,et,kernt);G2=gramx(xcoo,kt)$gram
  ycoo=xcoomat(y,tte,ttt,et,kernt);G1=gramx(ycoo,kt)$gram}
  if(which=="ex"&ytype=="categorical"){
    xcoo=xcoomat(x,tte,ttt,et,kernt)
    G2=gramx(xcoo,kt)$gram
    G1=gram.dis(y)
  }
  if(which=="ey"&ytype=="categorical"){
    xcoo=xcoomat(x,tte,ttt,et,kernt)
    G1=gramx(xcoo,kt)$gram
    G2=gram.dis(y)
  }
  G2inv=matpower(G2+eps*onorm(G2)*diag(n),-1)
  nu=sum((G1-G2%*%G2inv%*%G1)^2) 
  tr=function(a) return(sum(diag(a)))
  de=(1-tr(G2inv%*%G2)/n)^2   
  out = nu/de
  return(out)
}
###########################################################################
#      gcv tuning for et
###########################################################################
gcvt=function(x,ttt,et,kern){  
  n=dim(x)[1];nt=dim(x)[2] 
  px=dim(x)[3]
  
  kt=gramt(ttt[1,],ttt[1,],kern)# assume they are balanced
  ktinv=matpower(kt+et*onorm(kt)*diag(nt),-1)
  if(is.na(px)) px=1
    
  nuset=numeric();deset=numeric()
  if(px==1){
    for(i in 1:n){
      nuset=c(nuset,sum((x[i,]-kt%*%ktinv%*%x[i,])^2))
      tr=function(a) return(sum(diag(a)))
      if(kern=="brown") deset=c(deset,(1-tr(kt%*%ktinv)/(nt-1))^2)
      if(kern=="gauss") deset=c(deset,(1-tr(kt%*%ktinv)/(nt))^2)
    }   
  }else if(px>1){
    for(j in 1:px){
      for(i in 1:n){
        nuset=c(nuset,sum((x[i,,j]-kt%*%ktinv%*%x[i,,j])^2))
        tr=function(a) return(sum(diag(a)))
        if(kern=="brown") deset=c(deset,(1-tr(kt%*%ktinv)/(nt-1))^2)
        if(kern=="gauss") deset=c(deset,(1-tr(kt%*%ktinv)/(nt))^2)
      }   
    }
  }
  
  return(sum(nuset/deset))
}
###########################################################################
#      gcv tuning for ex when y is categorical 
###########################################################################
gcvx.categorical=function(x,y,ex,et,kernt){
  nt=dim(x)[2];n=dim(x)[1];Q=diag(n)-rep(1,n)%*%t(rep(1,n))/n
  tto=seq(from=0,to=1,length=nt);tte=seq(from=0,to=1,
                                         length=2*(nt-1)+1);ttt=t(matrix(tto,nt,n))
  kt=gramt(tto,tto,kernt)
  xcoo=xcoomat(x,tte,ttt,et,kernt)
  G2=gramx(xcoo,kt)$gram
  G2inv=matpower(G2+ex*onorm(G2)*diag(n),-1)
  G1=gram.dis(y)
  nu=sum((G1-G2%*%G2inv%*%G1)^2)
  tr=function(a) return(sum(diag(a)))
  de=(1-tr(G2inv%*%G2)/n)^2
  return(nu/de)
}
###########################################################################
#      gcv tuning for ex when y is categorical 
###########################################################################
gcvy.categorical=function(x,y,ey,et,kernt){
  nt=dim(x)[2];n=dim(x)[1];Q=diag(n)-rep(1,n)%*%t(rep(1,n))/n
  tto=seq(from=0,to=1,length=nt);tte=seq(from=0,to=1,
                                         length=2*(nt-1)+1);ttt=t(matrix(tto,nt,n))
  kt=gramt(tto,tto,kernt)
  xcoo=xcoomat(x,tte,ttt,et,kernt)
  G2=gramx(xcoo,kt)$gram
  G2inv=matpower(G2+ex*onorm(G2)*diag(n),-1)
  G1=gram.dis(y)
  nu=sum((G1-G2%*%G2inv%*%G1)^2)
  tr=function(a) return(sum(diag(a)))
  de=(1-tr(G2inv%*%G2)/n)^2
  return(nu/de)
}
############################################################### 
#   find minimzer of a vector 
############################################################### 
minimizer=function(x,y) return(x[order(y)[1]])

###############################################################
#     function:  multiple correlation
###############################################################
mcorr=function(u,v){
  u=as.matrix(u);v=as.matrix(v) 
  if(dim(u)[2]==1) return(c(abs(cor(u,v))))
  if(dim(u)[2]!=1){
    suu.nhalf=matpower(var(u),-1/2);svv.inv=solve(var(v))
    return(c(abs(sum(diag(suu.nhalf%*%cov(u,v)%*%svv.inv%*%cov(v,u)%*%suu.nhalf)))))}
}

gcvxy.fd=function(n,kx,ky,eps,which,ytype){
  if(which=='ey') G1=kx;G2=ky
  if(which=='ex') G1=ky;G2=kx
  G2inv=mppower2(G2+eps*onorm(G2)*diag(n),-1)
  nu=sum((G1-G2%*%G2inv%*%G1)^2) 
  
  de=(1-tr(G2%*%G2inv)/n)^2 
  return(nu/de)
}
tr=function(a) return(sum(diag(a)))
gcvxy.fd2=function(n,kx,ky,eps,which,kt=NULL,kt.half=NULL,kt.nhalf=NULL, xcoo=NULL, Q=NULL){
  if(which=='ey'){
    g2 = ky
    kyinv = mppower2(g2 + eps*onorm(g2)*diag(n),-1)
    mat = kx - kyinv %*% ky %*% kx
    
    nu=sum((mat)^2) 
    de=(1-tr(kyinv %*% ky)/n)^2 
  } 
  if(which=='ex'){
    g2 = kx
    kxinv = mppower2(g2 + eps*onorm(g2)*diag(n),-1)
    mat = ky -kxinv %*% kx %*% ky
    
    nu=sum((mat)^2) 
    de=(1-tr( kxinv%*%kx)/n)^2 
  } 
  
  return(nu/de)
}

gcvt.fd = function(n,xcoo,kt,et){
  ktinv=mppower2(kt+et*onorm(kt)*diag(nrow(kt)),-1)
  nuset=numeric();deset=numeric()
  nuset=c(nuset,sum((xcoo[i,]-kt%*%ktinv%*%xcoo[i,])^2))
  tr=function(a) return(sum(diag(a)))
  deset=c(deset,(1-tr(kt%*%ktinv)/(nt))^2)
  return(sum(nuset/deset))
}
