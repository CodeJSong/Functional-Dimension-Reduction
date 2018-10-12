nonilnear.sdr.fd = function(x,y,tt=NULL, method='fgsir',ytype='function',
                            et=0.05, ex=0.05, ey=0.05, 
                            x.tes=NULL){
  require(fda)
  test = TRUE
  if(is.null(x.tes)) {
    x.tes = x
    test = FALSE
  }
  if(ytype=="scalar") KY=gramy(y);py=1
  if(ytype=="categorical") KY=gram.dis(y);py=1

  if(is.fd(x)){
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
  }else{
    stop('x should be a functional object (fda package)')
  }
  if(is.fd(y)){
    py=dim(y$coef)[3]
    if(is.na(py)) py=1
    basis.name = y$basis$type
    nbasis = y$basis$nbasis
    ycoo = (y$coef)
    if(py>1){
      ycoo = aperm(ycoo,c(2,1,3))
      KY=gramx(ycoo,kt,py)$gram
    }else if(py==1){
      ycoo = t(ycoo) 
      KY=gramx(ycoo,kt,py)$gram
    }
  }else if(!is.fd(y) & ytype=='function'){
    stop('y should be a functional object (fda package)')
    #ycoo=xcoomat(y,tte,ttt,et,kern)
    #KY=gramx(ycoo,kt)$gram
  }
  
  kt.half=mppower(kt,1/2,10^(-8));Q=qmat(n)
  kt.nhalf=matpower(kt+et*onorm(kt)*diag(nbasis),-1/2)

  if(px>1){
    xcoo = aperm(xcoo, c(2,1,3))
    x.test.coo = aperm(x.test.coo,c(2,1,3))
  }else if(px==1) {
    xcoo = t(xcoo)
    x.test.coo = t(x.test.coo)
  }
  temp=gramx(xcoo,kt,px)
  KX = temp$gram
  KX.test = KX
  gamma = temp$gamma
  if(test){
    # if test set is different from the trainig set
    n.tot = n + n.test
    if(px>1){
      coo.all = array(0,c(n.tot, nbasis, px))
      coo.all[1:n,,] = xcoo
      coo.all[(n+1):n.tot,,] = x.test.coo
    }else if(px==1){
      coo.all = matrix(0,n.tot, nbasis)
      coo.all[1:n,] = xcoo
      coo.all[(n+1):n.tot,] = x.test.coo
    }
    KX.test = gramx(coo.all, kt, px, gamma=gamma)$gram
    KX.test = KX.test[(n+1):n.tot,1:n]
  }
  ## Estimate central class 
  M <- matrix(0, n, n)
  one <- matrix(1, n, 1)
  Q <- diag(n) - one %*% t(one) / n    
  GX <- Q %*% KX %*% Q    # centered gram matrix
  GY <- Q %*% KY %*% Q    # centered gram matrix
  if(method=="fgsir"){
    GX.inv <- matpower2(GX,-1,ex)
    GX.inv.half <- matpower2(GX,-1/2, ex)
    GX.inv.3half <- matpower2(GX,-3/2, ex)
    GX.inv.2 <- matpower2(GX,-2, ex)
    #GY.inv <- matpower2(GY%*%t(GY)+ey*diag(n),-1)%*% GY
    M <-  GX.inv.3half %*% GX %*% GY  %*% GX %*% GX.inv.3half 
    A <- eigen(sym(M), symmetric=TRUE)
    eigval <- A$values
    A <- A$vec
    Mat <-  Q %*% GX.inv.half %*% A
    pred <- KX.test %*% Mat
  }
  if(method=="fgsave"){
    Vfave <- matrix(0, n, n)
    GX.inv <- matpower2(GX,-1,ex)
    GX.inv.half <- matpower2(GX,-.5,ex)
    
    if(yclass) ey=0
    KY.inv <- matpower2(KY,-1,ey)
    for(i in 1:n){ 
      by <- KY[,i]
      Ay <- KX %*% diag(as.vector(KY.inv %*% by)) - KY.inv %*% by %*% t(by) %*% KY.inv %*% KX
      Vfave <- Vfave + (Q- GX.inv %*% Ay %*% GX.inv) %*% GX %*%  (Q- GX.inv %*% Ay %*% GX.inv)
    }
    Vfave = Q %*% GX.inv.half %*% Vfave %*% GX.inv.half %*% Q
    M = Vfave
    A <- eigen(sym(M), sym=TRUE)
    eigval <- A$values
    A <- A$vec
    Mat <- GX.inv.half %*% A
    pred <- KX.test %*% Mat
  }
  out <- list(eval = eigval, pred = pred)
  return(out)
}
gramy=function(y){
  n=length(y)
  k2=y%*%t(y);k1=t(matrix(diag(k2),n,n));k3=t(k1);k=k1-2*k2+k3
  sigma=sum(sqrt(k))/(2*choose(n,2));gamma=1/(2*sigma^2)
  return(exp(-gamma*(k1-2*k2+k3)))
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
###########################################################################
#       function: discrete kernel
###########################################################################
gram.dis=function(y){
  n=length(y);yy=matrix(y,n,n);diff=yy-t(yy);vecker=rep(0,n^2)
  vecker[c(diff)==0]=1;vecker[c(diff)!=0]=0 
  return(matrix(vecker,n,n))}
###########################################################################
#      gram matrix in the x domain, using rkhs gauss kernel 
###########################################################################
gramx=function(x,kt, p=1,gamma=NULL){
  n=dim(x)[1]
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
  return(list(gram=exp(-gamma*(k)),gam=gamma))
}
###########################################################################
#       function: Moor-Panrose power
###########################################################################
mppower=function(matrix,power,ignore){
  eig = eigen(matrix)
  eval = eig$values
  evec = eig$vectors
  m = length(eval[abs(eval)>ignore])
  tmp = evec[,1:m]%*%diag(eval[1:m]^power)%*%t(evec[,1:m])
  return(tmp)
}
###########################################################################
#       function: matpower
###########################################################################
matpower=function(a,alpha){
  a=(a+t(a))/2;tmp=eigen(a)
  return(tmp$vectors%*%diag((tmp$values)^alpha)%*%t(tmp$vectors))}
matpower2 = function(a, alpha, eps)
  matpower(a+eps*onorm(a)*diag(nrow(a)),alpha)
###########################################################################
#       function: operator norm of a matrix
###########################################################################
onorm=function(a) return(eigen(round((a+t(a))/2,8))$values[1])
###########################################################################
#       function: symmetrize a matrix 
###########################################################################
sym=function(a) return(round((a+t(a))/2,7))
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
get.fd = function(x,tt,kern,nbasis,ncv){
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

###############################################################
# make.quad(x)  
##############################################################
make.quad = function(x,inter=FALSE,deg=2){
  n= nrow(x)
  p=ncol(x)
  # add interaction term
  if(inter){
    p.add = p*(p+1)/2
    x.add = matrix(0,n, p.add)
    count = 0
    vars = rep("",p.add)
    for(i in 1:p){
      for(j in i:p){
        count = count + 1
        x.add[,count] = x[,i]*x[,j]
        vars[count] = paste0(i,j)
      }
    }  
  }
  # only polynomial term
  if(!inter){
    p.add = p*deg
    x.add = matrix(0,n, p.add)
    count = 0
    vars = rep("",p.add)
    for(i in 1:p){
      for(j in 2:deg){
        count = count + 1
        x.add[,count] = x[,i]^j
        vars[count] = paste0(i,"^",deg)
      }
    }
    
  }
  out = list(x=cbind(x,x.add),vname=vars)
  return(out)
}

###############################################################
# function: Q = I - J/n
###############################################################
qmat = function(n) return(diag(n)-rep(1,n)%*%t(rep(1,n))/n)