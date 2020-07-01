#------------------------------------------#
# Subcodes (nested-RKHS)   by Jun Song
#------------------------------------------#
# 
# 1. gram: inner product matrix of two random elements 
# 1-1. inprod.rkhs : inner product gram matrix of RKHS objects  (kern = gauss or brown)
# 1-1-1. gramt: gram-matrix for the 1st level space  (function space)
# 1-1-2. gramy: gram-matrix for a scalar Y // categorical Y
# 1-1-3. inprod.rkhs: inner product of RKHS objects
#        norm.rkhs   
# 1-2. gram : gram matrix for the second level space (kernel = gauss or poly)
#--------------------------------------------------------#
# 2. get.fd : get coefficients for RKHS (1st level)


#######################################################################
# 1. gram: inner product matrix of two random elements 
#-------------------------------------------------------------#
## 1-1. inprod.rkhs : inner product gram matrix of RKHS objects----
#-------------------------------------------------------------#
## 1-1-1. gramt: gram-matrix for the 1st level space  (function space)
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

## 1-1-2. gramy: gram-matrix for Y----
gramy=function(y){
  n=length(y)
  k2=y%*%t(y);k1=t(matrix(diag(k2),n,n));k3=t(k1);k=k1-2*k2+k3
  sigma=sum(sqrt(k))/(2*choose(n,2));gamma=1/(2*sigma^2)
  return(exp(-gamma*(k1-2*k2+k3)))
}

gram.dis=function(y){
  n=length(y);yy=matrix(y,n,n);diff=yy-t(yy);vecker=rep(0,n^2)
  vecker[c(diff)==0]=1;vecker[c(diff)!=0]=0 
  return(matrix(vecker,n,n))}

## 1-1-3. inprod.rkhs: inner product of RKHS objects ----
# out : <X_i,Y_j>
# x,y: coefficients of X and Y based on KT     
# dim(x) = (n, nbasis, p)
# is.null(y) => y=x
inprod.rkhs = function(x,kt,y=NULL){
  if(is.null(y)) y = x
  n1=dim(x)[1]; n2=dim(y)[1]
  p=dim(x)[3]
  if(is.na(p)) p=1
  if(p==1){
    xip = x %*% kt %*% t(y)
  }else{
    xip = array(0, c(n1, n2,p))
    for(j in 1:p){
      xip[,,j]=xip[,,j] + x[,,j]%*%kt%*%t(y[,,j])
    }  
  }
  return(xip)
}

## 1-1-3. norm.rkhs: norm of RKHS objects ----
# out : \|X_i-Y_j\|
norm.rkhs = function(x,kt,y=NULL){
  xx = inprod.rkhs(x,kt,x)
  if(is.null(y)){
    xy=xx;yy=xx
  }else{
    xy = inprod.rkhs(x,kt,y)
    yy = inprod.rkhs(y,kt,y) 
  }
  n1 = nrow(xx); n2=nrow(yy)
  xxrow = matrix(diag(xx),n1,n2)
  yycol = matrix(diag(yy),n1,n2,byrow=T)
  
  normxy = xxrow + yycol - 2*xy
  normxy[normxy<0]=0
  return(normxy)
}
# get norm from inner-product
norm.rkhs.xip = function(xip){
  if(is.matrix(xip)){
    k = t(-2 * xip + diag(xip)) + diag(xip)
    k[k<0]=0
  }else{  # xip is array of n X n x p
    n=dim(xip)[1]; p=dim(xip)[3]
    k = array(0,c(n,n,p))
    for(j in 1:p){
      k[,,j] =  t(-2 * xip[,,j] + diag(xip[,,j])) + diag(xip[,,j])  
    }
  }
  return(k)
}
gam = function(k){
  k[k<0]=0
  sigma=sum(sqrt(k))/(2*choose(n,2));
  gamma=1/(2*sigma^2)
  return(gamma)
}


#-------------------------------------------------------------#
## 1-2. gram : gram matrix for the second level space----
#-------------------------------------------------------------#
# kernel: gauss or poly
# parameters: gamma (for gauss), d (for poly)
gram = function(xip, kernel="gauss", gamma=NULL, d=3){
  if(kernel=="gauss"){
    normx = norm.rkhs.x(xip)
    if(is.null(gamma)){
      gamma=gam(normx)
    }
    out = exp(-gamma*(k))
  }else if(kernel=="poly"){
    out = (1+xip)^d
  }else{
    stop('kernel can be gauss or poly')
  }
  return(out)
}


#-------------------------------------------------------------#
## 2. get.fd.rkhs : get coefficients for RKHS----
#-------------------------------------------------------------#
# output : coeff, kt
get.fd.rkhs = function(x, tt, kern="gauss", et=0.05, tune=TRUE, ncv1=10){
  if(tune){
    epset=exp(seq(log(10^(-10)), log(10^(2)), len=ncv1))
    
    ttt=matrix(seq(from=0,to=1,length=nt),n,nt,byrow=T)
    out.gcvt=numeric();for(i in 1:ncv1) out.gcvt=c(out.gcvt,gcvt(x,tt,epset[i],kern))
    et=epset[which.min(out.gcvt)]
  }
  nt=dim(x)[2];n=dim(x)[1];Q=diag(n)-rep(1,n)%*%t(rep(1,n))/n
  
  kt = gramt(tt, tt, kern)
  scale = eigen(kt,sym=TRUE)$value[1]
  ktinv = matpower(kt+scale*et*diag(nrow(kt)),-1)

  p=dim(x)[3]
  if(is.na(p)) p=1
  if(p==1){
    xcoef = x %*% ktinv
    xeval = xcoef %*% kt
  }else if(p>1){
    xcoef = array(0,c(n,nt,p))
    xeval = array(0,c(n,nt,p))
    for(j in 1:p){
      xcoef[,,j] = x[,,j] %*% ktinv 
      xeval[,,j] = xcoef[,,j] %*% kt
    }  
  }
  return(list(coef=xcoef, kt=kt, kern=kern, et=et))
}
#-------------------------------------------------------------#
#   gcvt:   gcv tuning for et
#-------------------------------------------------------------#
gcvt=function(x,tt,et,kern){  
  n=dim(x)[1];nt=dim(x)[2] 
  px=dim(x)[3]
  
  kt = gramt(tt, tt, kern)
  scale = eigen(kt,sym=TRUE)$value[1]
  ktinv = matpower(kt+scale*et*diag(nrow(kt)),-1)
  
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