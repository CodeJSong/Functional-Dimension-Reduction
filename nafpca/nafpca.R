source('nafpca_subcodes.R')
source('../subcodes_general.R')
source('../subcodes_ftn_estimation.R')
source('../subcodes_nested_RKHS.R')
library(RcppArmadillo)
Rcpp::sourceCpp('nafpca_cpp_sub.cpp')
library(fda)
library(doParallel)
library(foreach)
library(rARPACK)  # faster eigen-decomposition

cl <- parallel::makeCluster(detectCores(logical = TRUE))
registerDoParallel(cl)
#on.exit(stopCluster(cl))

#####################################################
# Additive Functional PCA
#####################################################
# x : n X nt X p   (if p=1, n X nt matrix, else, array)
# t : time point vector
# p : dimension of functions

nafpca = function(x=NULL, tt=NULL, p=1,type="1dsup", unbalanced=FALSE,
                 basisname="bspline", nbasis=31,
                 c=1, d=3, kernel="gaussian",
                 ncv1=10, ncv2=15,  gamma.tune=TRUE, poly.tune=TRUE,
                 m1=10,m2=10,
                  shx=2,ex=0.05){
  start.time = Sys.time()
  if(!is.null(dim(x))){ # is.null is required if data are a list-class
    p = dim(x)[3]
    if(is.na(p)) p = 1  
  }
  if(type!="2dsup"){
    
    # get functional object
    # 1. basis expansion approach
    if(basisname=='bspline'|basisname=='fourier'){
      if(is.matrix(x)){
        x=t(x)   # nt  X n
      }else if(is.array(x)) x = aperm(x, c(2,1,3))  # nt X n X p
      if(!unbalanced) temp=get.fd(xraw=x,tt=tt,basisname=basisname,nbasis=nbasis,ncv1)
      if(unbalanced) temp=get.fd.sparse(x, tt, basisname, nbasis, ncv1)  
      
      
      xcoef = temp$coef;basis = temp$basis
      
      if(is.matrix(xcoef)){
        xcoef=t(xcoef)   # n X nbasis
      }else if(is.array(xcoef)) xcoef = aperm(xcoef, c(2,1,3))  # n X nbasis X p
      
      
      if(basisname=="bspline") GB = bsplinepen(basis,0)
      if(basisname=="fourier") GB = fourierpen(basis,0)
      
      xip = inprod.rkhs(xcoef,GB)

    }else{
      #2. rkhs based approach
      temp = get.fd.rkhs(x, tt, kern=basisname, et=0.05, tune=TRUE, ncv1=10)
      xcoef = temp$coef; kt = temp$kt
      xip = inprod.rkhs(xcoef,kt)
    }
    n = dim(xcoef)[1]
    gx=NULL
    cv.shx = NULL
    cv.cd=NULL
    

    
    if(p==1){
      
      one = matrix(1, n, 1)
      Q = diag(n) - one %*% t(one) / n
      
      if(kernel=="gaussian"){
        
        normX = norm.rkhs.xip(xip)
        gx = gam(normX)
        if(gamma.tune){
          shx.grid = c(exp(seq(log(10^(-10)),log(10^(2)),len=ncv2)))
          cv.shx =  foreach(i=1:ncv2, .combine='c', .export=c("cv_gamma", "matpower", "Gn")) %dopar%
            cv_gamma(gx, shx.grid, normX, Q)
          shx = shx.grid[which.max(cv.shx)]
        }
        gx = gx*shx
        KX = exp(-gx*normX)
      }else if(kernel=="poly"){
        if(poly.tune){
          d.grid = 1:3
          c.grid = 1
          cd.grid = expand.grid(c.grid, d.grid)
          ncv2 = nrow(cd.grid)
          cv.cd =  foreach(i=1:ncv2, .combine='c', .export=c("cv_poly", "matpower", "Gn")) %dopar%
            cv_poly(cd.grid, xip, Q)
          c = cd.grid[which.max(cv.cd),1]
          d = cd.grid[which.max(cv.cd),2]
        }
        KX = (c+xip)^d
      }
      
      GX = Q %*% KX %*% Q
      #GX.e = GX + ex*as.numeric(eigen_val(GX)) *diag(n)
      GX.inv.half = mppower(GX,-0.5)
      
      eigenM = eigen(GX, sym=TRUE)
      Mat = eigenM$vec
      pred = t(GX)%*% GX.inv.half %*% Mat
    }else if(p>1){
      
      one = matrix(1, n, 1)
      Q = diag(n) - one %*% t(one) / n
      
      if(kernel=="gaussian"){
        normxs = array(0,c(n,n,p))
        gxs = rep(0,p)
        
        for(j in 1:p){
          normxs[,,j] =  t(-2 * xip[,,j] + diag(xip[,,j])) + diag(xip[,,j])  
          gxs[j] = gam(normxs[,,j])
        }
        
        if(gamma.tune & p==2){
          ncv2=min(15,ncv2)
          shx.grid = c(exp(seq(log(10^(-10)),log(10^(2)),len=ncv2)))
          shx.grid = expand.grid(shx.grid,shx.grid)
          ## c function cannot be exported... to foreach
          cv.shx =  foreach(i=1:(ncv2^2), .combine='c', .export=c("cv_gamma.mult", "matpower", "Gn")) %dopar%
            cv_gamma.mult(gxs, shx.grid, normxs, n, p, Q)
          
          shxs = shx.grid[which.max(cv.shx),]
          
          shx=shxs
        }else if(gamma.tune & (p>2)){
          cat('tuning gamma will take so long for p>=3. individual tuning--')
          shx.grid = c(exp(seq(log(10^(-10)),log(10^(2)),len=ncv2)))
          shxs = rep(2,p)
          for(j in 1:p){
            normX = normxs[,,j]
            gx = gxs[j]
            cv.shx =  foreach(i=1:ncv2, .combine='c', .export=c("cv_gamma", "matpower", "Gn")) %dopar%
              cv_gamma(gx, shx.grid, normX, Q)  
            shx = shx.grid[which.max(cv.shx)]
            shxs[j] = shx
          }
          
          
        }
        
        KXs = matrix(0, n*p, n*p)
        D.inv.half = matrix(0, n*p, n*p)
        GXs = matrix(0, n*p, n)
        gxs = rep(0,2)
        for(j in 1:p){
          normxs[,,j] =  t(-2 * xip[,,j] + diag(xip[,,j])) + diag(xip[,,j])  
          gxs[j] = as.numeric(gam(normxs[,,j])*shxs[j])
          KXj = exp(-gxs[j]*normxs[,,j]); KXs[((j-1)*n+1):(j*n), ((j-1)*n+1):(j*n)] = KXj
          GXj = Q %*% KXj %*% Q
          D.inv.half[((j-1)*n+1):(j*n), ((j-1)*n+1):(j*n)] = mppower(GXj,-0.5);
          GXs[((j-1)*n+1):(j*n),]= GXj
        }
        gx = gxs
      }else if(kernel=="poly"){
        KXs = matrix(0, n*p, n*p)
        D.inv.half = matrix(0, n*p, n*p)
        GXs = matrix(0, n*p, n)
        for(j in 1:p){
          KXj = (c+xip[,,j])^d
          GXj = Q %*% KXj %*% Q
          D.inv.half[((j-1)*n+1):(j*n), ((j-1)*n+1):(j*n)] = mppower(GXj,-0.5);
          GXs[((j-1)*n+1):(j*n),]= GXj
        }
      }
      
      M = GXs %*% t(GXs)
      eigenM = eigen(D.inv.half %*% M %*% D.inv.half, sym=TRUE)
      Mat = eigenM$vec
      pred = t(GXs) %*% D.inv.half %*% Mat
      GX.inv.half = D.inv.half
    }
  }else if(type=="2dsup"){
    temp=NULL
    # input : vectorized   x = n * nt^2
    n = dim(x)[1]
    one = matrix(1, n, 1)
    Q = diag(n) - one %*% t(one) / n
    tmp=dim(x)[2]; n1 = n2 = sqrt(tmp[1])
    
    U = fourier.basis(1:n1/n1, m1)
    V = fourier.basis(1:n2/n2, m2)
    m1 = dim(U)[1]; m2= dim(V)[1]   # redefine m1, m2 just in case pre-defiend m1 and m2 are even numbers.
    
    UU.inv = mppower(U %*% t(U), -1)
    VV.inv = mppower(V %*% t(V), -1)
    
    coefmat =  kronecker(VV.inv %*% V, UU.inv %*% U)
    normX = normmat.tensor(x, U,V, UU.inv, VV.inv)
    
    gx = gam(normX)*shx
    KX = exp(-gx*normX)
    GX = Q %*% KX %*% Q
    GX.inv.half = mppower(GX,-0.5, 0)
    
    eigenM = eigen(GX, sym=TRUE)
    Mat = eigenM$vec
    pred = t(GX)%*% GX.inv.half %*% Mat
  }
  
  comp.time = Sys.time()-start.time
  out = list(pred=pred, comp.time=comp.time, eval = eigenM$val, dim=which.max(Gn(eigenM$val,n)), ftn=temp, gx=gx, GX.inv.half = GX.inv.half, Mat=Mat,
             cv.shx=cv.shx, shx=shx, c=c,d=d, cv.cd=cv.cd)
  return(out)
}
# cv_gamma = function(gam, shx.grid, normx, Q){
#   gx = gam*shx.grid[i];
#   KX = exp(-gx*normx);
#   GX = Q %*% KX %*% Q;
#   evals = eigen(GX,sym=TRUE)$val
#   return(evals[1])
# }
cv_poly = function(cd.grid, xip, Q){
  c = cd.grid[i,1]
  d = cd.grid[i,2]
  KX = (c+xip)^d
  GX = Q %*% KX %*% Q;
  eigenM = eigen(GX, sym=TRUE)
  Mat = eigenM$vec
  GX.inv.half = matpower(GX,-0.5)
  pred = t(GX)%*% GX.inv.half %*% Mat
  
  ev = eigenM$val
  n = nrow(KX)
  d = which.max(Gn(ev,n))
  if(d==1){
    return(var(pred[,1]))
  }else{
    return(sum(diag(var(pred[,1:d]))))  
  }
  #return(var(pred[,1]))
}
cv_gamma = function(gam, shx.grid, normx, Q){
  gx = gam*shx.grid[i];
  KX = exp(-gx*normx);
  GX = Q %*% KX %*% Q;
  eigenM = eigen(GX, sym=TRUE)
  Mat = eigenM$vec
  GX.inv.half = matpower(GX,-0.5)
  pred = t(GX)%*% GX.inv.half %*% Mat
  
  ev = eigenM$val
  n = nrow(KX)
  d = which.max(Gn(ev,n))
  if(d==1){
    return(var(pred[,1]))
  }else{
    return(sum(diag(var(pred[,1:d]))))  
  }
  #return(var(pred[,1]))
}
cv_gamma.mult = function(gxs, shx.grid, normxs, n, p, Q){
  KXs = matrix(0, n*p, n*p)
  D.inv.half = matrix(0, n*p, n*p)
  GXs = matrix(0, n*p, n)
  shxs = as.vector(shx.grid[i,])
  
  for(j in 1:p){
    gxs[j] = as.numeric(gxs[j]*shxs[j])
    KXj = exp(-gxs[j]*normxs[,,j]); KXs[((j-1)*n+1):(j*n), ((j-1)*n+1):(j*n)] = KXj
    GXj = Q %*% KXj %*% Q
    D.inv.half[((j-1)*n+1):(j*n), ((j-1)*n+1):(j*n)] = matpower(GXj,-0.5);
    GXs[((j-1)*n+1):(j*n),]= GXj
  }
  
  M = GXs %*% t(GXs)
  eigenM = eigen(D.inv.half %*% M %*% D.inv.half, sym=TRUE)
  Mat = eigenM$vec
  pred = t(GXs) %*% D.inv.half %*% Mat
  
  ev = eigenM$val
  d = which.max(Gn(ev,n))
  if(d==1){
    return(var(pred[,1]))
  }else{
    return(sum(diag(var(pred[,1:d]))))  
  }
  #return(var(pred[,1]))
}
predict.fapca = function(pcobj, x, tt){
  temp = pcobj$ftn
  n = temp$n

  one = matrix(1, n, 1)
  Q = diag(n) - one %*% t(one) / n
  xcoef = temp$coef
  xip = temp$xip;
  normX = temp$normx;
  GB = temp$binnerprod;
  basis = temp$basis
  gx = pcobj$gx
  GX.inv.half = pcobj$GX.inv.half
  Mat = pcobj$Mat

  bname = basis$type
  nbasis = basis$nbasis
  p = temp$p

  temp = get.ftn(x, tt, p=p, basisname=bname, nbasis=nbasis)
  test.coef = temp$coef

  if(p==1){
    x1ip <- t(xcoef) %*% GB %*% xcoef
    x2ip <- t(test.coef) %*% GB %*% test.coef
    x12ip <- t(xcoef) %*% GB %*% test.coef
    normX = t(diag(x2ip)- t(x12ip)) + diag(x1ip)-  x12ip
    KX1X2 = exp(-gx * normX)

    pred = t(KX1X2) %*% Q%*% GX.inv.half %*% Mat
  }else{
    stop('need coding')
  }
  return(pred)
}



###########################################################################
#    fpca: functional PCA 
###########################################################################
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
  B.half = mppower(GB,0.5)
  
  if(p==1){
    Sigma = B.half %*% xcoef %*% Q %*% t(xcoef) %*% B.half / n
    egn = eigen(Sigma, sym=TRUE)
    B.inv.half = mppower(GB,-0.5)
    pred = Q %*% t(xcoef) %*% GB %*% B.inv.half %*% egn$vec
    out = list(pred=pred, eval=egn$val, mat =  B.inv.half %*% egn$vec)
  }else if(p>1){ # BX is the same for now...
    M.half = B.half %*% xcoef[,,1] %*% Q
    B.inv.half = mppower(GB,-0.5)
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
    out = list(pred=pred, eval=egn$val, mat=D.inv.half %*% egn$vec)
  }
  return(out)
}


fpca.predict = function(fpcobj, ftn, x, tt){
  # p==1
  temp = ftn
  p = temp$p
  n = temp$n
  xcoef = temp$coef
  xip = temp$xip;
  normX = temp$normx;
  GB = temp$binnerprod;
  basis = temp$basis
  lambdaout = temp$lambdaout
  
  bname = basis$type
  nbasis = basis$nbasis
  
  mat = fpcobj$mat
  
  temp = get.ftn(x, tt, p=p, basisname=bname, nbasis=nbasis)
  test.coef = temp$coef
  n = ncol(test.coef)
  one = matrix(1, n, 1)
  Q = diag(n) - one %*% t(one) / n
  pred = Q %*% t(test.coef) %*% GB %*% mat
  return(pred)
}


kfpca.predict = function(xtest, tttest, xtrain,tttrain, train, nharm=10, ncv=15){  # only for balanced observed data/ 1 dimension
  # train : kfpca output - kfpca object
  ### when xtequal==TRUE
  x2=xtest
  ttt=sort(unique(as.numeric(unlist(tttest))))
  n2 = length(x2)
  x=xtrain
  tx = tttrain
  n1=length(x)
  n=n1+n2
  for(i in (n1+1):n){
    tx[[i]] = ttt
    x[[i]] = x2[[i-n1]]
  }
  xt.unique = as.vector(sort(unique(unlist(tx))))
  ntx = length(xt.unique)
  basis = train$basis ; traincoef = train$xcoef;
  X = eval.basis(xt.unique, basis)
  binnerprod = bsplinepen(basis,0)
  RR=bsplinepen(basis, 2)
  lambda_all= exp(seq(log(10^(-6)),log(10^2),length=ncv))

  xmat = makematrix(x,ntx)
  gcvf.mat = bsplinegcv_mat(xmat=X, rr=RR, ystar=t(as.matrix(xmat)), lambdagrid=lambda_all)
  xcoef = gcvf.mat$mucoef

  xip = t(xcoef) %*% binnerprod %*% xcoef
  normX = t(-2 * xip + diag(xip)) + diag(xip)
  gx = gam(normX)/train$shx
  KX = exp(-gx*normX)
  kfpca.pred = KX[(n1+1):n, 1:n1] %*% train$KX.inv.half %*% train$Mat

  x2coef = xcoef[,(n1+1):n]
  real2.fd = fd(x2coef,basis)
  fpca.pred = matrix(0,n2,nharm)
  for(j in 1:nharm){
    fpca.pred[,j] = inprod(train$fpca$harmonics[j], center.fd(real2.fd))
  }
  return(list(kfpca.pred=kfpca.pred, fpca.pred = fpca.pred))
}
cvx.ftn = function(normX, y, ncv2, M,Q,standard="lda", randomsearch=FALSE ){
  gx0 = gam(normX)   # make it default value
  nclass = length(unique(y)); n=length(y)
  nj = rep(0,nclass)
  for(j in 1:nclass){
    nj[j] = sum(y==j)
  }
  cvwinx= c(1,exp(seq(log(0.1),log(10),len=ncv2-1)))
  cvex =  c(n^(-1/4),exp(seq(log(n^(-1/4)/10),log(n^(-1/4)*10),len=ncv2-1)))
  if(randomsearch){
    cvwinx = runif(ncv2, 0.1,10)
  }
  cvgx = rep(0,ncv2); cvexv = rep(0,ncv2)
  for(i in 1:ncv2){
    gx = gx0 / cvwinx[i]
    KX = exp(-gx*normX)

    ## Estimate central class
    KX.inv.half = matpower2(KX,-0.5,ex)
    KX.half = matpower(KX,1/2)
    M = KX.inv.half %*% KX %*% Q %*% KX %*% KX.inv.half
    #       M = KX.half %*% Q %*% KX.half
    M = (M + t(M))/2
    eigenM = eigen(M, sym=TRUE)
    Mat = eigenM$vec
    elam = eigenM$val
    pred = KX %*% KX.inv.half %*% Mat
    d = which(cumsum(elam)/sum(elam) >0.999)[1]
    if(standard=="lda"){
      if(d>10) d=10
      preded = NULL
      for(k in 1:n){
        temp = pred[-k,1:d] ; yy = y[-k]
        preded = c(preded,predict(lda(temp, yy), pred[k,1:d])$class)
      }
      cvgx[i] = -sum(preded==y)
    }
    if(standard=="knn"){
      temp=NULL
      for(j in 1:n) temp=c(temp,knn(pred[-j,1:d],pred[j,1:d],y[-j]))
      cvgx[i] = -sum(temp==y)
    }
    if(standard=="condvar"){
      temp = pred[,1:d] %*% mppower(var(pred[,1:d]),-0.5)  ; yy =y
      v=0
      for(k in 1:n){
        temp2 = temp[-k,];
        for(j in 1:nclass){
          nj[j] = sum(yy[-k]==j)
        }
        for(j in 1:nclass){
          xj = temp2[yy[-k]==j,]
          v = v+ var(xj)*nj[j]
        }
        cvgx[i] = cvgx[i] + sum(diag(v))
      }
    }
    if(standard=="qda"){
      if(d>10) d=10
      preded = NULL
      for(k in 1:n){
        temp = pred[-k,1:d] ; yy = y[-k]
        preded = c(preded,predict(qda(temp, yy), pred[k,1:d])$class)
      }
      cvgx[i] = -sum(preded==y)
    }

  }
  shx = cvwinx[which.min(cvgx)]
  gx = gx0 / shx
  KX = exp(-gx*normX)
  return(list(KX=KX, shx=shx, cvgx=cvgx, gx=gx))
}


fpca.ex = function(xcoef,basis){
  xcoef = t(center(t(temp$xcoef)))
  B = binnerprod
  Sigma=mppower(B, 0.5) %*% xcoef %*% t(xcoef) %*% mppower(B, 0.5)/n
  FD=fd( mppower(B,-0.5) %*% eigen((DD+t(DD)/2),sym=TRUE)$vec, basis) ## functional object b_1(t) ...
  scores = t(xcoef) %*% mppower(B,0.5) %*% eigen((DD+t(DD)/2),sym=TRUE)$vec[,1:2]
}


Gn=function(evals, n){
  maxd=length(evals)
  cumsum(evals) - (1:maxd)*log(n+1)
}
