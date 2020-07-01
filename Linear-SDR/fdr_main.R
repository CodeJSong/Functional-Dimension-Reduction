#############################################
# functional dimension reduction
#############################################
# Linear  (WIRE, WAVE, WDR)
# Nonlinear (fgsir fgsave)
# FPCA
# FSIR FSAVE
#############################################
###########################################################################
#      linear functional sdr
# time domain is set to be "nt" equally spaced [0,1] 
###########################################################################
fsdr=function(x=x.tra,y=y.tra,kern,et=5*10^(-6),ex=0.05,ey=0.05,method,ytype,x.tes=NULL){
  n=dim(x)[1];nt=dim(x)[2]
  px=dim(x)[3]
  if(is.na(px)) px=1
  if(is.null(x.tes))x.tes=x
  tto=seq(from=0,to=1,length=nt);tte=seq(from=0,to=1,
                                         length=2*(nt-1)+1);ttt=t(matrix(tto,nt,n))
  n.test = dim(x.tes)[1];nt.test=dim(x.tes)[2]
  ttt.test=t(matrix(tto,nt.test,n.test))
  xcoo=xcoomat(x,tte,ttt,et,kern)
  x.test.coo=xcoomat(x.tes,tte,ttt.test,et,kern)
  
  kt=gramt(tto,tto,kern);kt.half=mppower(kt,1/2,10^(-8));Q=qmat(n)
  kt.nhalf=matpower(kt+et*onorm(kt)*diag(nt),-1/2)
  if(ytype=="scalar") ky=gramy(y)
  if(ytype=="function") {ycoo=xcoomat(y,tte,ttt,et,kern);ky=gramx(ycoo,kt)$gram}
  if(ytype=="categorical") ky=gram.dis(y)
  if(ytype=="scalar"|ytype=="function")
  {tau=(1/n)*matpower((1/n)*ky+ey*onorm((1/n)*ky)*diag(n),-1)%*%ky}
  if(ytype=="categorical") tau=mppower(ky,-1,10^(-7))%*%ky
  if(px==1) {
    tmp=(1/n)*kt.half%*%t(xcoo)%*%Q%*%xcoo%*%kt.half
    ax=matpower(tmp+ex*onorm(tmp)*diag(nt),-1/2)

  } 
  if(px>1){
    
    j=1
    tmp=(1/n)*kt.half%*%t(xcoo[,,j])%*%Q%*%xcoo[,,j]%*%kt.half
    axj = mppower2(tmp+ex*onorm(tmp)*diag(nt),-0.5)
    ax = axj
    for(j in 2:px){
      tmp=(1/n)*kt.half%*%t(xcoo[,,j])%*%Q%*%xcoo[,,j]%*%kt.half
      axj = mppower2(tmp+ex*onorm(tmp)*diag(nt),-0.5)
      ax = as.matrix(bdiag(ax, axj))
    }
    kt.nhalf = diag.mat(kt.nhalf,px)
    kt.half = diag.mat(kt.half,px)
    kt = diag.mat(kt,px)
    
    xcoo.arr = xcoo
    x.test.coo.arr = x.test.coo
    xcoo = matrix(0,n,nt*px)
    x.test.coo = matrix(0,n.test,nt.test*px)
    for(i in 1:n){
      xcoo[i,] = as.vector(xcoo.arr[i,,])
    }
    for(i in 1:n.test){
      x.test.coo[i,] = as.vector(x.test.coo.arr[i,,])
    }
    nt = nt*px
  }
  lcoo=kt.nhalf%*%ax%*%kt.half%*%t(xcoo)%*%tau
  taubar=apply(tau,2,mean)
  qcoo=array(0,c(nt,nt,n))
  for(i in 1:n){
    qcoo[,,i]=kt.nhalf%*%ax%*%kt.half%*%(t(xcoo)%*%diag(taubar[i]*rep(1,n)-tau[,i])%*%
                                           xcoo)%*%kt.half%*%ax%*%kt.half}
  if(method=="wire"){
    if(ytype=="scalar"|ytype=="function")
    {bxy=sym(t(xcoo)%*%Q%*%((1/n)*ky)%*%
               matpower((1/n)*ky+ey*onorm((1/n)*ky)*diag(n),-1)%*%Q%*%(xcoo))}
    if(ytype=="categorical")
    {bxy=sym(t(xcoo)%*%Q%*%((1/n)*ky)%*%mppower((1/n)*ky,-1,10^(-7))%*%Q%*%(xcoo))}
    v=eigen(sym(ax%*%kt.half%*%bxy%*%kt.half%*%ax))$vectors
    #pred=t(t(v)%*%kt.nhalf%*%kt.half%*%ax%*%kt.nhalf%*%kt%*%t(xcoo))
   
    pred.train=t(t(v)%*%kt.nhalf%*%kt.half%*%ax%*%kt.nhalf%*%kt%*%t(xcoo))
    pred.test=t(t(v)%*%kt.nhalf%*%kt.half%*%ax%*%kt.nhalf%*%kt%*%t(x.test.coo))
    return(list(pred.train=pred.train, pred.test=pred.test))}
  if(method=="wave"){
    swave1=0;swave2=0;swave3=0;swave4=0;
    for(i in 1:n){
      tmp=lcoo[,i]%*%t(lcoo[,i])%*%kt;swave1=swave1+tmp%*%tmp/n
      swave2=swave2+lcoo[,i]%*%t(lcoo[,i])%*%kt%*%qcoo[,,i]/n
      swave3=swave3+qcoo[,,i]%*%lcoo[,i]%*%t(lcoo[,i])%*%kt/n
      swave4=swave4+qcoo[,,i]%*%qcoo[,,i]/n};
    swave=swave1+swave2+swave3+swave4;smethod=swave}
  if(method=="wdr"){
    swdr1=0;swdr2=0;swdr3=0;swdr4=0
    for(i in 1:n){
      swdr1=swdr1+qcoo[,,i]%*%qcoo[,,i]/n
      swdr2=swdr2+qcoo[,,i]/n
      swdr3=swdr3+lcoo[,i]%*%t(lcoo[,i])%*%kt/n
      swdr4=swdr4+c(t(lcoo[,i])%*%kt%*%lcoo[,i]/n)}
    swdr=swdr1+swdr2%*%swdr2+swdr3%*%swdr3+swdr4*swdr3;smethod=swdr}
  v=eigen(sym(kt.nhalf%*%(kt%*%smethod)%*%kt.nhalf))$vectors;fcoo=kt.nhalf%*%v
  #pred=matrix(0,n,nt)
  #for(i in 1:nt) pred[,i]=c(t(fcoo[,i])%*%kt.half%*%ax%*%kt.half%*%t(xcoo))
  pred.test=matrix(0,n.test,nt)
  pred.train=matrix(0,n,nt)
  for(i in 1:nt) {
    pred.train[,i]=c(t(fcoo[,i])%*%kt.half%*%ax%*%kt.half%*%t(xcoo))
    pred.test[,i]=c(t(fcoo[,i])%*%kt.half%*%ax%*%kt.half%*%t(x.test.coo))
  }
  return(list(pred.train=pred.train, pred.test=pred.test))
}
###########################################################################
#       function: Tuning ftn
###########################################################################
tune.fsdr = function(x,y,n,nt,kern,ytype,ncv=6){
  
  epset=exp(seq(log(10^(-6)), log(10^(-1)), len=ncv))
  etset=exset=eyset=epset
  
  ttt=matrix(seq(from=0,to=1,length=nt),n,nt,byrow=T)
  out.gcvt=numeric();for(i in 1:ncv) out.gcvt=c(out.gcvt,gcvt(x,ttt,etset[i],kern))
  et=minimizer(etset,out.gcvt)
  out.gcvx=numeric();for(i in 1:ncv) out.gcvx=c(out.gcvx,gcvxy(x,y,exset[i],"ex",kern,et,ytype))
  ex=minimizer(exset,out.gcvx)
  out.gcvy=numeric();for(i in 1:ncv) out.gcvy=c(out.gcvy,gcvxy(x,y,eyset[i],"ey",kern,et,ytype))
  ey=minimizer(eyset,out.gcvy)
  return(list(et=et,ex=ex,ey=ey))
}
###############################################################
# linear fsdr with fd object (fda package)  function on function
###############################################################
fsdr.fd = function(x,y,tt=NULL,et=0.05,ex=0.05,ey=0.05,method,ytype,x.tes=NULL, tune=TRUE,ncv=10){
  require(fda)
  if(is.null(x.tes)) x.tes = x
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
  }
  if(is.fd(y)){
    py=dim(y$coef)[3]
    if(is.na(py)) py=1
    basis.name = y$basis$type
    nbasis = y$basis$nbasis
    ycoo = (y$coef)
  }else if(!is.fd(y)){
    py=0
    ycoo=y
    if(ytype=="scalar") ky=gramy(y)
    if(ytype=="categorical") ky=gram.dis(y)
  }
  
  if(py>1){
    ycoo.arr = ycoo
    ycoo = matrix(0, n, py*nbasis)
    for(i in 1:n){
      ycoo[i,] = as.vector(ycoo.arr[,i,])
    }
    ky=gramx(ycoo,kt,py)$gram
  }else if(py==1){
    ycoo = t(ycoo) 
    ky=gramx(ycoo,kt,py)$gram
  }
  
  
  if(px>1){
    xcoo.arr = xcoo ; x.test.coo.arr = x.test.coo
    xcoo = matrix(0, n, px*nbasis)
    x.test.coo = matrix(0, n.test, px*nbasis)
    for(i in 1:n){
      xcoo[i,] = as.vector(xcoo.arr[,i,])
    }
    for(i in 1:n.test){
      x.test.coo[i,] = as.vector(x.test.coo.arr[,i,])
    }
  }else if(px==1) {
    xcoo = t(xcoo)
    x.test.coo = t(x.test.coo)
  }
  kx=gramx(xcoo,kt,px)$gram
  
  # Tuning
  if(tune==TRUE){
    
    epset=exp(seq(log(10^(-6)), log(10^(-1)), len=ncv))
    
    
    
    #out.gcvt=numeric();for(i in 1:6) out.gcvt=c(out.gcvt,gcvt.fd(n,xcoo,kt,epset[i]))
    #et=minimizer(epset,out.gcvt)
    
    kt.half=mppower2(kt,0.5);Q=qmat(n)
    kt.nhalf=mppower2(kt+et*onorm(kt)*diag(nbasis),-0.5)
    
    out.gcvx=numeric();for(i in 1:ncv) out.gcvx=c(out.gcvx, gcvxy.fd2(n,kx,ky,epset[i],'ex',kt,kt.half,kt.nhalf, xcoo, Q))
    ex=minimizer(epset,out.gcvx)
    out.gcvy=numeric();for(i in 1:ncv) out.gcvy=c(out.gcvy, gcvxy.fd2(n,kx,ky,epset[i],'ey',kt,kt.half,kt.nhalf,xcoo, Q))
    ey=minimizer(epset,out.gcvy)
  }
  
  
  
  kt.half=mppower2(kt,0.5);Q=qmat(n)
  kt.nhalf=mppower2(kt+et*onorm(kt)*diag(nbasis),-0.5)
  
  # get ax
  if(px>1){
    j=1
    tmp=(1/n)*kt.half%*%xcoo.arr[,,j]%*%Q%*%t(xcoo.arr[,,j])%*%kt.half
    axj = mppower2(tmp+ex*onorm(tmp)*diag(nbasis),-0.5)
    ax = axj
    for(j in 2:px){
      tmp=(1/n)*kt.half%*%xcoo.arr[,,j]%*%Q%*%t(xcoo.arr[,,j])%*%kt.half
      axj = mppower2(tmp+ex*onorm(tmp)*diag(nbasis),-0.5)
      ax = as.matrix(bdiag(ax, axj))
    }
    kt.nhalf = diag.mat(kt.nhalf,px)
    kt.half = diag.mat(kt.half,px)
    kt = diag.mat(kt,px)
    nbasis = nbasis*px
  }else if(px==1) {
    tmp=(1/n)*kt.half%*%t(xcoo)%*%Q%*%xcoo%*%kt.half
    ax=mppower2(tmp+ex*onorm(tmp)*diag(nbasis),-0.5)
  }
  
  tau=(1/n)*mppower2((1/n)*ky+ey*onorm((1/n)*ky)*diag(n),-1)%*%ky
  lcoo=kt.nhalf%*%ax%*%kt.half%*%t(xcoo)%*%tau
  taubar=apply(tau,2,mean)
  qcoo=array(0,c(nbasis,nbasis,n))
  
  for(i in 1:n){
    qcoo[,,i]=kt.nhalf%*%ax%*%kt.half%*%(t(xcoo)%*%diag(taubar[i]*rep(1,n)-tau[,i])%*%
                                           xcoo)%*%kt.half%*%ax%*%kt.half}
  if(method=="wire"){
    if(ytype=="scalar"|ytype=="function")
    {bxy=t(xcoo)%*%Q%*%((1/n)*ky)%*%
      mppower2((1/n)*ky+ey*onorm((1/n)*ky)*diag(n),-1)%*%Q%*%(xcoo)}
    if(ytype=="categorical")
    {bxy=t(xcoo)%*%Q%*%((1/n)*ky)%*%mppower2((1/n)*ky)%*%Q%*%(xcoo)}
    v = svd(sym(ax%*%kt.half%*%bxy%*%kt.half%*%ax))$u
    #v=eigen(sym(ax%*%kt.half%*%bxy%*%kt.half%*%ax))$vectors
    #pred=t(t(v)%*%kt.nhalf%*%kt.half%*%ax%*%kt.nhalf%*%kt%*%t(xcoo))

    pred = xcoo %*% kt %*% kt.nhalf %*% ax %*% kt.half %*% kt.nhalf %*% v
    pred.test = x.test.coo %*% kt %*% kt.nhalf %*% ax %*% kt.half %*% kt.nhalf %*% v
    if(px>1){
      nt = nt/px
      tmp = matrix(0,n,nt)
      for(j in 1:px){
        tmp = tmp+pred.train[,((j-1)*nt+1):(j*nt) ]
      }
      pred.train = tmp
      tmp = matrix(0,n,nt.test)
      for(j in 1:px){
        tmp = tmp+pred.test[,((j-1)*nt.test+1):(j*nt.test) ]
      }
      pred.test = tmp
    }
    return(list(pred=pred, pred.test=pred.test, et=et,ex=ex,ey=ey,kx=kx))}
  if(method=="wave"){
    swave1=0;swave2=0;swave3=0;swave4=0;
    for(i in 1:n){
      tmp=lcoo[,i]%*%t(lcoo[,i])%*%kt;swave1=swave1+tmp%*%tmp/n
      swave2=swave2+lcoo[,i]%*%t(lcoo[,i])%*%kt%*%qcoo[,,i]/n
      swave3=swave3+qcoo[,,i]%*%lcoo[,i]%*%t(lcoo[,i])%*%kt/n
      swave4=swave4+qcoo[,,i]%*%qcoo[,,i]/n};
    swave=swave1+swave2+swave3+swave4;smethod=swave}
  if(method=="wdr"){
    swdr1=0;swdr2=0;swdr3=0;swdr4=0
    for(i in 1:n){
      swdr1=swdr1+qcoo[,,i]%*%qcoo[,,i]/n
      swdr2=swdr2+qcoo[,,i]/n
      swdr3=swdr3+lcoo[,i]%*%t(lcoo[,i])%*%kt/n
      swdr4=swdr4+c(t(lcoo[,i])%*%kt%*%lcoo[,i]/n)}
    swdr=swdr1+swdr2%*%swdr2+swdr3%*%swdr3+swdr4*swdr3;smethod=swdr}
  v=svd(kt.nhalf%*%(kt%*%smethod)%*%kt.nhalf)$u;fcoo=kt.nhalf%*%v
  #pred=matrix(0,n,nt)
  #for(i in 1:nt) pred[,i]=c(t(fcoo[,i])%*%kt.half%*%ax%*%kt.half%*%t(xcoo))
  pred=xcoo %*% kt.half%*%ax%*%kt.half %*% fcoo
  pred.test=x.test.coo %*% kt.half%*%ax%*%kt.half %*% fcoo
  if(px>1){
    nt = nt/px
    tmp = matrix(0,n,nt)
    for(j in 1:px){
      tmp = tmp+pred.train[,((j-1)*nt+1):(j*nt) ]
    }
    pred.train = tmp
    tmp = matrix(0,n,nt.test)
    for(j in 1:px){
      tmp = tmp+pred.test[,((j-1)*nt.test+1):(j*nt.test) ]
    }
    pred.test = tmp
  }
  return(list(pred=pred,pred.test=pred.test,et=et,ex=ex,ey=ey,kx=kx))
}
non.sdr.fd = function(x,y,tt=NULL, method='fgsir',ytype='function',
                      et=0.05, ex=0.05, ey=0.05, cvx=FALSE, ncv=10,
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
      tmp.gram = gramx(ycoo,kt,py)
      KY=tmp.gram$gram
      normY = tmp.gram$norm
    }else if(py==1){
      ycoo = t(ycoo) 
      tmp.gram = gramx(ycoo,kt,py)
      KY=tmp.gram$gram
      normY = tmp.gram$norm
    }
  }else if(!is.fd(y) & ytype=='function'){
    stop('y should be a functional object (fda package)')
    #ycoo=xcoomat(y,tte,ttt,et,kern)
    #KY=gramx(ycoo,kt)$gram
  }
  if(ytype=='categorical'){
    yclass=TRUE
  }else{ yclass=FALSE }
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
  gamma = temp$gamma
  KX = temp$gram
  normX = tmp.gram$norm
  KX.test = KX
  # Rcpp - tuning parameter search
  # if(cvx){
  #   cvwinx <- c(2,exp(seq(log(0.05),log(20),len=ncv-1)))
  #   cvwiny <- c(2,exp(seq(log(0.05),log(20),len=ncv-1)))
  #   cvex <-  c(n^(-1/4),exp(seq(log(n^(-1/4)/50),log(n^(-1/4)*50),len=ncv-1)))
  #   
  #   cvxout <- cvxsearch(cvwinx, cvwiny, cvex, normX, normY, as.numeric(yclass), KY) 
  #   #cvyout <- cvxsearch(cvwinx, cvwiny, cvex, normY, normX)
  #   #ey <- cvyout$ex; shyy <- cvyout$shy
  #   shx = cvxout$shx; shy=cvxout$shy; ex=cvxout$ex; ey = cvxout$ey;
  #   if(!yclass){
  #     gy <- gam(normY)/shy
  #     KY <- exp(-gy*normY)      
  #   }
  #   gamma <- gam(normX)/shx
  #   KX <- exp(-gamma*normX)
  # }
  # 
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
    # R version
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
    # Rcpp version
    # M <- faveloop(GX,GY,KY,KX,Q,ex,ey, as.numeric(yclass))
    A <- eigen(sym(M), sym=TRUE)
    eigval <- A$values
    A <- A$vec
    Mat <- GX.inv.half %*% A
    pred <- KX.test %*% Mat
  }
  out <- list(eval = eigval, pred = pred)
  return(out)
}

fpca = function(ftn){
  temp = ftn$coef
  n = dim(temp)[2]; p = dim(temp)[3]; nt = dim(temp[1])
  
  if(is.na(p)) p = 1
  xcoef = temp
  GB = bsplinepen(ftn$basis,0)
  
  basis =ftn$basis
  
  one = matrix(1, n, 1)
  Q = diag(n) - one %*% t(one) / n
  B.half = matpower(GB,0.5)
  
  if(p==1){
    Sigma = B.half %*% xcoef %*% Q %*% t(xcoef) %*% B.half / n
    egn = eigen(Sigma, sym=TRUE)
    B.inv.half = matpower(GB,-0.5)
    pred = Q %*% t(xcoef) %*% GB %*% B.inv.half %*% egn$vec
    out = list(pred=pred, eval=egn$val, mat =  B.inv.half %*% egn$vec)
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
    out = list(pred=pred, eval=egn$val, mat=D.inv.half %*% egn$vec)
  }
  return(out)
}



###########################################################################
#       function: fpc.loading  - principal component functions
###########################################################################
# Example
# tmp=fpc.loading(fd.sub, fpca.out,npt=20)
### 1st PC of 1st functional variable
# plot(tmp$tt, tmp$pcftn[,1,1], type='l')  
### 1st PC of 3rd functional variable
# plot(tmp$tt, tmp$pcftn[,1,3], type='l')  
fpc.loading = function(fd, fpc, npt=20){
  tmp.dm=dim(fd$coefs)
  nt=tmp.dm[1];n = tmp.dm[2];p=tmp.dm[3]
  
  
  basis=fd$basis
  nbasis = basis$nbasis
  argvals = basis$rangeval
  tseq = seq(argvals[1],argvals[2], len=npt)  # len: number of points for drawing the curve
  eval.t =eval.basis(tseq, basis) 
  if(p==1){
    out = eval.t %*% fpc$mat
  }else if(p>1){
    out = array(0,c(npt,nbasis,p))
    for(j in 1:p){
      ind = ((j-1)*nt+1) : (j*nt)
      out[,,j]=eval.t %*% fpc$mat[ind,1:nbasis]
    }
  }
  return(list(tt=tseq,pcftn=out))
}

#------------------------------------------------------------#
# FSIR (2003, Ferre and Yao), FSAVE (2014, Lian and Li)
# 11/03/2017
#------------------------------------------------------------#
# x : n X nt matrix  / y : n vector
#------------------------------------------------------------#

# FSIR

# d : number of edr directinos
fsir_fsave = function(x, y, tt, H=10, K=3, d=3,nknots=4, xtest=NULL){
  if(is.null(xtest)) xtest=x
  X = t(x); Y=y
  sortxy = sortY(X,Y)
  X = sortxy$X; Y= sortxy$Y
  ##################################
  ########  SIR.PCA ################
  ##################################
  p<-1 #p>1 only used for multiple functional study
  
  ninternal=4 # number basis is ninternal+4
  mybreaks<-c()
  for (j in 1:p){
    mybreaks<-c(mybreaks,c(seq(0,1,length=ninternal+2),1,1)+j-1)
  }
  mybreaks<-mybreaks[1:(length(mybreaks)-2)]
  mybasis<-create.bspline.basis(rangeval=c(0,p),breaks=mybreaks)
  mybasismatrix<-getbasismatrix(evalarg=tt,mybasis)
  mypen0<-eval.penalty(mybasis)
  
  myfd<-Data2fd(argvals=tt,y=X,mybasis) # number column of X is n
  myfd.test <-Data2fd(argvals=tt,y=t(xtest),mybasis) # number column of X is n 
  mypca<-pca.fd(myfd,nharm=mybasis$nbasis-1)
  #############################################
  
  myscores = mypca$score    # used in SAVE
  myEXgY<-getEXgY(H,myscores,Y) # used in SAVE
  
  Gamma<-diag(mypca$values)     # used in SAVE
  
  myVEXgY<-cov(myEXgY)
  
  
  invrootGamma<-diag(1/sqrt(mypca$values[1:K])) # used in SAVE
  temp<-eigen(invrootGamma%*%myVEXgY[1:K,1:K]%*%invrootGamma)$vectors[,1:min(d,K)]
  betascore<-invrootGamma%*%temp    # Coefficient of beta under PCA basis
  beta.fd = fd(mypca$harmonics$coef[,1:K]%*%betascore, mybasis)    
  
  pred_sir_pc = inprod(myfd, beta.fd)   # <X,Betas>
  pred_sir_pc.test = inprod(myfd.test, beta.fd)
  beta_sir_pc<-mybasismatrix%*%(mypca$harmonics$coef[,1:K]%*%betascore)  # Evaluated beta functions
  
  ##################################
  ########  SAVE.PCA ##############
  ##################################
  
  myVXgY<-getVXgY(H,myscores,Y)
  temp<-matrix(0,nrow=K,ncol=K)
  for (h in 1:H){
    temp<-temp+diag(K)-2*invrootGamma%*%myVXgY[1:K,1:K,h]%*%invrootGamma+
      invrootGamma%*%myVXgY[1:K,1:K,h]%*%diag(1/mypca$values[1:K])%*%myVXgY[1:K,1:K,h]%*%invrootGamma
  }
  betascore<-invrootGamma%*%(eigen(temp)$vectors[,1:min(d,K)])
  beta_save_pc<-mybasismatrix%*%(mypca$harmonics$coef[,1:K]%*%betascore)
  beta.fd = fd(mypca$harmonics$coef[,1:K]%*%betascore, mybasis)
  pred_save_pc = inprod(myfd, beta.fd)
  pred_save_pc.test = inprod(myfd.test, beta.fd)
  ##################################
  ########  SIR.SPLINE ##############
  ##################################
  ##################################################
  p=1
  mybreaks<-c()
  for (j in 1:p){
    mybreaks<-c(mybreaks,c(seq(0,1,length=nknots+2),1,1)+j-1)
  }
  mybreaks<-mybreaks[1:(length(mybreaks)-2)]
  mybasis<-create.bspline.basis(rangeval=c(0,p),breaks=mybreaks,norder=3)
  mybasismatrix<-getbasismatrix(tt,mybasis)
  mypen0<-eval.penalty(mybasis)
  rootmypen0<-squareroot(mypen0)
  myfd<-Data2fd(argvals=tt,y=X,mybasis)
  mycoefs<-t(myfd$coefs)
  if (nknots==0){
    mypen0<-t(mybasismatrix)%*%mybasismatrix/ntime
    rootmypen0<-squareroot(mypen0)
    mycoefs<-t(solve(t(mybasismatrix)%*%mybasismatrix)%*%t(mybasismatrix)%*%X)
  }
  #################################################
  
  myEXgY<-getEXgY(H,mycoefs,Y)   
  myVEXgY<-rootmypen0%*%cov(myEXgY)%*%rootmypen0
  Gamma<-cov(mycoefs)   # used in SAVE
  invrootGamma<-squareroot(solve(Gamma))  # used in SAVE
  
  temp<-eigen(invrootGamma%*%cov(myEXgY)%*%invrootGamma)$vectors[,1:min(d,mybasis$nbasis)]
  betacoef<-solve(mypen0)%*%invrootGamma%*%temp
  beta_sir_bsp<-mybasismatrix%*%betacoef
  
  beta.fd = fd(betacoef, mybasis)
  pred_sir_bsp = inprod(myfd, beta.fd)
  pred_sir_bsp.test = inprod(myfd.test, beta.fd)
  ##################################
  ########  SAVE.SPLINE ##############
  ##################################
  
  
  myVXgY<-getVXgY(H,mycoefs,Y)  
  
  temp<-matrix(0,nrow=mybasis$nbasis,ncol=mybasis$nbasis)
  for (h in 1:H){
    temp<-temp+diag(mybasis$nbasis)-2*invrootGamma%*%myVXgY[,,h]%*%invrootGamma+
      invrootGamma%*%myVXgY[,,h]%*%solve(Gamma)%*%myVXgY[,,h]%*%invrootGamma
  }
  betacoef<-solve(mypen0)%*%invrootGamma%*%(eigen(temp)$vectors[,1:min(d,mybasis$nbasis)])
  beta_save_bsp<-mybasismatrix%*%betacoef
  
  beta.fd = fd(betacoef, mybasis)
  pred_save_bsp = inprod(myfd, beta.fd)
  pred_save_bsp.test = inprod(myfd.test, beta.fd)
  #############################
  # Output
  #############################
  pred.all = array(0,c(4,length(Y),min(K,d)))
  pred.all[1,,]=pred_sir_pc
  pred.all[2,,]=pred_sir_bsp
  pred.all[3,,]=pred_save_pc
  pred.all[4,,]=pred_save_bsp
  pred.all=pred.all[,sortxy$recov_ind,]
  
  pred.all.test = array(0,c(4,nrow(xtest),min(K,d)))
  pred.all.test[1,,]=pred_sir_pc.test
  pred.all.test[2,,]=pred_sir_bsp.test
  pred.all.test[3,,]=pred_save_pc.test
  pred.all.test[4,,]=pred_save_bsp.test
  return(list(pred.all=pred.all, pred.all.test=pred.all.test))
}

sortY<-function(X,Y){
  index<-order(Y);recov_ind=rank(Y)
  Y<-Y[index]
  X<-X[,index]
  return(list(X=X,Y=Y,sort_ind=index,recov_ind=recov_ind))  # sort(Y)[rank(Y)] = Y
}

get.spline.coef = function(X, nknots){
  p=1
  mybreaks<-c()
  for (j in 1:p){
    mybreaks<-c(mybreaks,c(seq(0,1,length=nknots+2),1,1)+j-1)
  }
  mybreaks<-mybreaks[1:(length(mybreaks)-2)]
  mybasis<-create.bspline.basis(rangeval=c(0,p),breaks=mybreaks,norder=3)
  mybasismatrix<-getbasismatrix(seq(0,p,length=ntime*p),mybasis)
  mypen0<-eval.penalty(mybasis)
  rootmypen0<-squareroot(mypen0)
  myfd<-Data2fd(argvals=grid,y=X,mybasis)
  mycoefs<-t(myfd$coefs)
  if (nknots==0){
    mypen0<-t(mybasismatrix)%*%mybasismatrix/ntime
    rootmypen0<-squareroot(mypen0)
    mycoefs<-t(solve(t(mybasismatrix)%*%mybasismatrix)%*%t(mybasismatrix)%*%X)
  }
  return(list(mycoefs=mycoefs,mybasis=mybasis,mybasismatrix=mybasismatrix))
}
get.pca = function(X){
  
  return(list(mypca=mypca,mybasis=mybasis,mybasismatrix=mybasismatrix))
}

getEXgY<-function(H,X,Y){
  slicepos<-as.vector(round(quantile(1:n,seq(0,1,length=H+1))))
  nslice<-diff(slicepos); nslice[1]<-nslice[1]+1
  sX<-apply((X),2,cumsum)
  sX<-sX[slicepos[-1],]
  sX<-rbind(sX[1,], apply(sX,2,diff))
  sX<-sX/nslice
  return(sX) # H by dimension matrix
}

squareroot<-function(A){
  A<-(A+t(A))/2
  temp<-svd(A)
  return( temp$u%*%diag(sqrt(temp$d))%*%t(temp$u) )
}

getVXgY<-function(H,X,Y){
  slicepos<-as.vector(round(quantile(1:n,seq(0,1,length=H+1))))
  slicepos[1]<-0
  VXgY<-array(0,dim=c(dim(X)[2],dim(X)[2],H))
  for (h in 1:H){
    index<-(slicepos[h]+1):(slicepos[h+1])
    VXgY[,,h]<-cov(X[index,])
  }
  return(VXgY)
}





fsir <- function(x, y, H, k){
  n <- length(y)
  J <- ncol(x)
  z <- center(x)
  ### Make a slice  # of slice is H
  as.vector(round(quantile(1:n,seq(0,1,length=H+1))))
  
  d <- (max(y)-min(y))/H   ## length of slice
  m <- min(y)
  yindex <- rep(0,n)         # yy[i]=h implies i'th observation is in slice h out of H
  yindex <- (y<=m+d)*1
  yindex <- yindex+(y>m+(H-1)*d)*H
  for ( i in 2:(H-1))
  {
    yindex <- yindex+(m+(i-1)*d<y)*(y<=m+i*d)*i
  }
  ####  Calculate ph, mh, V , mh[h,] is a vector
  ph <- rep(0, H) #ps
  hh <- rep(0, J) # hs
  R <- t(z)%*%z/n # RnJn
  AA<-matrix(0,J,J)
  V <- matrix(0, J, J)
  for ( i in 1:H)
  {
    ph[i] <- sum(yindex==i)/n
    for (j in 1:J)
      if(ph[i]!=0) hh[j] <- mean(z[yindex==i,j])
      V <- V + ph[i]*hh%*%t(hh)
  }
  # k : reduced dimension
  
  p <- eigen(R)$vec
  Pk <- p[,1:k]%*%((solve(p%*%t(p))%*%t(p))[1:k,])
  Rk <- Pk %*% R %*% Pk
  Rk.inv.half <- matpower(Rk, -1/2)
  M <- Rk.inv.half %*% V %*% Rk.inv.half
  beta <- eigen(M,sym=TRUE)$vec[,1:k]
  xi <- x%*%beta
  return(list(beta=beta,pred=xi))
}
center <- function(x){
  return(t(t(x)-apply(x,2,mean)))}