########################################
# Models for Simulations
########################################

para.curve = function(tt, model=1){
  n = length(tt)

  if(model==1){
    out = matrix(c(cos(2*pi*tt),
                   sin(2*pi*tt)), n, 2)
  }else if(model==2){
    out = matrix(c(2*cos(2*pi*tt),
                   sin(4*pi*tt)), n, 2)
  }else if(model==3){
    out = matrix(c(sin(4*pi*tt),2*cos(2*pi*tt)), n, 2)
  }
  return(out)
}
############################################################################### 
#                     Generate beta and vt
#         beta is the coefficient of index in x; 8 beta are 
#    generated and are used in different models (not all are used)
#           vt is the coefficient of index in y; 7 vt are 
#      generated and are used in different models (not all are used
############################################################################### 
betavt=function(t,j){
  if(j==1) return(list(vt=sin(3*pi*t/2),beta=sin(1*pi*t/2)))
  if(j==2) return(list(vt=sin(5*pi*t/2),beta=sin(2*pi*t/2)))
  if(j==3) return(list(vt=(2*t-1)^3+1  ,beta=sin(3*pi*t/2)))
  if(j==4) return(list(vt=(2*t-1)^2-1  ,beta=sin(3*pi*t/2)+sin(1*pi*t/2)))
  if(j==5) return(list(vt=sin(3*pi*t/2),beta=cos(pi*(2*t-1))))
  if(j==6) return(list(vt=sin(pi*t/2)  ,beta=(2*t-1)^3+1))
  if(j==7) return(list(vt=cos(t)       ,beta=(2*t-1)^2-1))
  if(j==8) return(list(beta=sin(5*pi*t/2)))}
#################################################################################
#              generate standard Brownian motion
#       x : n X nt matrix x[1,] : 1st curve; tt : 0,1 interval
#################################################################################
generate.BM <- function(n=100, nt=20){
  nt.big=(nt-1)*10+1
  rand_norm = matrix(rnorm(n*(nt.big-1),0,sqrt(1/(nt.big-1))),nt.big-1,n)
  x.big=cbind(0,t(apply(rand_norm,2,cumsum)))
  tt=seq(from=0,to=1,length=nt);tt.big=seq(from=0,to=1,length=nt.big)
  ips=matrix(0,n,8);for(j in 1:8) for(i in 1:n) {
    ips[i,j]=sum(betavt(tt.big,j)$beta*x.big[i,])/nt.big}
  x=numeric();for(i in 1:nt) x=cbind(x,x.big[,(i-1)*10+1])
  tt=matrix(tt,n,nt,byrow=T)
  return(list(x=x,ips=ips))
}
#--------------------------------------------------------#
# Nolinear Additive Functional PCA
#--------------------------------------------------------#
###########
# Model I-1  -- nonlinear relation - 1-d functional data
###########
sim.model.11 = function(n=100, nt=20, p=1, balanced=TRUE){
  tt = 1:nt/nt
  y = rbinom(n,1,0.5)
  n1 = sum(y==0); n2 = n - n1
  rtemp = rnorm(n1, 1, .2)
  thtemp = runif(n1, 0, 360)
  xx = matrix(0, n, 2);


  xx[y==0,1] = rtemp * cos(thtemp)
  xx[y==0,2] = rtemp * sin(thtemp)

  rtemp <- rnorm(n2, 4, .5)
  thtemp <- runif(n2, 0, 360)

  xx[y==1,1] = rtemp * cos(thtemp)
  xx[y==1,2] = rtemp * sin(thtemp)

  x= matrix(0, n, nt)
  for(i in 1:n){
    x[i,] <-  xx[i,1] * (cos(tt*pi)) +  xx[i,2] * sin(tt)+ rnorm(nt,0 , .1)
  }
  out = list(n=n, nt=nt, x=x, tt=tt, y=y,p=p)
  return(out)
}

###########
# Model I-2   --  true function: linear -- just sum of a few eigenfunctions for each group
###########
sim.model.12 = function(n=100, nt=20, p=1, balanced=TRUE){
  require(fda)
  nbasis=6
  basis=create.bspline.basis(c(0,1),nbasis=nbasis)
  
  z=matrix(rnorm(n*2,0,2),n,2)
  tt = 1:nt/nt
  
  
  y = sample(2,n, replace=TRUE)
  n1 = sum(y==0); n2 = n - n1
  
  x= matrix(0, n, nt)
  
  
  y1 = z[,1]%*%t(eval.fd(tt, fd(c(1,1,rep(0,nbasis-2))))) + z[,2] %*% t(eval.fd(tt, fd(c(0,1,rep(0,nbasis-2)),basis))) + rnorm(n*nt, 0, .1)
  y2 = z[,1]%*%t(eval.fd(tt, fd(c(rep(0,nbasis-4),0,0,1,1)))) + z[,2] %*% t(eval.fd(tt, fd(c(rep(0,nbasis-5),0,0,0,1,1),basis))) + rnorm(n*nt, 0, .1)
  # y1 = z[,1]%*%t(cos(2*pi*tt)) + rnorm(n*nt, 0, .1)
  # y2 = z[,2]%*%t(sin(2*pi*tt)) + rnorm(n*nt, 0, .1)
  
  x[y==1,] = y1[y==1,]
  x[y==2,] = y2[y==2,]
  out = list(n=n, nt=nt, x=x, tt=tt, y=y,p=p)
  return(out)
  
}

############
# Model II-1
############
sim.model.21 = function(n=100, nt=20, p=2, balanced=TRUE){
  y = rbinom(n,1,0.5)
  if(balanced){
    nnt = nt
  }else{
    nnt = nt * 10
  }
  ttt = 1:nnt/nnt
  tt = matrix(0, n, nt)
  n1 = sum(y==0); n2 = n - n1
  for(i in 1:n){
    tt[i,]=sort(sample(ttt,nt))
  }

  tt1 = tt[y==0,] #+ rnorm(nt*n1,sd=0.1)
  tt2 = tt[y==1,]#+ rnorm(nt*n2,sd=0.1)

  rtemp <- rnorm(nt * n1, 1, .2)
  x11<-matrix(0,n1, nt )
  x12<-matrix(0,n1, nt )
  ftmp = para.curve(as.vector(t(tt1)), model=1)
  x11[,1:nt] <- matrix(rtemp * ftmp[,1], n1, nt, byrow=T)
  x12[,1:nt] <- matrix(rtemp * ftmp[,2], n1, nt, byrow=T)

  rtemp <- rnorm(nt * n2, 4, .4)
  x21<-matrix(0,n2, nt)
  x22<-matrix(0,n2, nt)
  ftmp = para.curve(as.vector(t(tt2)), model=1)
  x21[,1:nt] <- matrix(rtemp * ftmp[,1], n2, nt, byrow=T)
  x22[,1:nt] <- matrix(rtemp * ftmp[,2], n2, nt, byrow=T)

  x1 = matrix(0, n, nt);x2 = matrix(0, n, nt)
  x1[y==0,] = x11
  x1[y==1,] = x21
  x2[y==0,] = x12
  x2[y==1,] = x22

  if(balanced){
    x = array(0,c(n,nt,p))
    x[,,1] = x1
    x[,,2] = x2
    tl = ttt
  }else if(!balanced){
    x = NULL
    x$x1 = split(t(x1), rep(1:nrow(x1), each = ncol(x1)))
    x$x2 = split(t(x2), rep(1:nrow(x2), each = ncol(x2)))
    
    tl = NULL
    tl$x1 = split(t(tt), rep(1:nrow(tt), each = ncol(tt)))
    tl$x2 = split(t(tt), rep(1:nrow(tt), each = ncol(tt)))
  }
  out = list(n=n, nt=nt, x=x, tt=tl, y=y,p=p)
  return(out)
}

############
# Model II-2 (Customization)
############
sim.model.22 = function(n=100, nt=20, p=2, balanced=TRUE, m1=1, m2=2, r1=1, r2=2, sd1=.2, sd2=.4){
  y = rbinom(n,1,0.5)
  if(balanced){
    nnt = nt
  }else{
    nnt = nt * 10
  }
  ttt = 1:nnt/nnt
  tt = matrix(0, n, nt)
  n1 = sum(y==0); n2 = n - n1
  for(i in 1:n){
    tt[i,] = sort(sample(ttt,nt))
  }
  tt1 = tt[y==0,]
  tt2 = tt[y==1,]

  rtemp <- rnorm(nt * n1, r1, sd1)
  x11<-matrix(0,n1, nt )
  x12<-matrix(0,n1, nt )
  ftmp = para.curve(as.vector(t(tt1)), model=m1)
  x11[,1:nt] <- matrix(rtemp * ftmp[,1], n1, nt, byrow=T)
  x12[,1:nt] <- matrix(rtemp * ftmp[,2], n1, nt, byrow=T)

  rtemp <- rnorm(nt * n2, r2, sd2)
  x21<-matrix(0,n2, nt)
  x22<-matrix(0,n2, nt)
  ftmp = para.curve(as.vector(t(tt2)), model=m2)
  x21[,1:nt] <- matrix(rtemp * ftmp[,1], n2, nt, byrow=T)
  x22[,1:nt] <- matrix(rtemp * ftmp[,2], n2, nt, byrow=T)

  x1 = matrix(0, n, nt);x2 = matrix(0, n, nt)
  x1[y==0,] = x11
  x1[y==1,] = x21
  x2[y==0,] = x12
  x2[y==1,] = x22


  x = NULL
  x$x1 = split(t(x1), rep(1:nrow(x1), each = ncol(x1)))
  x$x2 = split(t(x2), rep(1:nrow(x2), each = ncol(x2)))

  tl = NULL
  tl$x1 = split(t(tt), rep(1:nrow(tt), each = ncol(tt)))
  tl$x2 = split(t(tt), rep(1:nrow(tt), each = ncol(tt)))

  out = list(n=n, nt=nt, x=x, tt=tl, y=y,p=p)
  return(out)
}

###########
# Model II-3   --  combination of linear & nonlinear 
###########
sim.model.23 = function(n=100, nt=20, p=2, balanced=TRUE){
  require(fda)
  
  nbasis=6
  basis=create.bspline.basis(c(0,1),nbasis=nbasis)
  
  z=matrix(rnorm(n*2,0,4),n,2)
  tt = 1:nt/nt
  
  
  y = sample(2,n, replace=TRUE)
  n1 = sum(y==1); n2 = sum(y==2); n3=sum(y==3); n4=sum(y==4)
  
  x= matrix(0, n, nt)
  
  
  y1 = z[,1]%*%t(eval.fd(tt, fd(c(1,1,rep(0,nbasis-2))))) + z[,2] %*% t(eval.fd(tt, fd(c(0,1,rep(0,nbasis-2)),basis))) + rnorm(n*nt, 0, .1)
  y2 = z[,1]%*%t(eval.fd(tt, fd(c(rep(0,nbasis-4),0,0,1,1)))) + z[,2] %*% t(eval.fd(tt, fd(c(rep(0,nbasis-5),0,0,0,1,1),basis))) + rnorm(n*nt, 0, .1)
  # y1 = z[,1]%*%t(cos(2*pi*tt)) + rnorm(n*nt, 0, .1)
  # y2 = z[,2]%*%t(sin(2*pi*tt)) + rnorm(n*nt, 0, .1)
  
  x[y==1,] = y1[y==1,]
  x[y==2,] = y2[y==2,]
  
  x.arr = array(0,c(n,nt,2))
  x.arr[,,1] = x
  ################################
  rtemp = rnorm(n1, 1, .2)
  thtemp = runif(n1, 0, 360)
  xx = matrix(0, n, 2);
  
  
  xx[y==1,1] = rtemp * cos(thtemp)
  xx[y==1,2] = rtemp * sin(thtemp)
  
  rtemp <- rnorm(n2, 4, .5)
  thtemp <- runif(n2, 0, 360)
  
  xx[y==2,1] = rtemp * cos(thtemp)
  xx[y==2,2] = rtemp * sin(thtemp)
  
  x= matrix(0, n, nt)
  for(i in 1:n){
    x[i,] <-  xx[i,1] * (cos(tt*pi)) +  xx[i,2] * sin(tt)+ rnorm(nt,0 , .1)
  }
  
  x.arr[,,2] = x
  out = list(n=n, nt=nt, x=x.arr, tt=tt, y=y,p=p)
  return(out)
}

############
# Model III-1
############
sim.model.31 = function(n=100, nt=20){
  
  y = sample(1:3,n, replace=TRUE)
  x = matrix(0, n, nt)
  tt = 1:nt/nt*21
  
  u = matrix(runif(2*n),n,2)
  
  # h1, h2, h3
  h1 = function(ttt){
    out = 6-abs(ttt-11)
    out[out<0] = 0
    return(out)
  }
  h2 = function(ttt){
    return(h1(ttt-4))
  }
  h3 = function(ttt){
    return(h1(ttt+4))
  }
  for(i in 1:n){
    if(y[i]==1) x[i,] = u[i,1] * h1(tt) + (1-u[i,2]) * h1(tt) + rnorm(1, sd=.1)
    if(y[i]==2) x[i,] = u[i,1] * h2(tt) + (1-u[i,2]) * h2(tt)+ rnorm(1, sd=.1)
    if(y[i]==3) x[i,] = u[i,1] * h3(tt) + (1-u[i,2]) * h3(tt)+ rnorm(1, sd=.1)
  }
  return(list(x=x,y=y, tt=tt, u=u))
}


############
# Model III-2
############
sim.model.32 = function(n=100, nt=20, p=1, balanced=TRUE){
  require(fda)
  
  nbasis=6
  basis=create.bspline.basis(c(0,1),nbasis=nbasis)
  
  z=matrix(rnorm(n*2,0,2),n,2)
  tt = 1:nt/nt
  
  
  y = sample(4,n, replace=TRUE)
  n1 = sum(y==1); n2 = sum(y==2); n3=sum(y==3); n4=sum(y==4)
  
  x= matrix(0, n, nt)
  
  
  y1 = z[,1]%*%t(eval.fd(tt, fd(c(1,1,rep(0,nbasis-2))))) + z[,2] %*% t(eval.fd(tt, fd(c(0,1,rep(0,nbasis-2)),basis))) + rnorm(n*nt, 0, .1)
  y2 = z[,1]%*%t(eval.fd(tt, fd(c(rep(0,nbasis-4),0,0,1,1)))) + z[,2] %*% t(eval.fd(tt, fd(c(rep(0,nbasis-5),0,0,0,1,1),basis))) + rnorm(n*nt, 0, .1)
  # y1 = z[,1]%*%t(cos(2*pi*tt)) + rnorm(n*nt, 0, .1)
  # y2 = z[,2]%*%t(sin(2*pi*tt)) + rnorm(n*nt, 0, .1)
  
  x[y==1,] = y1[y==1,] 
  x[y==2,] = y2[y==2,]
  x[y==3,] = y1[y==3,] * y2[y==3,]^2
  
  ################################
  rtemp = rnorm(n3, 1, .2)
  thtemp = runif(n3, 0, 360)
  xx = matrix(0, n, 2);
  
  
  xx[y==3,1] = rtemp * cos(thtemp)
  xx[y==3,2] = rtemp * sin(thtemp)
  
  rtemp <- rnorm(n4, 4, .5)
  thtemp <- runif(n4, 0, 360)
  
  xx[y==4,1] = rtemp * cos(thtemp)
  xx[y==4,2] = rtemp * sin(thtemp)
  
  
  for(i in c(which(y==4))){
    x[i,] <-  xx[i,1] * (cos(tt*pi)) +  xx[i,2] * sin(tt)+ rnorm(nt,0 , .1)
  }
  
  out = list(n=n, nt=nt, x=x, tt=tt, y=y,p=p)
  return(out)
}
############
# Model III-3 
############
sim.model.33 = function(n=100, nt=20, p=1, balanced=TRUE){
  require(fda)
  
  nbasis=6
  basis=create.bspline.basis(c(0,1),nbasis=nbasis)
  
  z=matrix(rnorm(n*2,0,2),n,2)
  tt = 1:nt/nt
  
  
  y = sample(3,n, replace=TRUE)
  n1 = sum(y==1); n2 = sum(y==2); n3=sum(y==3); n4=sum(y==4)
  
  x= matrix(0, n, nt)
  
  
  y1 = z[,1]%*%t(eval.fd(tt, fd(c(1,1,rep(0,nbasis-2))))) + z[,2] %*% t(eval.fd(tt, fd(c(0,1,rep(0,nbasis-2)),basis))) + rnorm(n*nt, 0, .2)
  y2 = z[,1]%*%t(eval.fd(tt, fd(c(rep(0,nbasis-4),0,0,1,1)))) + z[,2] %*% t(eval.fd(tt, fd(c(rep(0,nbasis-5),0,0,0,1,1),basis))) + rnorm(n*nt, 0, .2)
  # y1 = z[,1]%*%t(cos(2*pi*tt)) + rnorm(n*nt, 0, .1)
  # y2 = z[,2]%*%t(sin(2*pi*tt)) + rnorm(n*nt, 0, .1)
  
  x[y==1,] = y1[y==1,]*(y2[y==1,])
  
  ################################
  rtemp = rnorm(n2, 1, .2)
  thtemp = runif(n2, 0, 360)
  xx = matrix(0, n, 2);
  
  
  xx[y==2,1] = rtemp * cos(thtemp)
  xx[y==2,2] = rtemp * sin(thtemp)
  
  rtemp <- rnorm(n3, 2, .5)
  thtemp <- runif(n3, 0, 360)
  
  xx[y==3,1] = rtemp * cos(thtemp)
  xx[y==3,2] = rtemp * sin(thtemp)
  
  
  for(i in c(which(y==2),which(y==3))){
    x[i,] <-  xx[i,1] * (cos(tt*pi)) +  xx[i,2] * sin(tt)+ rnorm(nt,0 , .1)
  }
  
  out = list(n=n, nt=nt, x=x, tt=tt, y=y,p=p)
  return(out)
}



############
# Model III-4 nonadditive
############
sim.model.34 = function(n=100, nt=20, p=1, balanced=TRUE){
  require(fda)
  
  nbasis=6
  basis=create.bspline.basis(c(0,1),nbasis=nbasis)
  
  z=matrix(rnorm(n*2,0,2),n,2)
  tt = 1:nt/nt
  
  
  y = sample(3,n, replace=TRUE)
  n1 = sum(y==1); n2 = sum(y==2); n3=sum(y==3); n4=sum(y==4)
  
  x= matrix(0, n, nt)
  
  
  y1 = z[,1]%*%t(eval.fd(tt, fd(c(1,1,rep(0,nbasis-2))))) + z[,2] %*% t(eval.fd(tt, fd(c(0,1,rep(0,nbasis-2)),basis))) + rnorm(n*nt, 0, .2)
  y2 = z[,1]%*%t(eval.fd(tt, fd(c(rep(0,nbasis-4),0,0,1,1)))) + z[,2] %*% t(eval.fd(tt, fd(c(rep(0,nbasis-5),0,0,0,1,1),basis))) + rnorm(n*nt, 0, .2)
  

  
  ################################
  rtemp = rnorm(n, 1, .2)
  thtemp = runif(n, 0, 360)
  xx = matrix(0, n, 2);
  
  xx[,1] = rtemp * cos(thtemp)
  xx[,2] = rtemp * sin(thtemp)
  
  y3 = xx[,1] %*% t(cos(tt*pi)) + xx[,2] %*% t(sin(tt))+ matrix(rnorm(nt*n,0 , .1),n,nt)
  
  rtemp <- rnorm(n, 2, .5)
  thtemp <- runif(n, 0, 360)
  
  xx[,1] = rtemp * cos(thtemp)
  xx[,2] = rtemp * sin(thtemp)
  
  y4 =xx[,1] %*% t(cos(tt*pi)) + xx[,2] %*% t(sin(tt))+ matrix(rnorm(nt*n,0 , .1),n,nt)
  
  
  # change 
  x.non.add = matrix(0, n, nt)
  x[y==1,] = y3[y==1,]/(1+abs(y1[y==1]*y2[y==1]))
  x[y==2,] = sqrt(abs(y2[y==2,]))*y1[y==2,]^2
  x[y==3,] = y1[y==3,]*y2[y==3,]^2/(1+abs(y3[y==3,]))
  
  
  out = list(n=n, nt=nt, x=x, tt=tt, y=y,p=p)
  return(out)
}

############
# Model III-5  nonadditive
############
#x1 x2 brownian mition x1x3
sim.model.35 = function(n=100, nt=20, p=3, balanced=TRUE){
  require(fda)
  tt = 1:nt/nt
  
  nbasis=20
  basis=create.bspline.basis(c(0,1),nbasis=nbasis)
  
  
  
  x.arr = array(0, c(n, nt, p))
  #result=generate.BM(n,nt);x=result$x;ips=result$ips
  coef = matrix(rnorm(n*nbasis,0,1), nbasis,n)
  fd.x = fd(coef, basis)
  x = as.matrix(t(eval.fd(tt,fd.x))) + matrix(2*sin(pi*tt),n,nt,byrow=T)
  x.arr[,,1] = x
  
  #result=generate.BM(n,nt);x=result$x;ips=result$ips
  coef = matrix(rnorm(n*nbasis,0,1), nbasis,n)
  fd.x = fd(coef, basis)
  x = as.matrix(t(eval.fd(tt,fd.x))) + matrix(4*sin(pi*tt),n,nt,byrow=T)
  x.arr[,,2] = x
  
  
  
  
  
  x.arr[,,3] = (x.arr[,,1]*x.arr[,,2]) + matrix(rnorm(n*nt,0, .1), n,nt)
  
  x3fd.coef = get.fd(t(x.arr[,,3]), tt=tt, basis=basis)
  x3fd = fd(x3fd.coef$coef, basis)
  beta.fd = fd(c(0,0,rep(1,16),0,0),basis)
  y = inprod(beta.fd,x3fd)
  yy = rep(2,n)
  yy[y<4] = 1
  # class: 1,2,3,4
  # yy = rep(4,n)
  # ttmp = quantile(y)
  # yy[y<ttmp[2]] = 1
  # yy[(y>=ttmp[2]) &  (y<ttmp[3])] = 2
  # yy[(y>=ttmp[3]) &  (y<ttmp[4])] = 3
  # 
  # 
  # y = sample(2,n, replace=TRUE)
  # n1 = sum(y==1); n2=n-n1
  # 
  # x.arr[y==1,,3] = (matrix(rnorm(n1*nt,0,.5),n1,nt)*x.arr[y==1,,1]*x.arr[y==1,,2]) + matrix(rnorm(n1*nt,0, .1), n1,nt)
  # x.arr[y==2,,3] = (matrix(rnorm(n2*nt,0.5,.5),n2,nt)*x.arr[y==2,,1]*x.arr[y==2,,2]) + matrix(rnorm(n2*nt,0, .1), n2,nt)
  # 
  # 
  out = list(n=n, nt=nt, x=x.arr, tt=tt, y=yy,p=p)
  return(out)
}

############
# Model III-6  nonadditive  - truely
############
#x1 x2 brownian mition x1x3
sim.model.36 = function(n=100, nt=20, p=3, balanced=TRUE){
  require(fda)
  tt = 1:nt/nt
  
  nbasis=6
  basis=create.bspline.basis(c(0,1),nbasis=nbasis)
  
  
  
  x.arr = array(0, c(n, nt, p))
  #result=generate.BM(n,nt);x=result$x;ips=result$ips
  coef = matrix(rnorm(n*nbasis,0,2), nbasis,n)
  fd.x = fd(coef, basis)
  x = as.matrix(t(eval.fd(tt,fd.x))) +rnorm(n,0,1) %*% t(tt^3/2+tt^2)
  
  x.arr[,,1] = x
  
  #result=generate.BM(n,nt);x=result$x;ips=result$ips
  
  coef = matrix(rnorm(n*nbasis,0,1), nbasis,n)
  fd.x = fd(coef, basis)
  x = as.matrix(t(eval.fd(tt,fd.x))) + rnorm(n,0,1) %*% t(2*sqrt(tt)+(1-tt)^2)
  x.arr[,,2] = x
  
  
  
  
  
  x.arr[,,3] = (x.arr[,,1]*x.arr[,,2]) + matrix(rnorm(n*nt,0, .1), n,nt)
  
  x3fd.coef = get.fd(t(x.arr[,,3]), tt=tt, basis=basis)
  x3fd = fd(x3fd.coef$coef, basis)
  beta.fd = fd(c(0,0,1,1,0,0),basis)
  y = inprod(beta.fd,x3fd)
  y = as.vector(y)
  yy = rep(2,n)
  yy[y<median(y)] = 1
  p=2
  # class: 1,2,3,4
  # yy = rep(4,n)
  # ttmp = quantile(y)
  # yy[y<ttmp[2]] = 1
  # yy[(y>=ttmp[2]) &  (y<ttmp[3])] = 2
  # yy[(y>=ttmp[3]) &  (y<ttmp[4])] = 3
  # 
  # 
  # y = sample(2,n, replace=TRUE)
  # n1 = sum(y==1); n2=n-n1
  # 
  # x.arr[y==1,,3] = (matrix(rnorm(n1*nt,0,.5),n1,nt)*x.arr[y==1,,1]*x.arr[y==1,,2]) + matrix(rnorm(n1*nt,0, .1), n1,nt)
  # x.arr[y==2,,3] = (matrix(rnorm(n2*nt,0.5,.5),n2,nt)*x.arr[y==2,,1]*x.arr[y==2,,2]) + matrix(rnorm(n2*nt,0, .1), n2,nt)
  # 
  # 
  out = list(n=n, nt=nt, x=x.arr[,,1:2], tt=tt, y.org=y,y=yy,p=p)
  return(out)
}
