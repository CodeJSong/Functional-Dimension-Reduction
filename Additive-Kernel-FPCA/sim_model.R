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
###########
# Model I-1
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
# Model I-2
###########
sim.model.12 = function(n=100, nt=20, p=1, balanced=TRUE){
  tt = 1:nt/nt
  y = rbinom(n,1,0.5)+1
  n1 = sum(y==0); n2 = n - n1

  x= matrix(0, n, nt)
  for(i in 1:n){
    x[i,] <-   cos(tt * (y[i])^2)+ rnorm(nt,0 , 1)
  }
  out = list(n=n, nt=nt, x=x, tt=tt, y=y,p=p)
  return(out)
}


############
# Model II-1
############
sim.model.21 = function(n=100, nt=20, p=1, balanced=TRUE){
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


  x = NULL
  x$x1 = split(t(x1), rep(1:nrow(x1), each = ncol(x1)))
  x$x2 = split(t(x2), rep(1:nrow(x2), each = ncol(x2)))

  tl = NULL
  tl$x1 = split(t(tt), rep(1:nrow(tt), each = ncol(tt)))
  tl$x2 = split(t(tt), rep(1:nrow(tt), each = ncol(tt)))

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
