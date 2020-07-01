source('fdr_sub.R')
source('fdr_main.R')

bike.all=read.csv('hour.csv')  
ind = (bike.all$workingday==0) 
bike.tmp = bike.all[ind,]

# remove the dates which does not have full hours
tmp=split(bike.tmp[,c("hr")],bike.tmp$dteday)
rm.list = names(which(unlist(lapply(tmp,length))!=24))
bike = bike.tmp[!(bike.tmp$dteday %in% rm.list),]


#temp: Normalized temperature in Celsius. The values are divided to 41 (max)
bike$dteday = as.Date(bike$dteday,format='%Y-%m-%d')
dates = (unique(bike$dteday))
dates = sort(as.Date(dates, format="%m/%d/%Y"))


hr = unique(bike$hr)


n = length(dates)
nt = length(hr)


dat=array(0,c(nt,n,4))
dimnames(dat) = list(hr/24,dates, c('temp','count', 'humid','windspeed'))

for(i in 1:n){
  index = which(bike$dteday==dates[i])
  dat[,i,1] = bike$temp[index]*41
  dat[,i,2] = bike$cnt[index]
  dat[,i,3] = bike$hum[index]*100
  dat[,i,4] = bike$windspeed[index]*67
}

##############################
library(fda)
nbasis=31

tt = seq(0,1,len=nt)


## Get Original Data for dimension reduction
nt = length(tt)
tt = seq(0,1,len=nt)
x = (dat[,,1])  # nt X n
y = (dat[,,2])
y = log(y+1)

## Get functional object  for linear regerssion
tmp = get.fd(x,tt,kern="bspline",nbasis=nbasis,ncv=30)
x.fd = fd(tmp$coef,tmp$basis)
tmp = get.fd(y,tt,kern="bspline",nbasis=nbasis,ncv=30)
y.fd = fd(tmp$coef,tmp$basis)


# Basic setting for dimension reduction
ncv=30

ytype='function'
kern='gauss' 
#kern='brown'
# For Linear regression
basis=tmp$basis
betabasis = tmp$basis
betafdPar = fdPar(betabasis)
betaList = list(betafdPar, betafdPar)


###################################
## Training (2011) / Test  (2012)
###################################
set.seed(1)
ndim=10
method.names = c("linear", "nonlinear", "FPCA/Linear", "FPCA/Nonlinear", "WIRE","WAVE","WDR", "WIRE","WAVE","WDR", "fGSIR/Linear","fgsir/Non")
mse.1112 = array(0,c(length(method.names), ndim)); 
rownames(mse.1112)= method.names
mse.1112.train=mse.1112
dates[100];dates[101]


train.ind = 1:100;n.train=100
test.ind = (1:n)[-train.ind];n.test=n-n.train


# get x-y data
x.tr = x[,train.ind]
y.tr = y[,train.ind]
x.test = x[,test.ind]
y.test = y[,test.ind]

x.fd.tr = fd(x.fd$coefs[,train.ind],x.fd$basis)
y.fd.tr = fd(y.fd$coefs[,train.ind],y.fd$basis)

x.fd.test = fd(x.fd$coefs[,test.ind],x.fd$basis)
y.fd.test = fd(y.fd$coefs[,test.ind],y.fd$basis)



# Linear regression w/o dim. reduction   
# y: knee.test ~ x:hip.test   directly
f.lm <- fRegress(y.fd.tr ~ x.fd.tr)
beta0.tmp=fd(f.lm$betaestlist[[1]]$fd$coefs,basis)
beta.tmp=fd(f.lm$betaestlist[[2]]$fd$coefs,basis)

y.test.hat.fd = beta0.tmp + beta.tmp*x.fd.test
y.hat = eval.fd(tt,y.test.hat.fd)

r1 = (y.hat-y.test)^2

mse.1112[ 1,]=(sum(r1)/nt/n.test)

# train.
mse.1112.train[1,] = sum((eval.fd(tt,f.lm$yhatfdobj$fd) - y.tr)^2)/nt/n.train


# Nonlinear regression w/o dim. reduction
nonlinear.out=nonlinear.est(y=y.tr, x=t(x.tr), x.test=t(x.test), ytype=ytype, xtype='function',
                            kern=kern, h=.1,cvh=T , kfolds=5,ncv=20)
nonlinear.out=nonlinear.est(y=y.tr, x=x.fd.tr, x.test=x.fd.test, ytype=ytype, xtype='function',
                            kern=kern, h=.1,cvh=T , kfolds=5,ncv=20)

r2=(nonlinear.out$y.test.hat-y.test)^2
mse.1112[2,] = sum(r2)/nt/n.test

r2=(nonlinear.out$y.hat-y.tr)^2
mse.1112.train[2,] = sum(r2)/nt/n.train
# with PCA using training
pca.train=pca.fd(x.fd.tr,nharm=ndim)



# suff. dim. reduction using training data
tmp=tune.fsdr(t(x.tr),t(y.tr), n.train, nt, kern=kern, ytype=ytype, ncv=6)

et=tmp[[1]];ex=tmp[[2]];ey=tmp[[3]]

tmp1 = fsdr(x=t(x.tr),y=t(y.tr),kern=kern, et=et,ex=ex,ey=ey,method='wire',ytype='function',x.tes=t(x.test))
tmp2 = fsdr(x=t(x.tr),y=t(y.tr),kern=kern, et=et,ex=ex,ey=ey,method='wave',ytype='function',x.tes=t(x.test))
tmp3 = fsdr(x=t(x.tr),y=t(y.tr),kern=kern, et=et,ex=ex,ey=ey,method='wdr',ytype='function',x.tes=t(x.test))
# 
# tmp1 = fsdr.fd(x=x.fd.tr,y=y.fd.tr,method='wire',ytype='function',x.tes=x.fd.test)
# tmp2 = fsdr.fd(x=x.fd.tr,y=y.fd.tr, method='wave',ytype='function',x.tes=x.fd.test)
# tmp3 = fsdr.fd(x=x.fd.tr,y=y.fd.tr, method='wdr',ytype='function',x.tes=x.fd.test)
tmp=NULL
tmp[[1]] = tmp1
tmp[[2]] = tmp2
tmp[[3]] = tmp3

for(p in 1:ndim){
  
  # Linear with PCA 
  
  pca.train.pred=pca.train$score[,1:p]
  pca.pred = inprod(x.fd.test, pca.train$harmonics[1:p,])
  
  lm.fpca <- fRegress(y.fd.tr ~ pca.train.pred)
  
  betafd=fd(sapply(lm.fpca$betaestlist,coef),basis)
  betamat = eval.fd(tt,betafd)
  
  
  xtmp = cbind(1,pca.pred)
  y.hat=betamat %*% t(xtmp)
  r1pca = (y.hat-y.test)^2
  mse.1112[3,p]=(sum(r1pca)/nt/n.test)
  
  r2=(eval.fd(tt,lm.fpca$yhatfdobj$fd)-y.tr)^2
  mse.1112.train[3,p]=(sum(r2)/nt/n.train)
  # 
  # lm.fpca <- fRegress(y.fd.test ~ pca.pred)
  # y.hat = eval.fd(tt,lm.fpca$yhatfdobj$fd)
  # r1pca = (y.hat-y.test)^2
  # mse.1112[3,p]=(sum(r1pca)/nt/n.test)
  
  # Nonlinear with PCA
  nonlinear.out=nonlinear.est(y=y.tr, x=pca.train.pred, x.test=pca.pred, ytype=ytype, xtype='scalar',
                              kern=kern, h=1,cvh=T , kfolds=5)
  r2=(nonlinear.out$y.test.hat-y.test)^2
  mse.1112[4,p] = sum(r2)/nt/n.test
  
  r2=(nonlinear.out$y.hat-y.tr)^2
  mse.1112.train[4,p] = sum(r2)/nt/n.train
  
  # Linear with fgsir
  non.pred.train = non.sdr.fd(x.fd.tr,y.fd.tr,method='fgsir', ytype='function')$pred[,1:p]
  non.pred.test = non.sdr.fd(x.fd.tr,y.fd.tr,method='fgsir', ytype='function', x.tes = x.fd.test)$pred[,1:p]
  lm.fgsir <- fRegress(y.fd.tr ~ non.pred.train)
  
  betafd=fd(sapply(lm.fgsir$betaestlist,coef),basis)
  betamat = eval.fd(tt,betafd)
  
  
  xtmp = cbind(1,non.pred.test)
  y.hat=betamat %*% t(xtmp)
  r2 = (y.hat-y.test)^2
  mse.1112[11,p] = sum(r2)/nt/n.test
  
  r2=(eval.fd(tt,lm.fgsir$yhatfdobj$fd)-y.tr)^2
  mse.1112.train[11,p]=(sum(r2)/nt/n.train)
  
  # Nonlinear with fgsir
  nonlinear.out=nonlinear.est(y=y.tr, x=non.pred.train, x.test=non.pred.test, ytype=ytype, xtype='scalar',
                              kern=kern, h=1,cvh=T , kfolds=5)
  r2=(nonlinear.out$y.test.hat-y.test)^2
  mse.1112[12,p] = sum(r2)/nt/n.test
  
  r2=(nonlinear.out$y.hat-y.tr)^2
  mse.1112.train[12,p] = sum(r2)/nt/n.train
  
  
  # Linear with SDR
  for(j in 1:3){
    train.sdr=tmp[[j]][[1]][,1:p]
    newx=tmp[[j]][[2]][,1:p]
    
    
    # Linear regression w/ dim. reduction
    sdr.regress = fRegress(y.fd.tr ~ train.sdr)
    betafd=fd(sapply(sdr.regress$betaestlist,coef),basis)
    betamat = eval.fd(tt,betafd)
    
    
    xtmp = cbind(1,newx)
    y.hat=betamat %*% t(xtmp)
    
    
    r2 = (y.hat-y.test)^2
    
    mse.1112[4+j,p]=(sum(r2)/nt/n.test)
    
    
    r2=(eval.fd(tt,sdr.regress$yhatfdobj$fd)-y.tr)^2
    mse.1112.train[4+j,p]=(sum(r2)/nt/n.train)
  }
  
  # NonLinear with SDR
  for(j in 1:3){
    train.sdr=tmp[[j]][[1]][,1:p]
    newx=tmp[[j]][[2]][,1:p]
    
    # Linear regression w/ dim. reduction
    nonlinear.out=nonlinear.est(y=y.tr, x=train.sdr, x.test=newx, ytype=ytype, xtype='scalar',
                                kern=kern, h=1,cvh=T , kfolds=5)
    r2=(nonlinear.out$y.test.hat-y.test)^2
    
    mse.1112[7+j,p]=(sum(r2)/nt/n.test)
    
    r2=(nonlinear.out$y.hat-y.tr)^2
    
    mse.1112.train[7+j,p]=(sum(r2)/nt/n.train)
  }
  
  
}
round(mse.1112.train,3)
round(mse.1112,3)


save(mse.1112,mse.1112.train,file='bikeshare_11_12.RData')
###################################
## Training / Test    randomly 100 times
###################################
nsim =100
ndim = 10
method.names = c("linear", "nonlinear", "FPCA/Linear", "FPCA/Nonlinear", "WIRE","WAVE","WDR", "WIRE","WAVE","WDR")
mse.all = array(0,c(nsim, length(method.names), ndim)); 
dimnames(mse.all)[[2]]= method.names
mse.all.train=mse.all
## PCA


# setting
n.train = 100
n.test = n-n.train

set.seed(1)
# foreach(isim = 1:nsim)%do%{
for(isim in 1:nsim){
  
  start=Sys.time()
  train.ind = sample(1:n,n.train)
  test.ind = (1:n)[-train.ind]
  
  # get x-y data
  x.tr = x[,train.ind]
  y.tr = y[,train.ind]
  x.test = x[,test.ind]
  y.test = y[,test.ind]
  
  x.fd.tr = fd(x.fd$coefs[,train.ind],x.fd$basis)
  y.fd.tr = fd(y.fd$coefs[,train.ind],y.fd$basis)
  
  x.fd.test = fd(x.fd$coefs[,test.ind],x.fd$basis)
  y.fd.test = fd(y.fd$coefs[,test.ind],y.fd$basis)
  
  
  # Linear regression w/o dim. reduction   
  # y: knee.test ~ x:hip.test   directly
  f.lm <- fRegress(y.fd.tr ~ x.fd.tr)
  beta0.tmp=fd(f.lm$betaestlist[[1]]$fd$coefs,basis)
  beta.tmp=fd(f.lm$betaestlist[[2]]$fd$coefs,basis)
  
  y.test.hat.fd = beta0.tmp + beta.tmp*x.fd.test
  y.hat = eval.fd(tt,y.test.hat.fd)
  
  r1 = (y.hat-y.test)^2
  
  mse.all[isim, 1,]=(sum(r1)/nt/n.test)
  
  # train.
  mse.all.train[isim,1,] = sum((eval.fd(tt,f.lm$yhatfdobj$fd) - y.tr)^2)/nt/n.train
  
  
  # Nonlinear regression w/o dim. reduction
  nonlinear.out=nonlinear.est(y=y.tr, x=t(x.tr), x.test=t(x.test), ytype=ytype, xtype='function',
                              kern=kern, h=.1,cvh=T , kfolds=5,ncv=20)
  r2=(nonlinear.out$y.test.hat-y.test)^2
  mse.all[isim,2,] = sum(r2)/nt/n.test
  
  r2=(nonlinear.out$y.hat-y.tr)^2
  mse.all.train[isim,2,] = sum(r2)/nt/n.train
  # with PCA using training
  pca.train=pca.fd(x.fd.tr,nharm=ndim)
  
  
  # suff. dim. reduction using training data
  tmp=tune.fsdr(t(x.tr),t(y.tr), n.train, nt, kern=kern, ytype=ytype, ncv=6)
  
  et=tmp[[1]];ex=tmp[[2]];ey=tmp[[3]]
  
  tmp1 = fsdr(x=t(x.tr),y=t(y.tr),kern=kern, et=et,ex=ex,ey=ey,method='wire',ytype='function',x.tes=t(x.test))
  tmp2 = fsdr(x=t(x.tr),y=t(y.tr),kern=kern, et=et,ex=ex,ey=ey,method='wave',ytype='function',x.tes=t(x.test))
  tmp3 = fsdr(x=t(x.tr),y=t(y.tr),kern=kern, et=et,ex=ex,ey=ey,method='wdr',ytype='function',x.tes=t(x.test))
  tmp=NULL
  tmp[[1]] = tmp1
  tmp[[2]] = tmp2
  tmp[[3]] = tmp3
  
  for(p in 1:ndim){
    
    # Linear with PCA 
    
    pca.train.pred=pca.train$score[,1:p]
    pca.pred = inprod(x.fd.test, pca.train$harmonics[1:p,])
    
    lm.fpca <- fRegress(y.fd.tr ~ pca.train.pred)
    
    betafd=fd(sapply(lm.fpca$betaestlist,coef),basis)
    betamat = eval.fd(tt,betafd)
    
    
    xtmp = cbind(1,pca.pred)
    y.hat=betamat %*% t(xtmp)
    r1pca = (y.hat-y.test)^2
    mse.all[isim,3,p]=(sum(r1pca)/nt/n.test)
    
    r2=(eval.fd(tt,lm.fpca$yhatfdobj$fd)-y.tr)^2
    mse.all.train[isim,3,p]=(sum(r2)/nt/n.train)
    # 
    # lm.fpca <- fRegress(y.fd.test ~ pca.pred)
    # y.hat = eval.fd(tt,lm.fpca$yhatfdobj$fd)
    # r1pca = (y.hat-y.test)^2
    # mse.all[isim,3,p]=(sum(r1pca)/nt/n.test)
    
    # Nonlinear with PCA
    nonlinear.out=nonlinear.est(y=y.tr, x=pca.train.pred, x.test=pca.pred, ytype=ytype, xtype='scalar',
                  kern=kern, h=1,cvh=T , kfolds=5)
    r2=(nonlinear.out$y.test.hat-y.test)^2
    mse.all[isim,4,p] = sum(r2)/nt/n.test
    
    r2=(nonlinear.out$y.hat-y.tr)^2
    mse.all.train[isim,4,p] = sum(r2)/nt/n.train
    
    
    # Linear with SDR
    for(j in 1:3){
      train.sdr=tmp[[j]][[1]][,1:p]
      newx=tmp[[j]][[2]][,1:p]
      
      
      # Linear regression w/ dim. reduction
      sdr.regress = fRegress(y.fd.tr ~ train.sdr)
      betafd=fd(sapply(sdr.regress$betaestlist,coef),basis)
      betamat = eval.fd(tt,betafd)
      
      
      xtmp = cbind(1,newx)
      y.hat=betamat %*% t(xtmp)
      
      
      r2 = (y.hat-y.test)^2
      
      mse.all[isim,4+j,p]=(sum(r2)/nt/n.test)
      
      
      r2=(eval.fd(tt,sdr.regress$yhatfdobj$fd)-y.tr)^2
      mse.all.train[isim,4+j,p]=(sum(r2)/nt/n.train)
    }
    
    # NonLinear with SDR
    for(j in 1:3){
      train.sdr=tmp[[j]][[1]][,1:p]
      newx=tmp[[j]][[2]][,1:p]
      
      # Linear regression w/ dim. reduction
      nonlinear.out=nonlinear.est(y=y.tr, x=train.sdr, x.test=newx, ytype=ytype, xtype='scalar',
                                  kern=kern, h=1,cvh=T , kfolds=5)
      r2=(nonlinear.out$y.test.hat-y.test)^2
      
      mse.all[isim,7+j,p]=(sum(r2)/nt/n.test)
      
      r2=(nonlinear.out$y.hat-y.tr)^2
      
      mse.all.train[isim,7+j,p]=(sum(r2)/nt/n.train)
    }
    
    
  }
  cat(isim,'-th sample is done with ',Sys.time()-start,'\n')
}


apply(mse.all,c(2,3),mean)
round(apply(mse.all,c(2,3),mean),3)
apply(mse.all.train,c(2,3),mean)
round(apply(mse.all.train,c(2,3),mean))

apply(mse.all,c(2,3),sd)

# do 1:100 - 2011 to predict 2012 !! upto dimeniosn 7
#boxplot(t(apply(mse.all,c(2,3),mean))[3:10,])

save(mse.all,mse.all.train,file="bike_nonlinear_gauss_01242019.RData")
#save.image('bike_test_all.RData')



###################################
## Learning the data continuously
###################################
nsim=112
ndim = 10
method.names = c("linear", "nonlinear", "FPCA/Linear", "FPCA/Nonlinear", "WIRE","WAVE","WDR", "WIRE","WAVE","WDR")
mse.all = array(0,c(nsim, length(method.names), ndim)); 
dimnames(mse.all)[[2]]= method.names
mse.all.train=mse.all
## PCA


# setting
n.train = 100
n.test = n-n.train

# foreach(isim = 1:nsim)%do%{
for(isim in 1:nsim){
  
  start=Sys.time()
  train.ind = 1:(n.train+isim-1)
  test.ind = (1:n)[-train.ind]
  
  n.train = length(train.ind);n.test=length(test.ind)
  
  # get x-y data
  x.tr = x[,train.ind]
  y.tr = y[,train.ind]
  x.test = x[,test.ind]
  y.test = y[,test.ind]
  
  x.fd.tr = fd(x.fd$coefs[,train.ind],x.fd$basis)
  y.fd.tr = fd(y.fd$coefs[,train.ind],y.fd$basis)
  
  x.fd.test = fd(x.fd$coefs[,test.ind],x.fd$basis)
  y.fd.test = fd(y.fd$coefs[,test.ind],y.fd$basis)
  
  
  # Linear regression w/o dim. reduction   
  # y: knee.test ~ x:hip.test   directly
  f.lm <- fRegress(y.fd.tr ~ x.fd.tr)
  beta0.tmp=fd(f.lm$betaestlist[[1]]$fd$coefs,basis)
  beta.tmp=fd(f.lm$betaestlist[[2]]$fd$coefs,basis)
  
  y.test.hat.fd = beta0.tmp + beta.tmp*x.fd.test
  y.hat = eval.fd(tt,y.test.hat.fd)
  
  r1 = (y.hat-y.test)^2
  
  mse.all[isim, 1,]=(sum(r1)/nt/n.test)
  
  # train.
  mse.all.train[isim,1,] = sum((eval.fd(tt,f.lm$yhatfdobj$fd) - y.tr)^2)/nt/n.train
  
  
  # Nonlinear regression w/o dim. reduction
  nonlinear.out=nonlinear.est(y=y.tr, x=t(x.tr), x.test=t(x.test), ytype=ytype, xtype='function',
                              kern=kern, h=.1,cvh=T , kfolds=5,ncv=20)
  r2=(nonlinear.out$y.test.hat-y.test)^2
  mse.all[isim,2,] = sum(r2)/nt/n.test
  
  r2=(nonlinear.out$y.hat-y.tr)^2
  mse.all.train[isim,2,] = sum(r2)/nt/n.train
  # with PCA using training
  pca.train=pca.fd(x.fd.tr,nharm=ndim)
  
  
  # suff. dim. reduction using training data
  tmp=tune.fsdr(t(x.tr),t(y.tr), n.train, nt, kern=kern, ytype=ytype, ncv=6)
  
  et=tmp[[1]];ex=tmp[[2]];ey=tmp[[3]]
  
  tmp1 = fsdr(x=t(x.tr),y=t(y.tr),kern=kern, et=et,ex=ex,ey=ey,method='wire',ytype='function',x.tes=t(x.test))
  tmp2 = fsdr(x=t(x.tr),y=t(y.tr),kern=kern, et=et,ex=ex,ey=ey,method='wave',ytype='function',x.tes=t(x.test))
  tmp3 = fsdr(x=t(x.tr),y=t(y.tr),kern=kern, et=et,ex=ex,ey=ey,method='wdr',ytype='function',x.tes=t(x.test))
  tmp=NULL
  tmp[[1]] = tmp1
  tmp[[2]] = tmp2
  tmp[[3]] = tmp3
  
  for(p in 1:ndim){
    
    # Linear with PCA 
    
    pca.train.pred=pca.train$score[,1:p]
    pca.pred = inprod(x.fd.test, pca.train$harmonics[1:p,])
    
    lm.fpca <- fRegress(y.fd.tr ~ pca.train.pred)
    
    betafd=fd(sapply(lm.fpca$betaestlist,coef),basis)
    betamat = eval.fd(tt,betafd)
    
    
    xtmp = cbind(1,pca.pred)
    y.hat=betamat %*% t(xtmp)
    r1pca = (y.hat-y.test)^2
    mse.all[isim,3,p]=(sum(r1pca)/nt/n.test)
    
    r2=(eval.fd(tt,lm.fpca$yhatfdobj$fd)-y.tr)^2
    mse.all.train[isim,3,p]=(sum(r2)/nt/n.train)
    # 
    # lm.fpca <- fRegress(y.fd.test ~ pca.pred)
    # y.hat = eval.fd(tt,lm.fpca$yhatfdobj$fd)
    # r1pca = (y.hat-y.test)^2
    # mse.all[isim,3,p]=(sum(r1pca)/nt/n.test)
    
    # Nonlinear with PCA
    nonlinear.out=nonlinear.est(y=y.tr, x=pca.train.pred, x.test=pca.pred, ytype=ytype, xtype='scalar',
                                kern=kern, h=1,cvh=T , kfolds=5)
    r2=(nonlinear.out$y.test.hat-y.test)^2
    mse.all[isim,4,p] = sum(r2)/nt/n.test
    
    r2=(nonlinear.out$y.hat-y.tr)^2
    mse.all.train[isim,4,p] = sum(r2)/nt/n.train
    
    
    # Linear with SDR
    for(j in 1:3){
      train.sdr=tmp[[j]][[1]][,1:p]
      newx=tmp[[j]][[2]][,1:p]
      
      
      # Linear regression w/ dim. reduction
      sdr.regress = fRegress(y.fd.tr ~ train.sdr)
      betafd=fd(sapply(sdr.regress$betaestlist,coef),basis)
      betamat = eval.fd(tt,betafd)
      
      
      xtmp = cbind(1,newx)
      y.hat=betamat %*% t(xtmp)
      
      
      r2 = (y.hat-y.test)^2
      
      mse.all[isim,4+j,p]=(sum(r2)/nt/n.test)
      
      
      r2=(eval.fd(tt,sdr.regress$yhatfdobj$fd)-y.tr)^2
      mse.all.train[isim,4+j,p]=(sum(r2)/nt/n.train)
    }
    
    # NonLinear with SDR
    for(j in 1:3){
      train.sdr=tmp[[j]][[1]][,1:p]
      newx=tmp[[j]][[2]][,1:p]
      
      # Linear regression w/ dim. reduction
      nonlinear.out=nonlinear.est(y=y.tr, x=train.sdr, x.test=newx, ytype=ytype, xtype='scalar',
                                  kern=kern, h=1,cvh=T , kfolds=5)
      r2=(nonlinear.out$y.test.hat-y.test)^2
      
      mse.all[isim,7+j,p]=(sum(r2)/nt/n.test)
      
      r2=(nonlinear.out$y.hat-y.tr)^2
      
      mse.all.train[isim,7+j,p]=(sum(r2)/nt/n.train)
    }
    
    
  }
  cat(isim,'-th sample is done with ',Sys.time()-start,'\n')
}


save(mse.all,mse.all.train,file="bike_nonlinear_gauss_realtime_learning.RData")
#save.image('bike_test_all.RData')

