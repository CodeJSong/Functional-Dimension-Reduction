source('nonlinearfdr_fd.R')
load("hand_fd.RData")
library(fda)
clx=as.character(c(2,3,6))
cly=c('a','b','c')

xtr.ind = which(sepx=='train')
cl236 = which(xcl %in% clx)
xtr.ind=intersect(xtr.ind,cl236)

ytr.ind = which(sepy=='train')
clabc = which(ycl %in% cly)
ytr.ind=intersect(ytr.ind,clabc)

xtr = fd(digit.fd.coef$coef[,xtr.ind,],digit.fd.coef$basis)
ytr = fd(letter.fd.coef$coef[,ytr.ind,],letter.fd.coef$basis)
plot(ytr[,1])

tmp=nonilnear.sdr.fd(xtr,ytr)
pred = tmp$pred
plot(pred[,1:2],col=as.numeric(xcl[xtr.ind]))


