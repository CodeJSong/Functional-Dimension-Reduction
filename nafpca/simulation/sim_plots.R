source('simulation/sim_model.R')
#install.packages(c("MFPCA","funData"))
#install.packages("funData")
source('nafpca.R')
source('../../classification.R')

library(MFPCA);library(funData)

library(dplyr)
library(ggplot2)

scaling <- function(x,k) {
  #x = abs(x)
  out = (x - min(x)) / (max(x) - min(x)) / k
  return(out)
}
##################################################################
# Model I-1 : 1-d functional data - nonlinear relation in clustering
##################################################################
n=200; nt=20
set.seed(1)
dat = sim.model.11(n=n, nt=nt)
x = dat$x; tt = dat$tt; y = dat$y

# Original Data
plot(tt,x[1,], type='l', ylim=c(-5,5), xlab="t", ylab="X(t)", main = "Model I", col=y+1)
for(i in 1:n)lines(tt,x[i,], col=(y[i]+1))

# function estimation 
# 1. fda-based
# bspline: most natural
coef.bspline = get.fd(t(x),tt ,basisname='bspline')
fd.bsp = fd(coef.bspline$coef, coef.bspline$basis)
plot(fd.bsp, col=y+1)
evals = t(eval.fd(tt,fd.bsp))
sum((evals - x)^2) / n /nt

# fourier: tend to get x[0] = 0 , x[1] =0
coef.fourier = get.fd(t(x),tt ,basisname='fourier', nbasis=53)
fd.fourier=fd(coef.fourier$coef, coef.fourier$basis)
plot(fd.fourier, col=y+1)
evals = t(eval.fd(tt,fd.fourier))
sum((evals - x)^2) / n /nt

# 2. rkhs-based
# gaussian: natural
coef.gauss = get.fd.rkhs(x, tt, kern="gauss")
evals = coef.gauss$coef %*% coef.gauss$kt
plot(tt, evals[1,], col=y[1]+1, type='l', ylim=c(min(evals),max(evals)))
for(j in 1:n){
  lines(tt, evals[j,], col=y[j]+1)
}
sum((evals - x)^2) / n /nt
# brownian: natural  --- the best in this situation
coef.brown = get.fd.rkhs(x, tt, kern="brown")
evals = coef.brown$coef %*% coef.brown$kt
plot(tt, evals[1,], col=y[1]+1, type='l', ylim=c(min(evals),max(evals)))
for(j in 1:n){
  lines(tt, evals[j,], col=y[j]+1)
}
sum((evals - x)^2) / n /nt

# Nonlinear FPCA
tmp=nafpca(x,tt, shx=11, nbasis=21,gamma.tune = TRUE,
           basisname='bspline')
tmp=nafpca(x,tt, shx=11, nbasis=30,gamma.tune = TRUE,
           basisname='fourier')
tmp=nafpca(x,tt, shx=11, nbasis=30,gamma.tune = TRUE,ncv1=30,ncv2=30,
           basisname='gauss')
tmp=nafpca(x,tt, shx=11, nbasis=30,gamma.tune = TRUE,ncv1=30,ncv2=20,
           basisname='brown')
tmp=nafpca(x,tt, shx=11, nbasis=30,gamma.tune = TRUE,ncv1=30,ncv2=20,
           basisname='bspline', kernel='poly',
           c=.1,d=2)
tmp=nafpca(x,tt, shx=11, nbasis=30,gamma.tune = TRUE,
           basisname='fourier', kernel='poly',
           c=.5,d=2)
tmp=nafpca(x,tt, shx=11, nbasis=30,gamma.tune = TRUE,
           basisname='gauss', kernel='poly',
           c=.5,d=2)
tmp=nafpca(x,tt, shx=11, nbasis=30,gamma.tune = TRUE,
           basisname='brown', kernel='poly',
           c=.5,d=2)
tmp$eval

gk = Gn(tmp$eval, n)
plot(gk)
which.max(gk)

pred=tmp$pred
vec=c(1,2)
plotpred(pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Nonlinear")
tmp$dim

pairs(pred[,1:4],col=y+1)

# 2ndRKHS: gauss 
plot(tmp$cv.shx)
max(tmp$cv.shx)
tmp$shx

# 2ndRKHS: poly
plot(tmp$cv.cd[-(139:140)])
max(tmp$cv.cd)
tmp$c;tmp$d

########################## 
# Classification
########################## 
set.seed(0)
train = sample(1:n, 100)
c.methods=c("lda", "qda", "svm")
d = tmp$dim
d=1
imethod=3

pred.class = pred
#pred.class = fpc.pred
out = classify(x=pred.class[train,1:d], y=y, x.test=pred.class[-train,1:d], method=c.methods[imethod])
mean(y[-train] != out)


# Linear FPCA
fpc = fpca(tmp$ftn)
fpc.pred = fpc$pred
pred.class = fpc.pred
out = classify(x=pred.class[train,1:d], y=y, x.test=pred.class[-train,1:d], method=c.methods[imethod])
mean(y[-train] != out)


# Linear FPCA (my code)
fpc = fpca(tmp$ftn)
fpc.pred = fpc$pred
vec=c(1,2)
plotpred(fpc.pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Linear/FPCA")


#####################
# Visualization
#####################

out=data.frame(PC1=pred[,1], PC2 = pred[,2], Y=as.factor(y))
pdf('nfpc2d.pdf')
ggplot(out,aes(x=PC1,y=PC2, colour=Y))  +
  geom_point()+theme_bw() +
  scale_colour_manual(name="Y",  
                      values = c("0"="black", "1"="red"))+
  theme(plot.title = element_text(hjust = 0.5,size=20),
        legend.position="right")+
  theme(plot.title = element_text(hjust = 0.5,size=20))+
  labs(title = paste0('NAFPC'))
dev.off()


# Linear FPCA
fpc = fpca(tmp$ftn)
fpc.pred = fpc$pred

out=data.frame(PC1=fpc.pred[,1], PC2 = fpc.pred[,2], Y=as.factor(y))
pdf('fpc2d.pdf')
ggplot(out,aes(x=PC1,y=PC2, colour=Y))  +
  geom_point()+theme_bw() +
  scale_colour_manual(name="Y",  
                      values = c("0"="black", "1"="red"))+
  theme(plot.title = element_text(hjust = 0.5,size=20),
        legend.position="right")+
  theme(plot.title = element_text(hjust = 0.5,size=20))+
  labs(title = paste0('FPC'))
dev.off()

# curves - visualization of PC Scores

curve_df <- data.frame(id=1, time=tt, x = x[1,], nfpc1=pred[1,1], nfpc2=pred[1,2], nfpc3=pred[1,3],
                       fpc1=fpc.pred[1,1], fpc2=fpc.pred[1,2], fpc3=fpc.pred[1,3], y=y[i])
for(i in 2:n){
  curve_df <- rbind.data.frame(curve_df, data.frame(id=i, time=tt, x = x[i,], nfpc1=pred[i,1], nfpc2=pred[i,2], nfpc3=pred[i,3],
                                                    fpc1=fpc.pred[i,1], fpc2=fpc.pred[i,2], fpc3=fpc.pred[i,3], y=y[i]))
}

curve_df$id <- factor(curve_df$id)
curve_df$y <- factor(curve_df$y)

require(ggplot2)
p <- ggplot(curve_df, aes(x=time, y=x,group=id,colour=y))
pp <- p + geom_line() + theme_bw() +
  scale_colour_manual(name="Y",  
                      values = c("0"="black", "1"="red"))+
  theme(plot.title = element_text(hjust = 0.5,size=20),
        legend.position="right")+
  theme(plot.title = element_text(hjust = 0.5,size=20))+
  labs(title = paste0('Model I-1'), y='X(t)')
pdf('mi_description.pdf')
pp
dev.off()

p <- ggplot(curve_df, aes(x=time, y=x,group=id))
pp <- p + geom_line() + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('Model I (Original)'), y= 'X(t)')
#pdf('m1.pdf')
pp
#dev.off()

hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(nfpc1,1), 0,0))
p <- ggplot(hex.df, aes(x=time, y=x,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="none") +
  labs(title = paste0('1st NAFPC'), y= 'X(t)')

# pdf('nfpc1.pdf')
pp
# dev.off()


hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(nfpc2,1), 0,0))
p <- ggplot(hex.df, aes(x=time, y=x,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('2nd NAFPC'), y= 'X(t)')

# pdf('nfpc2.pdf')
pp
# dev.off()




hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(fpc1,1), 0,0))
p <- ggplot(hex.df, aes(x=time, y=x,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('1st FPC'), y= 'X(t)')

pdf('fpc1.pdf')
pp
dev.off()



hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(fpc2,1), 0,0))
p <- ggplot(hex.df, aes(x=time, y=x,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('2nd FPC'), y= 'X(t)')

pdf('fpc2.pdf')
pp
dev.off()


hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(nfpc1,1), scaling(nfpc2,1), scaling(nfpc3,1)))
hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(fpc1,1), scaling(fpc2,1), scaling(fpc3,1)))
hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(nfpc1,1), scaling(nfpc2,1), 0))
hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(nfpc2,1), 0,0))


hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(nfpc1,1), scaling(nfpc2,1), scaling(nfpc3,1)))
p <- ggplot(hex.df, aes(x=time, y=x,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('First 3 NAFPC'))
pdf('nfpc123.pdf')
pp
dev.off()

hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(fpc1,1), scaling(fpc2,1), scaling(fpc3,1)))
p <- ggplot(hex.df, aes(x=time, y=x,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('First 3 FPC'))
pdf('fpc123.pdf')
pp
dev.off()

p <- ggplot(curve_df, aes(x=time, y=x,group=id,colour=fpc1))
pp <- p + geom_line() + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('1st (linear) FPC'))

pp


p <- ggplot(curve_df, aes(x=time, y=x,group=id,colour=fpc2))
pp <- p + geom_line() + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('2nd (linear) FPC'))

pp

hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(fpc1), scaling(fpc2), scaling(fpc3)))


p <- ggplot(hex.df, aes(x=time, y=x,group=id,colour=hex))
pp <- p + geom_line() + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('Sparsely observed'))

pp

#######################################################
# Model I-2  true function: linear -- just sum of a few eigenfunctions for each group
#######################################################
set.seed(0)
nt=20
n=200
nbasis=21
dat = sim.model.12(n=n, nt=nt)
x = dat$x; tt = dat$tt; y = dat$y

# Original Data
plot(tt,x[1,], type='l', ylim=c(-5,5), xlab="t", ylab="X(t)", main = "Model I", col=y)
for(i in 1:n)lines(tt,x[i,], col=(y[i]))

# Nonlinear FPCA

nafpca.m2=nafpca(x,tt, basisname="bspline", ex=0,nbasis=nbasis, gamma.tune=TRUE,ncv1=50,ncv2=50)  # shx works

nafpca.m2=nafpca(x,tt, shx=11, nbasis=30,gamma.tune = TRUE, ncv2=100,
           basisname='bspline')
nafpca.m2=nafpca(x,tt, nbasis=30,gamma.tune = TRUE, ncv1=30,ncv2=30,
           basisname='fourier')
nafpca.m2=nafpca(x,tt, shx=11, nbasis=30,gamma.tune = TRUE,ncv1=30,ncv2=30,
           basisname='gauss')
nafpca.m2=nafpca(x,tt, shx=11, nbasis=30,gamma.tune = TRUE,ncv1=30,ncv2=20,
           basisname='brown')
nafpca.m2=nafpca(x,tt, shx=11, nbasis=30,gamma.tune = TRUE,ncv1=30,ncv2=20,
           basisname='bspline', kernel='poly',
           c=.1,d=2)
nafpca.m2=nafpca(x,tt, shx=11, nbasis=30,gamma.tune = TRUE,
           basisname='fourier', kernel='poly',
           c=.5,d=2)
nafpca.m2=nafpca(x,tt, shx=11, nbasis=30,gamma.tune = TRUE,
           basisname='gauss', kernel='poly',
           c=.5,d=2)
nafpca.m2=nafpca(x,tt, shx=11, nbasis=30,gamma.tune = TRUE,
           basisname='brown', kernel='poly',
           c=.5,d=2)


pred=nafpca.m2$pred
pairs(pred[,1:4], col=y)
vec=c(1,2)
plotpred(pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Nonlinear")
ncv2=length(nafpca.m2$cv.shx)
shx.grid = c(exp(seq(log(10^(-10)),log(10^(2)),len=ncv2)))
plot(shx.grid, nafpca.m2$cv.shx)


plot(Gn(nafpca.m2$eval,n))
which.max(Gn(nafpca.m2$eval,n))
cumsum(nafpca.m2$eval)/sum(nafpca.m2$eval)
nafpca.m2$shx
nafpca.m2$d
nafpca.m2$cv.cd


# Linear FPCA (my code)

x.ftn = get.fd(xraw=t(x),tt=tt,basisname='bspline',nbasis=nbasis,ncv=10)
fpc = fpca(x.ftn)
fpc.pred = fpc$pred
vec=c(1,2)
plotpred(fpc.pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Linear/FPCA")



#####################
# Visualization
#####################
#y=y-1
out=data.frame(PC1=pred[,1], PC2 = pred[,2], Y=as.factor(y))
pdf('mi2_nfpc2d.pdf')
ggplot(out,aes(x=PC1,y=PC2, colour=Y))  +
  geom_point()+theme_bw() +
  scale_colour_manual(name="Y",  
                      values = c("1"="black", "2"="red"))+
  theme(plot.title = element_text(hjust = 0.5,size=20),
        legend.position="right")+
  theme(plot.title = element_text(hjust = 0.5,size=20))+
  labs(title = paste0('NAFPCA'))
dev.off()


# Linear FPCA
out=data.frame(PC1=fpc.pred[,1], PC2 = fpc.pred[,2], Y=as.factor(y))
pdf('mi2_fpc2d.pdf')
ggplot(out,aes(x=PC1,y=PC2, colour=Y))  +
  geom_point()+theme_bw() +
  scale_colour_manual(name="Y",  
                      values = c("1"="black", "2"="red"))+
  theme(plot.title = element_text(hjust = 0.5,size=20),
        legend.position="right")+
  theme(plot.title = element_text(hjust = 0.5,size=20))+
  labs(title = paste0('FPC'))
dev.off()

# curves - visualization of PC Scores

curve_df <- data.frame(id=1, time=tt, x = x[1,], nfpc1=pred[1,1], nfpc2=pred[1,2], nfpc3=pred[1,3],
                       fpc1=fpc.pred[1,1], fpc2=fpc.pred[1,2], fpc3=fpc.pred[1,3], y=y[i])
for(i in 2:n){
  curve_df <- rbind.data.frame(curve_df, data.frame(id=i, time=tt, x = x[i,], nfpc1=pred[i,1], nfpc2=pred[i,2], nfpc3=pred[i,3],
                                                    fpc1=fpc.pred[i,1], fpc2=fpc.pred[i,2], fpc3=fpc.pred[i,3], y=y[i]))
}

require(ggplot2)
curve_df$id <- factor(curve_df$id)
curve_df$y <- factor(curve_df$y)

p <- ggplot(curve_df, aes(x=time, y=x,group=id,colour=y))
pp <- p + geom_line() + theme_bw() +
  scale_colour_manual(name="Y",  
                      values = c("0"="black", "1"="red"))+
  theme(plot.title = element_text(hjust = 0.5,size=20),
        legend.position="right")+
  theme(plot.title = element_text(hjust = 0.5,size=20))+
  labs(title = paste0('Model I-2 with Y'), y='X(t)')
pdf('mi2_description.pdf')
pp
dev.off()

p <- ggplot(curve_df, aes(x=time, y=x,group=id))
pp <- p + geom_line() + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('Model I (Original)'), y= 'X(t)')
#pdf('mi2.pdf')
pp
#dev.off()

hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(nfpc1,1), 0,0))
p <- ggplot(hex.df, aes(x=time, y=x,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="none") +
  labs(title = paste0('1st NAFPC'), y= 'X(t)')

pdf('nfpc1.pdf')
pp
dev.off()


hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(nfpc2,1), 0,0))
p <- ggplot(hex.df, aes(x=time, y=x,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('2nd NAFPC'), y= 'X(t)')

pdf('nfpc2.pdf')
pp
dev.off()




hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(fpc1,1), 0,0))
p <- ggplot(hex.df, aes(x=time, y=x,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('1st FPC'), y= 'X(t)')

pdf('fpc1.pdf')
pp
dev.off()



hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(fpc2,1), 0,0))
p <- ggplot(hex.df, aes(x=time, y=x,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('2nd FPC'), y= 'X(t)')

pdf('fpc2.pdf')
pp
dev.off()





hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(nfpc1,1), scaling(nfpc2,1), scaling(nfpc3,1)))
hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(fpc1,1), scaling(fpc2,1), scaling(fpc3,1)))
# hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(nfpc1,1), scaling(nfpc2,1), 0))
# hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(nfpc2,1), 0,0))


hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(nfpc1,1), scaling(nfpc2,1), scaling(nfpc3,1)))
p <- ggplot(hex.df, aes(x=time, y=x,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('First 3 NAFPC'))
pdf('nfpc123.pdf')
pp
dev.off()

hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(fpc1,1), scaling(fpc2,1), scaling(fpc3,1)))
p <- ggplot(hex.df, aes(x=time, y=x,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('First 3 FPC'))
pdf('fpc123.pdf')
pp
dev.off()



######################
# Model II-1 (unbalanced)
######################
set.seed(1)
n=200
nt=20
dat = sim.model.21(n=n, nt=nt, balanced = FALSE)
x = dat$x; tt = dat$tt; y = dat$y

i=which(y==0)[1]
plot(x[[1]][[i]], x[[2]][[i]], xlim=c(-5,5),ylim=c(-5,5), pch=16)
for(i in which(y==0)[1:10]){
  points(x[[1]][[i]], x[[2]][[i]], xlim=c(-5,5),ylim=c(-5,5), pch=16)
}

for(j in which(y==1)[1:10]){
  points(x[[1]][[j]], x[[2]][[j]], col=2, pch=16)
}

# Nonlinear FPCA
tmp=nafpca(x,tt, p=2, nbasis=21,unbalanced=TRUE, ncv1=10, ncv2=10)
pred=tmp$pred
vec=c(1,2)
plotpred(pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Nonlinear")
plot(tmp$cv.shx)
tmp$shx
# Linear FPCA
fpc = fpca(tmp$ftn)
fpc.pred = fpc$pred
plotpred(fpc.pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Linear")
# Linear FPCA (MFPCA/ PACE)
out = convert.fd(fd=x,tt=tt, p=2,n=n)
mfpc = MFPCA(out, M = 5, uniExpansions = list(list(type = "splines1D", k = 10),
                                              list(type = "splines1D", k = 10)))
fpc.pred = mfpc$scores

plotpred(fpc.pred, vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Linear/PACE")



curve_df <- data.frame(id=1, time=tt, x1 = x[[1]][[1]], x2=x[[2]][[1]],nfpc1=pred[1,1], nfpc2=pred[1,2], nfpc3=pred[1,3],
                       fpc1=fpc.pred[1,1], fpc2=fpc.pred[1,2], fpc3=fpc.pred[1,3], y=y[i])
for(i in 2:n){
  curve_df <- rbind.data.frame(curve_df, data.frame(id=i, time=tt, x1 = x[[1]][[i]], x2=x[[2]][[i]], nfpc1=pred[i,1], nfpc2=pred[i,2], nfpc3=pred[i,3],
                                                    fpc1=fpc.pred[i,1], fpc2=fpc.pred[i,2], fpc3=fpc.pred[i,3], y=y[i]))
}

curve_df$id <- factor(curve_df$id)
curve_df$y <- factor(curve_df$y)

require(ggplot2)
p <- ggplot(curve_df, aes(x=x1, y=x2,group=id,colour=y))
pp <- p + geom_point() + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('Sparsely observed'))

pp

p <- ggplot(curve_df, aes(x=x1, y=x2,group=id,colour=nfpc1))
pp <- p + geom_point() + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('Sparsely observed'))

pp


p <- ggplot(curve_df, aes(x=x1, y=x2,group=id,colour=fpc1))
pp <- p + geom_point() + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('Sparsely observed'))

pp


######################
# Model II-1 (balanced)
######################
n=100
nt=10
dat = sim.model.21(n=n, nt=nt, balanced = TRUE)
x = dat$x; tt = dat$tt; y = dat$y
y

i=which(y==0)[1]
plot(x[i,,1], x[i,,2], col=1, xlim=c(-5,5),ylim=c(-5,5), pch=16)
for(i in which(y==0)[1:10]){
  points(x[i,,1], x[i,,2], col=1, xlim=c(-5,5),ylim=c(-5,5), pch=16)
}
for(j in which(y==1)[1:10]){
  points(x[j,,1], x[j,,2], col=2, pch=16)
}

# Nonlinear FPCA
tmp=nafpca(x,tt,nbasis=30,gamma.tune = TRUE, ncv2=100,
                 basisname='bspline')
tmp=nafpca(x,tt,nbasis=30,gamma.tune = TRUE, ncv2=100,
           basisname='fourier')
tmp=nafpca(x,tt,nbasis=30,gamma.tune = TRUE, ncv2=100,
           basisname='brown')
tmp=nafpca(x,tt,nbasis=30,gamma.tune = TRUE, ncv2=100,
           basisname='gauss')
tmp=nafpca(x,tt,nbasis=30,gamma.tune = TRUE, ncv2=100, kernel = 'poly', d=1,
           basisname='bspline')
tmp=nafpca(x,tt,nbasis=30,gamma.tune = TRUE, ncv2=100, kernel = 'poly', d=3,
           basisname='fourier')
tmp=nafpca(x,tt,nbasis=30,gamma.tune = TRUE, ncv2=100, kernel = 'poly', d=3,
           basisname='brown')
tmp=nafpca(x,tt,nbasis=30,gamma.tune = TRUE, ncv2=100, kernel = 'poly', d=3,
           basisname='gauss')
pred=tmp$pred
vec=c(1,2)
plotpred(pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Nonlinear")
#plot(tmp$cv.shx)
max(tmp$cv.shx)

gk = Gn(tmp$eval, n)
#plot(gk)
which.max(gk)


plot(tmp$cv.cd)
max(tmp$cv.cd)

# Linear FPCA
fpc = fpca(tmp$ftn)
fpc.pred = fpc$pred
plotpred(fpc.pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Linear")
# Linear FPCA (MFPCA/ PACE)
# out = convert.fd(fd=x,tt=tt, p=2,n=n)
# mfpc = MFPCA(out, M = 5, uniExpansions = list(list(type = "splines1D", k = 10),
#                                               list(type = "splines1D", k = 10)))
# fpc.pred = mfpc$scores
# 
# plotpred(fpc.pred, vec=vec,ind.inf=y,
#          xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
#          main="Linear/PACE")



curve_df <- data.frame(id=1, time=tt, x1 = x[[1]][[1]], x2=x[[2]][[1]],nfpc1=pred[1,1], nfpc2=pred[1,2], nfpc3=pred[1,3],
                       fpc1=fpc.pred[1,1], fpc2=fpc.pred[1,2], fpc3=fpc.pred[1,3], y=y[i])
for(i in 2:n){
  curve_df <- rbind.data.frame(curve_df, data.frame(id=i, time=tt, x1 = x[[1]][[i]], x2=x[[2]][[i]], nfpc1=pred[i,1], nfpc2=pred[i,2], nfpc3=pred[i,3],
                                                    fpc1=fpc.pred[i,1], fpc2=fpc.pred[i,2], fpc3=fpc.pred[i,3], y=y[i]))
}

curve_df$id <- factor(curve_df$id)
curve_df$y <- factor(curve_df$y)

require(ggplot2)
p <- ggplot(curve_df, aes(x=x1, y=x2,group=id,colour=y))
pp <- p + geom_point() + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('Sparsely observed'))

pp

p <- ggplot(curve_df, aes(x=x1, y=x2,group=id,colour=nfpc1))
pp <- p + geom_point() + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('Sparsely observed'))

pp


p <- ggplot(curve_df, aes(x=x1, y=x2,group=id,colour=fpc1))
pp <- p + geom_point() + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('Sparsely observed'))

pp

######################
# Model II-2 (unbalanced)
######################
n=200
nt=20
set.seed(0)
dat = sim.model.22(n=n, nt=nt, p=2, balanced=FALSE, m1=1, m2=2, r1=1, r2=0.7, sd1=.4, sd2=.2)
x = dat$x; tt = dat$tt; y = dat$y

i=which(y==0)[1]
plot(x[[1]][[i]], x[[2]][[i]], xlim=c(-5,5),ylim=c(-5,5), pch=16)
for(i in which(y==0)[1:10]){
  points(x[[1]][[i]], x[[2]][[i]], xlim=c(-5,5),ylim=c(-5,5), pch=16)
}

for(j in which(y==1)[1:10]){
  points(x[[1]][[j]], x[[2]][[j]], col=2, pch=16)
}

# Nonlinear FPCA
tmp=nafpca(x,tt, p=2, nbasis=8,unbalanced=TRUE,ncv1=20, ncv2=20, basisname='bspline', gamma.tune=TRUE,)
pred=tmp$pred
vec=c(1,2)
plotpred(pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Nonlinear")
plot(tmp$cv.shx)
# Linear FPCA
fpc = fpca(tmp$ftn)
fpc.pred = fpc$pred
plotpred(fpc.pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Linear")

# Linear FPCA (MFPCA/ PACE)
# out = convert.fd(fd=x,tt=tt, p=2,n=n)
# mfpc = MFPCA(out, M = 5, uniExpansions = list(list(type = "splines1D", k = 10),
#                                               list(type = "splines1D", k = 10)))
# fpc.pred = mfpc$scores
# 
# plotpred(fpc.pred, vec=vec,ind.inf=y,
#          xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
#          main="Linear/PACE")



curve_df <- data.frame(id=1, time=tt, x1 = x[[1]][[1]], x2=x[[2]][[1]],nfpc1=pred[1,1], nfpc2=pred[1,2], nfpc3=pred[1,3],
                       fpc1=fpc.pred[1,1], fpc2=fpc.pred[1,2], fpc3=fpc.pred[1,3], y=y[i])
for(i in 2:n){
  curve_df <- rbind.data.frame(curve_df, data.frame(id=i, time=tt, x1 = x[[1]][[i]], x2=x[[2]][[i]], nfpc1=pred[i,1], nfpc2=pred[i,2], nfpc3=pred[i,3],
                                                    fpc1=fpc.pred[i,1], fpc2=fpc.pred[i,2], fpc3=fpc.pred[i,3], y=y[i]))
}

curve_df$id <- factor(curve_df$id)
curve_df$y <- factor(curve_df$y)
curve_df$ry = rgb(scaling(as.numeric(curve_df$y),1),0,0)
require(ggplot2)
p <- ggplot(curve_df, aes(x=x1, y=x2,group=id,colour=y))
pp <- p + geom_point(aes(colour=curve_df$y)) + theme_bw() + 
  scale_colour_manual(name="Y",  
                      values = c("0"="black", "1"="red"))+
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="right")+
  labs(title = paste0('Model II with Y'),
       x =  expression({X^{1}}(t)),
       y =  expression({X^{2}}(t)))
pdf('m2.pdf')
pp
dev.off()

p <- ggplot(curve_df, aes(x=x1, y=x2,group=id))
pp <- p + geom_point() + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="right")+
  labs(title = paste0('Model II (Original)'),
       x =  expression({X^{1}}(t)),
       y =  expression({X^{2}}(t)))
pdf('m2_original.pdf')
pp
dev.off()

out = data.frame(PC1=pred[,1], PC2=pred[,2], Y=as.factor(y))
p <- ggplot(out, aes(x=PC1, y=PC2, colour=Y))
pp <- p + geom_point(aes(colour=Y)) + theme_bw() + 
  scale_colour_manual(name="Y",  
                      values = c("0"="black", "1"="red"))+
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="right")+
  labs(title = paste0('NAFPC'))
pdf('m2nfpc2d.pdf')
pp
dev.off()


out = data.frame(PC1=fpc.pred[,1], PC2=fpc.pred[,2], Y=as.factor(y))
p <- ggplot(out, aes(x=PC1, y=PC2, colour=Y))
pp <- p + geom_point(aes(colour=Y)) + theme_bw() + 
  scale_colour_manual(name="Y",  
                      values = c("0"="black", "1"="red"))+
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="right")+
  labs(title = paste0('FPC'))
pdf('m2fpc2d.pdf')
pp
dev.off()



p <- ggplot(curve_df, aes(x=x1, y=x2,group=id,colour=nfpc1))
pp <- p + geom_point(aes(colour=nfpc1)) + theme_bw() +scale_colour_gradient(low="black",high="red")+
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="right")+
  labs(title = paste0('Model II with NAFPC'),
                      x =  expression({X^{1}}(t)),
                      y =  expression({X^{2}}(t)))
pdf('m2nfpc1.pdf')
pp
dev.off()

p <- ggplot(curve_df, aes(x=x1, y=x2,group=id,colour=fpc1))
pp <- p + geom_point(aes(colour=fpc1)) + theme_bw() +scale_colour_gradient(low="black",high="red")+
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="right")+
  labs(title = paste0('Model II with FPC'),
                      x =  expression({X^{1}}(t)),
                      y =  expression({X^{2}}(t)))
pdf('m2fpc1.pdf')
pp
dev.off()

pp



hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(fpc1,1),0,0))
p <- ggplot(hex.df, aes(x=x1, y=x2,group=id,colour=hex))
pp <- p + geom_point(colour=hex.df$hex) + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('First 3 linear FPC'))
pp


hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(nfpc1,1), 0, 0))
p <- ggplot(hex.df, aes(x=x1, y=x2,group=id,colour=hex))
pp <- p + geom_point(colour=hex.df$hex) + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('First 3 linear FPC'))
pp




#############################################
# Model I-3 : classification
#######################################################
set.seed(0)
nt=20
n=200
nbasis=21
dat = sim.model.31(n=n, nt=nt)
x = dat$x; tt = dat$tt; y = dat$y

# Original Data
plot(tt,x[1,], type='l', ylim=c(min(x),max(x)), xlab="t", ylab="X(t)", main = "Model I", col=y)
for(i in 1:n)lines(tt,x[i,], col=(y[i]))

# Nonlinear FPCA

tmp=nafpca(x,tt, basisname="bspline",nbasis=nbasis, gamma.tune=TRUE,ncv1=20,ncv2=20)  # shx works
pred=tmp$pred
vec=c(1,2)
plotpred(pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Nonlinear")

tmp$dim

aa=fpca(nafpca.m2$ftn)

pairs(pred[,1:10],col=y)
pairs(aa$pred[,1:10],col=y)



######################
# Model II-3 
######################
set.seed(0)
n=200
nt=20
nbasis=21
dat = sim.model.23(n=n, nt=nt)
x = dat$x; tt = dat$tt; y = dat$y

fdplot(x,y,tt, n.point=100)
i=which(y==1)[1]
plot(x[i,,1], x[i,,2], col=1, xlim=c(-5,5),ylim=c(-5,5), pch=16)
for(i in which(y==1)[1:10]){
  points(x[i,,1], x[i,,2], col=1, xlim=c(-5,5),ylim=c(-5,5), pch=16)
}
for(j in which(y==2)[1:10]){
  points(x[j,,1], x[j,,2], col=2, pch=16)
}



# Nonlinear FPCA
tmp=nafpca(x,tt, p=2, nbasis=nbasis, ncv1= 30, ncv2=30)
pred=tmp$pred
vec=c(1,2)
plotpred(pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Nonlinear")

pairs(pred[,1:4], col=y)
fpca.out = fpca(tmp$ftn)
fpc.pred = fpca.out$pred
vec=c(1,2)
plotpred(fpc.pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Nonlinear")


######################
# Model III-2 
######################
set.seed(0)
n=200
nt=20
nbasis=21
dat = sim.model.32(n=n, nt=nt)
x = dat$x; tt = dat$tt; y = dat$y

fdplot(x,y,tt)

# Nonlinear FPCA
tmp=nafpca(x,tt, p=1, nbasis=nbasis, ncv1= 50, ncv2=50)
pred=tmp$pred
vec=c(1,2)
plotpred(pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Nonlinear")

pairs(pred[,1:4], col=y)
fpca.out = fpca(tmp$ftn)
fpc.pred = fpca.out$pred
vec=c(1,2)
plotpred(fpc.pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Nonlinear")


######################
# Model III-3
######################
set.seed(0)
n=200
nt=20
nbasis=21
dat = sim.model.33(n=n, nt=nt)
x = dat$x; tt = dat$tt; y = dat$y

fdplot(x,y,tt)

# Nonlinear FPCA
tmp=nafpca(x,tt, p=1, nbasis=nbasis, ncv1= 50, ncv2=50)
pred=tmp$pred
vec=c(3,4)
plotpred(pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Nonlinear")
pairs(pred[,1:6], col=y)

#pairs(pred[,1:4], col=y)
fpca.out = fpca(tmp$ftn)
fpc.pred = fpca.out$pred
vec=c(1,2)
plotpred(fpc.pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Nonlinear")

######################
# Model III-4
######################
set.seed(0)
n=200
nt=20
nbasis=31
dat = sim.model.34(n=n, nt=nt)
x = dat$x; tt = dat$tt; y = dat$y

fdplot(x,y,tt)

# Nonlinear FPCA
tmp=nafpca(x,tt, p=1, nbasis=nbasis, ncv1= 50, ncv2=50)
pred=tmp$pred
vec=c(1,2)
plotpred(pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Nonlinear")
pairs(pred[,1:6], col=y)

#pairs(pred[,1:4], col=y)
fpca.out = fpca(tmp$ftn)
fpc.pred = fpca.out$pred
vec=c(1,2)
plotpred(fpc.pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Nonlinear")

######################
# Model III-5
######################
set.seed(0)
n=200
nt=20
nbasis=31
dat = sim.model.35(n=n, nt=nt)
x = dat$x; tt = dat$tt; y = dat$y
y=as.factor(y)

# Nonlinear FPCA
tmp=nafpca(x,tt, p=3, nbasis=nbasis, ncv1= 50, ncv2=50)
pred=tmp$pred

fpca.out = fpca(tmp$ftn)
fpc.pred = fpca.out$pred

summary(glm(y~pred[,1:5], family = binomial(link="logit")))
summary(glm(y~fpc.pred[,1:5],family = binomial(link="logit")))

vec=c(1,2)
plotpred(pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Nonlinear")
plotpred(fpc.pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Linear")

i=1
curve_df <- data.frame(id=1, time=tt, x1 = x[1,,1], x2=x[1,,2], x3=x[1,,3],nfpc1=pred[1,1], nfpc2=pred[1,2], nfpc3=pred[1,3],
                       fpc1=fpc.pred[1,1], fpc2=fpc.pred[1,2], fpc3=fpc.pred[1,3], y=y[i])
for(i in 2:n){
  curve_df <- rbind.data.frame(curve_df, data.frame(id=i, time=tt, x1 = x[i,,1], x2=x[i,,2], x3=x[i,,3],nfpc1=pred[i,1], nfpc2=pred[i,2], nfpc3=pred[i,3],
                                                    fpc1=fpc.pred[i,1], fpc2=fpc.pred[i,2], fpc3=fpc.pred[i,3], y=y[i]))
}


p <- ggplot(curve_df, aes(x=nfpc1, y=nfpc2,group=id,colour=y))
pp <- p + geom_point(aes(colour=y)) + theme_bw() +scale_colour_gradient(low="black",high="red")+
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="right")+
  labs(title = paste0('Model II with NAFPC'),
       x =  expression({X^{1}}(t)),
       y =  expression({X^{2}}(t)))
pp

p <- ggplot(curve_df, aes(x=fpc1, y=fpc2,group=id,colour=y))
pp <- p + geom_point(aes(colour=y)) + theme_bw() +scale_colour_gradient(low="black",high="red")+
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="right")+
  labs(title = paste0('Model II with NAFPC'),
       x =  expression({X^{1}}(t)),
       y =  expression({X^{2}}(t)))
pp



hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(nfpc1,1), 0,0))
p <- ggplot(hex.df, aes(x=time, y=x3,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="none") +
  labs(title = paste0('1st NAFPC'), y= 'X(t)')

pp


hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(fpc3,1), 0,0))
p <- ggplot(hex.df, aes(x=time, y=x3,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="none") +
  labs(title = paste0('1st FPC'), y= 'X(t)')

pp

hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(nfpc1,1), scaling(nfpc2,1),scaling(nfpc3,1),))
p <- ggplot(hex.df, aes(x=time, y=x3,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="none") +
  labs(title = paste0('1st NAFPC'), y= 'X(t)')

pp

hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(fpc1,1), scaling(fpc2,1),scaling(fpc3,1),))
p <- ggplot(hex.df, aes(x=time, y=x3,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="none") +
  labs(title = paste0('1st NAFPC'), y= 'X(t)')

pp


######################
# Model III-6
######################
set.seed(0)
n=200
nt=20
nbasis=21
dat = sim.model.36(n=n, nt=nt)
x = dat$x; tt = dat$tt; y = dat$y
y=as.factor(y)

# Nonlinear FPCA
tmp=nafpca(x,tt, p=2, nbasis=nbasis, ncv1= 10, ncv2=10)
pred=tmp$pred

plot(x[6,,1])
fpca.out = fpca(tmp$ftn)
fpc.pred = fpca.out$pred

summary(glm(y~pred[,1:5], family = binomial(link="logit")))
summary(glm(y~fpc.pred[,1:5],family = binomial(link="logit")))

vec=c(1,2)
plotpred(pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Nonlinear")
plotpred(fpc.pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Linear")

i=1
curve_df <- data.frame(id=1, time=tt, x1 = x[1,,1], x2=x[1,,2], nfpc1=pred[1,1], nfpc2=pred[1,2], nfpc3=pred[1,3],
                       fpc1=fpc.pred[1,1], fpc2=fpc.pred[1,2], fpc3=fpc.pred[1,3], y=y[i])
for(i in 2:n){
  curve_df <- rbind.data.frame(curve_df, data.frame(id=i, time=tt, x1 = x[i,,1], x2=x[i,,2],nfpc1=pred[i,1], nfpc2=pred[i,2], nfpc3=pred[i,3],
                                                    fpc1=fpc.pred[i,1], fpc2=fpc.pred[i,2], fpc3=fpc.pred[i,3], y=y[i]))
}


p <- ggplot(curve_df, aes(x=x1, y=x2,group=id,colour=nfpc3))
pp <- p + geom_point(aes(colour=nfpc3)) + theme_bw() +#scale_colour_gradient(low="black",high="red")+
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="right")+
  labs(title = paste0('Model II with NAFPC'),
       x =  expression({X^{1}}(t)),
       y =  expression({X^{2}}(t)))
pp

p <- ggplot(curve_df, aes(x=fpc1, y=fpc2,group=id,colour=y))
pp <- p + geom_point(aes(colour=y)) + theme_bw() +scale_colour_gradient(low="black",high="red")+
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="right")+
  labs(title = paste0('Model II with NAFPC'),
       x =  expression({X^{1}}(t)),
       y =  expression({X^{2}}(t)))
pp



hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(nfpc1,1), 0,0))
p <- ggplot(hex.df, aes(x=time, y=x3,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="none") +
  labs(title = paste0('1st NAFPC'), y= 'X(t)')

pp


hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(fpc3,1), 0,0))
p <- ggplot(hex.df, aes(x=time, y=x3,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="none") +
  labs(title = paste0('1st FPC'), y= 'X(t)')

pp

hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(nfpc1,1), scaling(nfpc2,1),scaling(nfpc3,1),))
p <- ggplot(hex.df, aes(x=time, y=x3,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="none") +
  labs(title = paste0('1st NAFPC'), y= 'X(t)')

pp

hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(fpc1,1), scaling(fpc2,1),scaling(fpc3,1),))
p <- ggplot(hex.df, aes(x=time, y=x3,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="none") +
  labs(title = paste0('1st NAFPC'), y= 'X(t)')

pp

