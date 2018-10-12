source('nafpca.R')
source('sim_model.R')
require(ggplot2)
require(dplyr)

scaling <- function(x,k) {
  #x = abs(x)
  out = (x - min(x)) / (max(x) - min(x)) / k
  return(out)
}
n=200; nt=20
set.seed(1)
###########
# Model I-1
###########
dat = sim.model.11(n=n, nt=nt)
x = dat$x; tt = dat$tt; y = dat$y

# Original Data
plot(tt,x[1,], type='l', ylim=c(-5,5), xlab="t", ylab="X(t)", main = "Model I", col=y+1)
for(i in 1:n)lines(tt,x[i,], col=(y[i]+1))

# Nonlinear FPCA
tmp=nafpca(x,tt, ex=0)
pred=tmp$pred
vec=c(1,2)
plotpred(pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Nonlinear FPC")

out=data.frame(PC1=pred[,1], PC2 = pred[,2], Y=as.factor(y))
ggplot(out,aes(x=PC1,y=PC2, colour=Y))  +
  geom_point()+theme_bw() +
  scale_colour_manual(name="Y",
                      values = c("0"="black", "1"="red"))+
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="right")+
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('Nonlinear FPC'))


# Linear FPCA
fpc = fpca(tmp$ftn)
fpc.pred = fpc$pred

out=data.frame(PC1=fpc.pred[,1], PC2 = fpc.pred[,2], Y=as.factor(y))

ggplot(out,aes(x=PC1,y=PC2, colour=Y))  +
  geom_point()+theme_bw() +
  scale_colour_manual(name="Y",
                      values = c("0"="black", "1"="red"))+
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="right")+
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('Linear FPC'))


# curves - visualization of PC Scores

curve_df <- data.frame(id=1, time=tt, x = x[1,], nfpc1=pred[1,1], nfpc2=pred[1,2], nfpc3=pred[1,3],
                       fpc1=fpc.pred[1,1], fpc2=fpc.pred[1,2], fpc3=fpc.pred[1,3], y=y[i])
for(i in 2:n){
  curve_df <- rbind.data.frame(curve_df, data.frame(id=i, time=tt, x = x[i,], nfpc1=pred[i,1], nfpc2=pred[i,2], nfpc3=pred[i,3],
                                                    fpc1=fpc.pred[i,1], fpc2=fpc.pred[i,2], fpc3=fpc.pred[i,3], y=y[i]))
}

curve_df$id <- factor(curve_df$id)
curve_df$y <- factor(curve_df$y)
p <- ggplot(curve_df, aes(x=time, y=x,group=id,colour=y))
pp <- p + geom_line() + theme_bw() +
  scale_colour_manual(name="Y",
                      values = c("0"="black", "1"="red"))+
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="right")+
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('Model I with Y'), y='X(t)')

pp

p <- ggplot(curve_df, aes(x=time, y=x,group=id))
pp <- p + geom_line() + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('Model I (Original)'), y= 'X(t)')

pp

#################################################
# color with principal components
#################################################
hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(nfpc1,1), 0,0))
p <- ggplot(hex.df, aes(x=time, y=x,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="none") +
  labs(title = paste0('1st Nonlinear FPC'), y= 'X(t)')

pp


hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(nfpc2,1), 0,0))
p <- ggplot(hex.df, aes(x=time, y=x,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('2nd Nonlinear FPC'), y= 'X(t)')

pp




hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(fpc1,1), 0,0))
p <- ggplot(hex.df, aes(x=time, y=x,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('1st Linear FPC'), y= 'X(t)')

pp



hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(fpc2,1), 0,0))
p <- ggplot(hex.df, aes(x=time, y=x,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('2nd Linear FPC'), y= 'X(t)')

pp


# hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(nfpc1,1), scaling(nfpc2,1), scaling(nfpc3,1)))
# hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(fpc1,1), scaling(fpc2,1), scaling(fpc3,1)))
# hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(nfpc1,1), scaling(nfpc2,1), 0))
# hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(nfpc2,1), 0,0))


hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(nfpc1,1), scaling(nfpc2,1), scaling(nfpc3,1)))
p <- ggplot(hex.df, aes(x=time, y=x,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('First 3 nonlinear FPC'))
pp

hex.df <-  curve_df %>%  mutate(hex = rgb(scaling(fpc1,1), scaling(fpc2,1), scaling(fpc3,1)))
p <- ggplot(hex.df, aes(x=time, y=x,group=id,colour=hex))
pp <- p + geom_line(colour=hex.df$hex) + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  labs(title = paste0('First 3 linear FPC'))
pp





######################
# Model II-2 (unbalanced)
######################
set.seed(1)
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
tmp=nafpca(x,tt, p=2)
pred=tmp$pred
vec=c(1,2)
plotpred(pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Nonlinear")

# Linear FPCA
fpc = fpca(tmp$ftn)
fpc.pred = fpc$pred
plotpred(fpc.pred,vec=vec,ind.inf=y,
         xlab=paste("PC ",vec[1]), ylab=paste("PC ",vec[2]),
         main="Linear")



curve_df <- data.frame(id=1, time=tt, x1 = x[[1]][[1]], x2=x[[2]][[1]],nfpc1=pred[1,1], nfpc2=pred[1,2], nfpc3=pred[1,3],
                       fpc1=fpc.pred[1,1], fpc2=fpc.pred[1,2], fpc3=fpc.pred[1,3], y=y[i])
for(i in 2:n){
  curve_df <- rbind.data.frame(curve_df, data.frame(id=i, time=tt, x1 = x[[1]][[i]], x2=x[[2]][[i]], nfpc1=pred[i,1], nfpc2=pred[i,2], nfpc3=pred[i,3],
                                                    fpc1=fpc.pred[i,1], fpc2=fpc.pred[i,2], fpc3=fpc.pred[i,3], y=y[i]))
}

curve_df$id <- factor(curve_df$id)
curve_df$y <- factor(curve_df$y)
curve_df$ry = rgb(scaling(as.numeric(curve_df$y),1),0,0)

p <- ggplot(curve_df, aes(x=x1, y=x2,group=id,colour=y))
pp <- p + geom_point(aes(colour=curve_df$y)) + theme_bw() +
  scale_colour_manual(name="Y",
                      values = c("0"="black", "1"="red"))+
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="right")+
  labs(title = paste0('Model II (with Y)'),
       x =  expression({X^{1}}(t)),
       y =  expression({X^{2}}(t)))
pp

p <- ggplot(curve_df, aes(x=x1, y=x2,group=id))
pp <- p + geom_point() + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="right")+
  labs(title = paste0('Model II (Original)'),
       x =  expression({X^{1}}(t)),
       y =  expression({X^{2}}(t)))
pp

out = data.frame(PC1=pred[,1], PC2=pred[,2], Y=as.factor(y))
p <- ggplot(out, aes(x=PC1, y=PC2, colour=Y))
pp <- p + geom_point(aes(colour=Y)) + theme_bw() +
  scale_colour_manual(name="Y",
                      values = c("0"="black", "1"="red"))+
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="right")+
  labs(title = paste0('Nonlinear FPC'))
pp


out = data.frame(PC1=fpc.pred[,1], PC2=fpc.pred[,2], Y=as.factor(y))
p <- ggplot(out, aes(x=PC1, y=PC2, colour=Y))
pp <- p + geom_point(aes(colour=Y)) + theme_bw() +
  scale_colour_manual(name="Y",
                      values = c("0"="black", "1"="red"))+
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="right")+
  labs(title = paste0('Linear FPC'))
pp



p <- ggplot(curve_df, aes(x=x1, y=x2,group=id,colour=nfpc1))
pp <- p + geom_point(aes(colour=nfpc1)) + theme_bw() +scale_colour_gradient(low="black",high="red")+
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="right")+
  labs(title = paste0('Model II (with Nonlinear FPC)'),
                      x =  expression({X^{1}}(t)),
                      y =  expression({X^{2}}(t)))
pp

p <- ggplot(curve_df, aes(x=x1, y=x2,group=id,colour=fpc1))
pp <- p + geom_point(aes(colour=fpc1)) + theme_bw() +scale_colour_gradient(low="black",high="red")+
  theme(plot.title = element_text(hjust = 0.5,size=22),
        legend.position="right")+
  labs(title = paste0('Model II (with Linear FPC)'),
                      x =  expression({X^{1}}(t)),
                      y =  expression({X^{2}}(t)))
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
  labs(title = paste0('First 3 nonlinear FPC'))
pp
