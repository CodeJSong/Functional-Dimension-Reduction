################################
# Simulation 
################################
setwd("~/CloudStation/Codes/GitHub/Functional-Dimension-Reduction")
source('nonlinear/nonlinear.R')
################################################################
# X:Function, y:Euclidean   // with Gaussian Kernel
################################################################
set.seed(0)
# 3.2.1 Model II-1---------------------------
model.xfyr1 <- function(n, nt, N, gamma=1,sd=0.1){
  # beta(t) = 2k( - t1)+1k(-t2)+ 1k(-t3)
  t1 <- 0.1#;t2<-0.6;t3<-0.9
  t11 <- 0.5;t22<-0.5;t33<-0.8
  tt <- (1:nt)/nt
  ###################################
  # Generate X
  ###################################
  x <- NULL   # observed x 
  ttt <- NULL # it is for the function x
  a <- NULL
  for(i in 1:n){
    xi <- rep(0, nt)
    a[[i]] <- rnorm(N)
    ttt[[i]] <- runif(N)
    for(k in 1:N){
      xi <- xi + a[[i]][k] * exp(-gamma*(ttt[[i]][k]-tt)^2)
    }
    x[[i]] <- xi
  }
  # X=\sum_k a[[i]][k] k(.-ttt[[i]][k])
  
  ###################################
  # Generate Y depending on model   #
  ###################################
  b1 <- rep(0, n) # b1[i] <- <beta , X_i>
  b2 <- rep(0, n)
  true.pred <- matrix(0, 2, n)
  for(i in 1:n){
    b1[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t1)^2)) 
    b2[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t11)^2))# + sum(a[[i]]*exp(-(ttt[[i]]-t3)^2))
  }
  tt <- matrix(rep(tt, n),n,nt,byrow=T)
  y <- b1 + b2 + rnorm(n,sd=sd)
  beta.function <- function(t, k){
    if(k==1) return(exp(-(t-0.1)^2))
    if(k==2) return(exp(-(t-0.5)^2))
  }
  
  result <- list(true.pred=(b1+b2), beta.function=beta.function,true.beta = cbind(b1,b2), x=x, y=y, tt=tt)
}


################################################################
# Do nonlinear SDR
################################################################
# generate model
xfyr1 <- model.xfyr1(n=100, nt=10, N=5, gamma=7,sd=.1)
# do f-GSIR
xfyr1fir <- fdr(x=xfyr1$x,y=xfyr1$y,tx=tt,xftn=TRUE, xtequal=TRUE, 
                kernel="gaussian")
# do f-GSAVE
xfyr1fave <- fdr(x=xfyr1$x,y=xfyr1$y,tx=tt,xftn=TRUE, xtequal=TRUE, option="FAVE", 
                 kernel="gaussian")

plot(xfyr1$true.pred,xfyr1fir$pred[,1], 
     main="True vs. FIR", ylab="1st FIR predictor",
     xlab=paste("spearman ",round(spearman(xfyr1$true.pred,xfyr1fir$pred[,1]),4)))

plot(xfyr1$true.pred,xfyr1fave$pred[,1], 
     main="True vs. FAVE", ylab="1st FAVE predictor",
     xlab=paste("spearman ",round(spearman(xfyr1$true.pred,xfyr1fave$pred[,1]),4)))

