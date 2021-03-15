
################################################################
# 3.2. X:Function, y:Euclidean   // with Kernel Thing
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

# 3.2.1 Model II-2---------------------------
model.xfyr2 <- function(n, nt, N, gamma=1, sd=0.1){
  # beta1(t) = k( - t1)#+1k(-t2)+ 1k(-t3)
  # beta2(t) = k(  - t11)#+ 1k(-t22) + 1k(-t33)
  t1 <- 0.1#;t2<-0.6;t3<-0.9
  t11 <- 0.5#;t22<-0.5;t33<-0.8
  t111 <- 0.9
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
  b1 <- rep(0, n) # b1[i] <- <beta_1 , X_i>
  b2 <- rep(0, n) # b2[i] <- <beta_2 , X_i>
  b3 <- rep(0, n) # b3[i] <- <beta_3 , X_i>
  b1sq <- rep(0, n) # b1sq[i] <- <beta_1 , X_i^2>
  b2sq <- rep(0, n) # b2sq[i] <- <beta_2 , X_i^2>
  true.pred <- matrix(0, 2, n)
  for(i in 1:n){
    b1[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t1)^2)) #+ sum(a[[i]]*exp(-(ttt[[i]]-t2)^2)) + sum(a[[i]]*exp(-(ttt[[i]]-t3)^2))
    b2[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t11)^2)) #+ sum(a[[i]]*exp(-(ttt[[i]]-t22)^2)) + sum(a[[i]]*exp(-(ttt[[i]]-t33)^2))
    b3[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t111)^2)) #+ sum(a[[i]]*exp(-(ttt[[i]]-t22)^2)) + sum(a[[i]]*exp(-(ttt[[i]]-t33)^2))
    #b1sq[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t1)^2))^2# + sum(a[[i]]*exp(-(ttt[[i]]-t2)^2))^2 + sum(a[[i]]*exp(-(ttt[[i]]-t3)^2))^2
    #b2sq[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t11)^2))^2 #+ sum(a[[i]]*exp(-(ttt[[i]]-t22)^2))^2 + sum(a[[i]]*exp(-(ttt[[i]]-t33)^2))^2
  }
  tt <- matrix(rep(tt, n),n,nt,byrow=T)
  y <- b1/(1+exp(b2)) + b3 + rnorm(n,sd=sd)
  beta.function <- function(t, k){
    if(k==1) return(exp(-(t-0.1)^2))
    if(k==2) return(exp(-(t-0.5)^2))
    if(k==3) return(exp(-(t-0.9)^2))
  }
  
  result <- list(true.pred=b1/(1+exp(b2)) +b3, true.beta= cbind(b1,b2,b3),beta.function=beta.function,x=x, y=y, tt=tt)
}


# 3.2.1 Model II-3---------------------------
model.xfyr3 <- function(n, nt, N, gamma=1, sd=0.1){
  # beta1(t) = k( - t1)+1k(-t2)+ 1k(-t3)
  # beta2(t) = k(  - t11)+ 1k(-t22) + 1k(-t33)
  t1 <- 0.1#;t2<-0.6;t3<-0.9
  t11 <- 0.9#;t22<-0.5;t33<-0.8
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
  b1 <- rep(0, n) # b1[i] <- <beta_1 , X_i>
  b2 <- rep(0, n) # b2[i] <- <beta_2 , X_i>
  b1sq <- rep(0, n) # b1sq[i] <- <beta_1 , X_i^2>
  b2sq <- rep(0, n) # b2sq[i] <- <beta_2 , X_i^2>
  true.pred <- matrix(0, 2, n)
  for(i in 1:n){
    b1[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t1)^2)) #+ sum(a[[i]]*exp(-(ttt[[i]]-t2)^2)) + sum(a[[i]]*exp(-(ttt[[i]]-t3)^2))
    b2[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t11)^2)) #+ sum(a[[i]]*exp(-(ttt[[i]]-t22)^2)) + sum(a[[i]]*exp(-(ttt[[i]]-t33)^2))
    b1sq[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t1)^2))^2 #+ sum(a[[i]]*exp(-(ttt[[i]]-t2)^2))^2 + sum(a[[i]]*exp(-(ttt[[i]]-t3)^2))^2
    b2sq[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t11)^2))^2 #+ sum(a[[i]]*exp(-(ttt[[i]]-t22)^2))^2+ sum(a[[i]]*exp(-(ttt[[i]]-t33)^2))^2
  }
  tt <- matrix(rep(tt, n),n,nt,byrow=T)
  #y <- log(b1sq+(1+exp(b2)))+ rnorm(n,sd=sd)
  y <- (b1 + b2) * rnorm(n, sd=sd)
  beta.function <- function(t, k){
    if(k==1) return(exp(-(t-0.1)^2))
    if(k==2) return(exp(-(t-0.9)^2))
  }
  result <- list(true.pred=(b1+b2), beta.function=beta.function,true.beta = cbind(b1,b2),x=x, y=y, tt=tt)
}

model.xfyr2.extra <- function(n, nt, N, gamma=1, sd=0.1){
  # beta1(t) = k( - t1)+1k(-t2)+ 1k(-t3)
  # beta2(t) = k(  - t11)+ 1k(-t22) + 1k(-t33)
  t1 <- 0.1#;t2<-0.6;t3<-0.9
  t11 <- 0.5#;t22<-0.5;t33<-0.8
  tt <- (1:nt)/nt
  ###################################
  # Generate X
  ###################################
  x <- NULL   # observed x 
  ttt <- NULL # it is for the function x
  a <- NULL
  xsq <- NULL
  for(i in 1:n){
    xi <- rep(0, nt)
    a[[i]] <- rnorm(N)
    ttt[[i]] <- runif(N)
    for(k in 1:N){
      xi <- xi + a[[i]][k] * exp(-gamma*(ttt[[i]][k]-tt)^2)
    }
    x[[i]] <- xi
    xsq[[i]] <- xi^2
  }
  xsqmat <- makematrix(xsq,nt)
  KT <- exp(-gamma*normmat(as.matrix(tt)))
  KT.inv <- matpower2(KT, -1, 0.3)
  if(t1 %in% tt) {
    b1sq <- t(KT[which(t1==tt),] %*% KT.inv %*% KT %*% KT.inv %*% t(xsqmat))
  }
  if(!(t1 %in% tt)) stop('need coding')
  # X=\sum_k a[[i]][k] k(.-ttt[[i]][k])
  
  ###################################
  # Generate Y depending on model   #
  ###################################
  b2 <- rep(0, n) # b1[i] <- <beta_1 , X_i>
  true.pred <- matrix(0, 2, n)
  for(i in 1:n){
    b2[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t11)^2)) #+ sum(a[[i]]*exp(-(ttt[[i]]-t22)^2)) + sum(a[[i]]*exp(-(ttt[[i]]-t33)^2))
  }
  tt <- matrix(rep(tt, n),n,nt,byrow=T)
  #y <- log(b1sq+(1+exp(b2)))+ rnorm(n,sd=sd)
  y <- log(b1sq + exp(b2) + 1 )+ rnorm(n, sd=sd)
  beta.function <- function(t, k){
    if(k==1) return(exp(-(t-0.1)^2))
    if(k==2) return(exp(-(t-0.9)^2))
  }
  result <- list(true.pred=(b1sq+exp(b2)+1),x=x, y=y, tt=tt)
}

model.xfyr1.extra <- function(n, nt, N, gamma=1, sd=0.1){
  # beta1(t) = k( - t1)#+1k(-t2)+ 1k(-t3)
  # beta2(t) = k(  - t11)#+ 1k(-t22) + 1k(-t33)
  t1 <- 0.1#;t2<-0.6;t3<-0.9
  t11 <- 0.5#;t22<-0.5;t33<-0.8
  t111 <- 0.9
  tt <- (1:nt)/nt
  KT <- exp(-gamma*normmat(as.matrix(tt)))
  KT.inv <- matpower2(KT, -1, 0.3)
  ###################################
  # Generate X
  ###################################
  x <- NULL   # observed x 
  xnormsq<-rep(0,n)
  ttt <- NULL # it is for the function x
  a <- NULL
  for(i in 1:n){
    xi <- rep(0, nt)
    xnormsqi <- 0
    a[[i]] <- rnorm(N)
    ttt[[i]] <- runif(N)
    for(k in 1:N){
      xi <- xi + a[[i]][k] * exp(-gamma*(ttt[[i]][k]-tt)^2)
    }
    x[[i]] <- xi
    xnormsq[i] <- t(xi) %*% KT.inv %*% KT %*% KT.inv %*% xi
  }
  # X=\sum_k a[[i]][k] k(.-ttt[[i]][k])
  # xnorm =  \sum_jk  a[[i]][j][k] k(ttt[[i]][j] - ttt[[i]][k] )
  ###################################
  # Generate Y depending on model   #
  ###################################
  b1 <- rep(0, n) # b1[i] <- <beta_1 , X_i>
  b2 <- rep(0, n) # b2[i] <- <beta_2 , X_i>
  true.pred <- matrix(0, 2, n)
  for(i in 1:n){
    b1[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t1)^2)) #+ sum(a[[i]]*exp(-(ttt[[i]]-t2)^2)) + sum(a[[i]]*exp(-(ttt[[i]]-t3)^2))
    b2[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t11)^2)) #+ sum(a[[i]]*exp(-(ttt[[i]]-t22)^2)) + sum(a[[i]]*exp(-(ttt[[i]]-t33)^2))
  }
  tt <- matrix(rep(tt, n),n,nt,byrow=T)
  y <- b1/(1+exp(b2)) + 0.2* xnormsq + rnorm(n,sd=sd)
  beta.function <- function(t, k){
    if(k==1) return(exp(-(t-0.1)^2))
    if(k==2) return(exp(-(t-0.5)^2))
    if(k==3) return(exp(-(t-0.9)^2))
  }
  
  result <- list(true.pred=b1/(1+exp(b2)) +0.2 *xnormsq, xnorm=xnormsq,x=x, y=y, tt=tt)
}

# model.brown.pca-----
set.seed(0)
# Functional principal component of Wiener Process 
fpcawiener <- function(t, j){
  sqrt(2) * sin((j-1/2)*pi*t)
}

generate.fpca <- function(n, nt, p, sd=1){
  nnt <- nt * 5
  tt <-(1:nnt)/nnt
  BM <- matrix(0, n, nnt)
  a <- matrix(rnorm(n*p,sd=sd),n,p)
  for(j in 1: p) {
    a[,j] <- a[,j]/((j-1/2)*pi)
    for(i in 1:n){
      BM[i,] <- BM[i,] + a[i,j]* fpcawiener(tt,j)
    }    
  }
  #unbalanced
  x.ub <- matrix(0, n, nt)
  t.ub <- matrix(0, n, nt)
  for(i in 1:n){
    t.ub[i, ] <- sort(sample(1:nnt, nt))
    x.ub[i, ] <- BM[i, t.ub[i,]]
  }
  t.ub <- t.ub / nnt
  #balanced
  t.b <- round(seq(1, nnt, len=nt))
  x.b <- BM[,t.b]
  t.b <- matrix(t.b, nrow=n, ncol=nt, byrow=T)/ nnt
  result <- list(x.ub=x.ub, x.b=x.b, t.ub=t.ub, t.b=t.b, a=a)
}
###################################
# Generate Y depending on model   #
###################################
model.pcaxfyr.1 <- function(n, nt, p, sd=1,ssq=0.1){
  result <- generate.fpca(n=n, nt=nt, p=p, sd=sd)
  a <- result$a
  y.ub <- rep(0, n)
  y.b <- rep(0, n)
  b1 <- matrix(0, 2, nt)  # row 1 : unbalnced, row 2 : balanced
  b2 <- matrix(0, 2, nt)  # row 1 : unbalnced, row 2 : balanced
  true.pred <- matrix(0, 2, n)
  for(i in 1:n){
    #b1[1, ] <- fpcawiener(t.ub[i,],1)
    #b1[2, ] <- fpcawiener(t.b[i,],1)
    true.pred[1,i] <- (a[i,1]) + a[i,2]
    true.pred[2,i] <- (a[i,1]) + a[i,2]
    
    y.ub[i] <- (true.pred[1,i]) + ssq*rnorm(1)
    y.b[i]  <- (true.pred[2,i]) + ssq*rnorm(1)
  }
  result$true.pred <- true.pred
  result$y.ub <- y.ub
  result$y.b <- y.b
  return(result)
}
model.pcaxfyr.2 <- function(n, nt, p, sd=1,ssq=0.1){
  result <- generate.fpca(n=n, nt=nt, p=p, sd=sd)
  a <- result$a
  y.ub <- rep(0, n)
  y.b <- rep(0, n)
  b1 <- matrix(0, 2, nt)  # row 1 : unbalnced, row 2 : balanced
  b2 <- matrix(0, 2, nt)  # row 1 : unbalnced, row 2 : balanced
  true.pred <- matrix(0, 2, n)
  for(i in 1:n){
    #b1[1, ] <- fpcawiener(t.ub[i,],1)
    #b1[2, ] <- fpcawiener(t.b[i,],1)
    true.pred[1,i] <- a[i,1]/(1+exp(a[i,2])) + a[i,3]
    true.pred[2,i] <- a[i,1]/(1+exp(a[i,2])) + a[i,3]
    
    y.ub[i] <- (true.pred[1,i]) + ssq*rnorm(1)
    y.b[i]  <- (true.pred[2,i]) + ssq*rnorm(1)
  }
  result$true.pred <- true.pred
  result$y.ub <- y.ub
  result$y.b <- y.b
  return(result)
}

model.pcaxfyr.3 <- function(n, nt, p, sd=1,ssq=0.1){
  result <- generate.fpca(n=n, nt=nt, p=p, sd=sd)
  a <- result$a
  y.ub <- rep(0, n)
  y.b <- rep(0, n)
  b1 <- matrix(0, 2, nt)  # row 1 : unbalnced, row 2 : balanced
  b2 <- matrix(0, 2, nt)  # row 1 : unbalnced, row 2 : balanced
  true.pred <- matrix(0, 2, n)
  for(i in 1:n){
    #b1[1, ] <- fpcawiener(t.ub[i,],1)
    #b1[2, ] <- fpcawiener(t.b[i,],1)
    true.pred[1,i] <- (a[i,1] + a[i,2])
    true.pred[2,i] <- (a[i,1] + a[i,2])
#     true.pred[1,i] <- (a[i,1])*(1+exp(a[i,2]))
#     true.pred[2,i] <- (a[i,1])*(1+exp(a[i,2]))
    
    y.ub[i] <- (true.pred[1,i]) * ssq*rnorm(1)
    y.b[i]  <- (true.pred[2,i]) * ssq*rnorm(1)
  }
  result$true.pred <- true.pred
  result$y.ub <- y.ub
  result$y.b <- y.b
  return(result)
}

model.pcaxfyr.2.extra <- function(n, nt, p, sd=1,ssq=0.1){
  result <- generate.fpca(n=n, nt=nt, p=p, sd=sd)
  xsq.ub <- (result$x.ub)^2;xsq.b <- (result$x.b)^2;
  tub.unique <- sort(unique(as.vector(result$t.ub)))
  tb.unique <- sort(unique(as.vector(result$t.b)))
  ntub <- length(tub.unique);ntb <- length(tb.unique)
  KTub <- matrix(0, ntub, ntub)
  KTb <- matrix(0, ntb, ntb)
  for(i in 1:ntub)  KTub[i,] <- c(tub.unique[1:i],rep(tub.unique[i],(ntub-i)))
  for(i in 1:ntb)  KTb[i,] <- c(tb.unique[1:i],rep(tb.unique[i],(ntb-i)))
  
  KTb.inv <- matpower2(KTb, -1, 0.3)
  KTub.inv <- matpower2(KTub, -1, 0.3)
  b1b <- fpcawiener(tb.unique,1);b1ub <- fpcawiener(tub.unique,1);
  b1sq.b <- t(t(b1b) %*% KTb.inv %*% KTb %*% KTb.inv %*% t(xsq.b))
  b1sq.ub <- rep(0,n)
  
  a <- result$a
  y.ub <- rep(0, n)
  y.b <- rep(0, n)
  b1 <- matrix(0, 2, nt)  # row 1 : unbalnced, row 2 : balanced
  b2 <- matrix(0, 2, nt)  # row 1 : unbalnced, row 2 : balanced
  true.pred <- matrix(0, 2, n)
  for(i in 1:n){
    index.i <- rep(0,nt)
    for(k in 1:nt){
      index.i[k] <- which(result$t.ub[i,k]==tub.unique)
    }
    KTi <- KTub[,index.i]
    KTiub.inv <- matpower2(KTub[index.i,index.i], -1, 0.3)
    b1sq.ub[i] <- t(b1ub) %*% KTub.inv %*% KTi %*% KTiub.inv %*% (xsq.ub[i,])
    true.pred[1,i] <- b1sq.ub[i] + exp(a[i,2])
    true.pred[2,i] <- b1sq.b[i] + exp(a[i,2])
    y.ub[i] <- log(true.pred[1,i] + 1) + ssq*rnorm(1)
    y.b[i]  <- log(true.pred[2,i] + 1) + ssq*rnorm(1)
  }
  result$true.pred <- true.pred
  result$y.ub <- y.ub
  result$y.b <- y.b
  return(result)
}

model.pcaxfyr.1.extra <- function(n, nt, p, sd=1,ssq=0.1){
  result <- generate.fpca(n=n, nt=nt, p=p, sd=sd)
  x.b <- (result$x.b)
  x.ub <- result$x.ub
  tub.unique <- sort(unique(as.vector(result$t.ub)))
  tb.unique <- sort(unique(as.vector(result$t.b)))
  ntub <- length(tub.unique);ntb <- length(tb.unique)
  KTub <- matrix(0, ntub, ntub)
  KTb <- matrix(0, ntb, ntb)
  for(i in 1:ntub)  KTub[i,] <- c(tub.unique[1:i],rep(tub.unique[i],(ntub-i)))
  for(i in 1:ntb)  KTb[i,] <- c(tb.unique[1:i],rep(tb.unique[i],(ntb-i)))
  
  KTb.inv <- matpower2(KTb, -1, 0.3)
  KTub.inv <- matpower2(KTub, -1, 0.3)
  b1b <- fpcawiener(tb.unique,1);b1ub <- fpcawiener(tub.unique,1);
  xnormsq.b <- diag(t(x.b %*% KTb.inv %*% KTb %*% KTb.inv %*% t(x.b)))
  xnormsq.ub <- rep(0,n)
  
  a <- result$a
  y.ub <- rep(0, n)
  y.b <- rep(0, n)
  b1 <- matrix(0, 2, nt)  # row 1 : unbalnced, row 2 : balanced
  b2 <- matrix(0, 2, nt)  # row 1 : unbalnced, row 2 : balanced
  true.pred <- matrix(0, 2, n)
  for(i in 1:n){
    index.i <- rep(0,nt)
    for(k in 1:nt){
      index.i[k] <- which(result$t.ub[i,k]==tub.unique)
    }
    KTi <- KTub[index.i,index.i]
    KTiub.inv <- matpower2(KTub[index.i,index.i], -1, 0.3)
    xnormsq.ub[i] <- t(x.ub[i,]) %*% KTiub.inv %*% KTi %*% KTiub.inv %*% (x.ub[i,])
    true.pred[1,i] <- 0.2 * xnormsq.ub[i] +  a[i,1]/(1+exp(a[i,2]))
    true.pred[2,i] <- 0.2 * xnormsq.b[i] +  a[i,1]/(1+exp(a[i,2]))
    y.ub[i] <- true.pred[1,i] + ssq*rnorm(1)
    y.b[i]  <- true.pred[2,i] + ssq*rnorm(1)
  }
  result$true.pred <- true.pred
  result$y.ub <- y.ub
  result$y.b <- y.b
  result$xnormsq<- xnormsq.b
  return(result)
}

################################################################
# 3.3. X:Euclidean, Y:Function-----
################################################################
mu_gauss <- function(fun,t){
  if(fun==0) mu <- exp(-(t-0.1)^2)
  if(fun==1) mu <- exp(-(t-0.3)^2)
  if(fun==2) mu <- exp(-(t-0.5)^2)
  if(fun==3) mu <- exp(-(t-0.7)^2)
  if(fun==4) mu <- exp(-(t-0.9)^2)
  return(mu)
}
mu_brown <- function(fun,t){
  if(fun==0)  mu <- cos(tt)
}
f <- function(x, k){  
  # x is p dimensional vector
  # fk : R^p -> R
  if(k==1) fk <- x[1] + x[2]
  if(k==2) fk <-   (1+x[1])/(1+ exp(x[2])) + x[3]
  if(k==4) fk <- cos(x[1])+sin(x[2])
  return(fk)
}
# 3.3.1 Model -------
xryf1 <- function(n, nt, p=10, sd,pbm=100,  ssq){
  tt <- seq(0,1,len=nt)
  BM <- matrix(0, n, nt)
  a <- matrix(rnorm(n*pbm,sd=sd),n,pbm)
  for(j in 1: pbm) {
    a[,j] <- a[,j]/((j-1/2)*pi)
    for(i in 1:n){
      BM[i,] <- BM[i,] + a[i,j]* fpcawiener(tt,j)
    }    
  }
  x <- matrix(rnorm(n*p,sd=sd),n,p)
  tt <- (1:nt)/nt
  y <- NULL
  true.pred <- rep(0, n)
  for(i in 1:n){
    true.pred[i] <- f(x[i,],1)
    y[[i]] <-  sin(true.pred[i] *tt^2 * pi) + ssq * BM[i,]
  }
  tt <- matrix(rep(tt, n),n,nt,byrow=T)
  return(list(y=y, x=x, true.pred=true.pred,tt=tt))
}
xryf1.new <- function(n, nt, p=10, mu, sd, N=7, gamma=7,sigma=.1){  #sd is for x , sigma is for Y
  
  t1 <- 0.2; t2 <- 0.5; t3 <- 0.8; t4<- 0.7; t5<-0.9#;t2<-0.6;t3<-0.9
  t11 <- 0.5;t22<-0.5;t33<-0.8
  tt <- (1:nt)/nt
  ###################################
  # Generate X
  ###################################
  
  x <- matrix(rnorm(n*p,mean=mu, sd=sd),n,p)
  true.pred <- rep(0,n)
  tt <- (1:nt)/nt
  y <- NULL
  ttt <- NULL # it is for the function x
  a <- NULL
  xx <- NULL # which is for epsilon(t)
  for(i in 1:n){
    xi <- rep(0, nt)
    a[[i]] <- rnorm(N)
    ttt[[i]] <- runif(N)
    for(k in 1:N){
      xi <- xi + a[[i]][k] * exp(-gamma*(ttt[[i]][k]-tt)^2)
    }
    xx[[i]] <- xi  ## epsilon(t)
    true.pred[i] <- f(x[i,],1) ;
    y[[i]] <-  true.pred[i] * ( exp(-gamma*(t1-tt)^2) - exp(-gamma*(t2-tt)^2) + exp(-gamma*(t3-tt)^2)) + sigma * xx[[i]]
  }
  # XX=\sum_k a[[i]][k] k(.-ttt[[i]][k])
  tt <- matrix(rep(tt, n),n,nt,byrow=T)
  return(list(y=y, x=x, true.pred=true.pred,tt=tt))
}
xryf2 <- function(n, nt, p=10, sd,pbm=100,  ssq){
  tt <- seq(0,1,len=nt)
  BM <- matrix(0, n, nt)
  a <- matrix(rnorm(n*pbm,sd=sd),n,pbm)
  for(j in 1: pbm) {
    a[,j] <- a[,j]/((j-1/2)*pi)
    for(i in 1:n){
      BM[i,] <- BM[i,] + a[i,j]* fpcawiener(tt,j)
    }    
  }
  x <- matrix(rnorm(n*p,sd=sd),n,p)
  tt <- (1:nt)/nt
  y <- NULL
  true.pred <- rep(0, n)
  for(i in 1:n){
    true.pred[i] <- f(x[i,],2)
    y[[i]] <-   true.pred[i] * tt^2 * exp(-tt)  + ssq * BM[i,]
  }
  tt <- matrix(rep(tt, n),n,nt,byrow=T)
  return(list(y=y, x=x, true.pred=true.pred,tt=tt))
}
xryf2.new <- function(n, nt, p=10, mu,sd, N=7, gamma=7,sigma=.1){  #sd is for x , sigma is for Y
  
  t1 <- 0.2; t2 <- 0.5; t3 <- 0.8#;t2<-0.6;t3<-0.9
  t11 <- 0.5;t22<-0.5;t33<-0.8
  tt <- (1:nt)/nt
  ###################################
  # Generate X
  ###################################
  
  x <- matrix(rnorm(n*p,mean=mu, sd=sd),n,p)
  true.pred <- rep(0, n)
  tt <- (1:nt)/nt
  y <- NULL
  ttt <- NULL # it is for the function x
  a <- NULL
  xx <- NULL # which is for epsilon(t)
  for(i in 1:n){
    xi <- rep(0, nt)
    a[[i]] <- rnorm(N)
    ttt[[i]] <- runif(N)
    for(k in 1:N){
      xi <- xi + a[[i]][k] * exp(-gamma*(ttt[[i]][k]-tt)^2)
    }
    xx[[i]] <- xi  ## epsilon(t)
    true.pred[i] <- f(x[i,],2) ;
    y[[i]] <- true.pred[i] * ( exp(-gamma*(t1-tt)^2) - exp(-gamma*(t2-tt)^2) + exp(-gamma*(t3-tt)^2)) + sigma * xx[[i]]
  }
  # XX=\sum_k a[[i]][k] k(.-ttt[[i]][k])
  tt <- matrix(rep(tt, n),n,nt,byrow=T)
  return(list(y=y, x=x, true.pred=true.pred,tt=tt))
}
xryf3 <- function(n, nt, p=10, pbm=100, sd,ssq){
  tt <- seq(0,1,len=nt)
  BM <- matrix(0, n, nt)
  a <- matrix(rnorm(n*pbm,sd=sd),n,pbm)
  for(j in 1: pbm) {
    a[,j] <- a[,j]/((j-1/2)*pi)
    for(i in 1:n){
      BM[i,] <- BM[i,] + a[i,j]* fpcawiener(tt,j)
    }    
  }
  x <- matrix(rnorm(n*p,sd=sd),n,p)
  tt <- (1:nt)/nt
  y <- NULL
  true.pred <- rep(0, n)  ## squared x is not good for FAVE - 
  for(i in 1:n){
#     true.pred[i] <- 1+10 * sum(x[i,]^2) 
    true.pred[i] <- sum(x[i,]^2) * x[i,1]
    y[[i]] <-  true.pred[i]  *  BM[i,] 
  }
  tt <- matrix(rep(tt, n),n,nt,byrow=T)
  return(list(y=y, x=x, true.pred=true.pred,tt=tt))
}
xryf1.nnew <- function(n, nt, p=10, pbm=100, sd,ssq){
  tt <- (1:nt)/nt
  BM <- matrix(0, n, nt)
  a <- matrix(rnorm(n*pbm,sd=sd),n,pbm)
  for(j in 1: pbm) {
    a[,j] <- a[,j]/((j-1/2)*pi)
    for(i in 1:n){
      BM[i,] <- BM[i,] + a[i,j]* fpcawiener(tt,j)
    }    
  }
  x <- matrix(rnorm(n*p,sd=sd),n,p)
  tt <- seq(0,1,len=nt)
  y <- NULL
  true.pred <- rep(0, n)  #
  rhotb <- rep(0,nt)
  for(j in 1:3){
    rhotb <-  rhotb +  fpcawiener(tt,j) 
  }
  for(i in 1:n){
   
    true.pred[i] <- x[i,1] / (1+ exp(x[i,2])) + x[i,3]
    
    y[[i]] <- true.pred[i] * (rhotb) + ssq * BM[i,]
  }
  tt <- matrix(rep(tt, n),n,nt,byrow=T)
  return(list(y=y, x=x, true.pred=true.pred,tt=tt))
}
xryf3.new <- function(n, nt, p=10, mu,sd, N=7, gamma=7,sigma=.1){  #sd is for x , sigma is for Y
  
  t1 <- 0.2; t2 <- 0.5; t3 <- 0.8#;t2<-0.6;t3<-0.9
  t11 <- 0.5;t22<-0.5;t33<-0.8
  tt <- (1:nt)/nt
  ###################################
  # Generate X
  ###################################
  
  x <- matrix(rnorm(n*p,mean=mu, sd=sd),n,p)
  true.pred <- rep(0, n)
  tt <- (1:nt)/nt
  y <- NULL
  ttt <- NULL # it is for the function x
  a <- NULL
  xx <- NULL # which is for epsilon(t)
  for(i in 1:n){
    xi <- rep(0, nt)
    a[[i]] <- rnorm(N)
    ttt[[i]] <- runif(N)
    for(k in 1:N){
      xi <- xi + a[[i]][k] * exp(-gamma*(ttt[[i]][k]-tt)^2)
    }
    xx[[i]] <- xi  ## epsilon(t)
    true.pred[i] <- f(x[i,],1) ;
    y[[i]] <- true.pred[i] *  sigma * xx[[i]] +  ( exp(-gamma*(t1-tt)^2) - exp(-gamma*(t2-tt)^2) + exp(-gamma*(t3-tt)^2))
  }
  # XX=\sum_k a[[i]][k] k(.-ttt[[i]][k])
  tt <- matrix(rep(tt, n),n,nt,byrow=T)
  return(list(y=y, x=x, true.pred=true.pred,tt=tt))
}
xryf4 <- function(n, nt, p=10, sd, ssq){  #sd is for x , sigma is for Y
  x <- matrix(rnorm(n*p,sd=sd),n,p)
  tt <- seq(0,1,len=nt)
  y <- NULL
  true.pred <- rep(0, n)
  BM <- matrix(0, n, nt)
  a <- matrix(rnorm(n*p,sd=sd),n,100)
  for(j in 1: 100) {
    a[,j] <- a[,j]/((j-1/2)*pi)
    for(i in 1:n){
      BM[i,] <- BM[i,] + a[i,j]* fpcawiener(tt,j)
    }    
  }
  for(i in 1:n){
    rhoi <- rep(0, nt)
    for(k in 1:5){
      rhoi <- rhoi + fpcawiener(tt,k)
    }
    true.pred[i] <- f(x[i,],1)
    y[[i]] <- fpcawiener(tt,0) + rhoi*true.pred[i] + ssq * BM[i,]
  }
  tt <- matrix(rep(tt, n),n,nt,byrow=T)
  return(list(y=y, x=x, true.pred=true.pred,tt=tt))
}
xryf5 <- function(n, nt, p=10, sd, ssq){
  x <- matrix(rnorm(n*p,sd=sd),n,p)
  tt <- seq(0,1,len=nt)
  y <- NULL
  true.pred <- rep(0, n)
  BM <- matrix(0, n, nt)
  a <- matrix(rnorm(n*p,sd=sd),n,100)
  for(j in 1: 100) {
    a[,j] <- a[,j]/((j-1/2)*pi)
    for(i in 1:n){
      BM[i,] <- BM[i,] + a[i,j]* fpcawiener(tt,j)
    }    
  }
  for(i in 1:n){
    rhoi <- rep(0, nt)
    for(k in 1:5){
      rhoi <- rhoi +fpcawiener(tt,k)
    }
    true.pred[i] <- x[i,1]/(1+exp(x[i,2])) + x[i,3]
    y[[i]] <-  fpcawiener(tt,0) + rhoi*true.pred[i] + ssq * BM[i,]
  }
  tt <- matrix(rep(tt, n),n,nt,byrow=T)
  return(list(y=y, x=x, true.pred=true.pred,tt=tt))
}
xryf1.extra <- function(n, nt, p=10, sd, ssq){  #sd is for x , sigma is for Y
  x <- matrix(rnorm(n*p,sd=sd),n,p)
  tt <- seq(0,1,len=nt)
  y <- NULL
  true.pred <- rep(0, n)
  BM <- matrix(0, n, nt)
  a <- matrix(rnorm(n*p,sd=sd),n,100)
  for(j in 1: 100) {
    a[,j] <- a[,j]/((j-1/2)*pi)
    for(i in 1:n){
      BM[i,] <- BM[i,] + a[i,j]* fpcawiener(tt,j)
    }    
  }
  for(i in 1:n){
    rhoi <- rep(0, nt)
    for(k in 1:5){
      rhoi <- rhoi + fpcawiener(tt,k)
    }
    true.pred[i] <- (x[i,1]^2 + x[i,2]^2)
    y[[i]] <- fpcawiener(tt,0) + rhoi*true.pred[i] + ssq * BM[i,]
  }
  tt <- matrix(rep(tt, n),n,nt,byrow=T)
  return(list(y=y, x=x, true.pred=true.pred,tt=tt))
}
################################################################
# 3.4. X:Function, y:Function----------
################################################################
# Brownian motion------
# 3.2.1 Model IV-1---------------------------
mu2 <- function(j, tt){
  exp(-(tt-j/10)^2)*tt^j
}
model.xfyf1.new <- function(n, ntx, nty, N=5, N2=20,  sigma=.1){
  # beta1(t) = k( - t1)+1k(-t2)+ 1k(-t3)
  # beta2(t) = k(  - t11)+ 1k(-t22) + 1k(-t33)
  t1 <- 0.2;t2<-0.5;t3<-0.8
  t11 <- 0.9#;t22<-0.5;t33<-0.8
  tt <- (1:nt)/nt
  ###################################
  # Generate X
  ###################################
  x <- NULL   # observed x 
  ttt <- NULL # it is for the function x
  a <- NULL
  a2 <- NULL
  xx <- NULL
  ttt2<- NULL
  for(i in 1:n){
    xi <- rep(0, nt)
    xi2 <- rep(0,nt)
    a[[i]] <- rnorm(N)
    a2[[i]] <- rnorm(N2)
    ttt[[i]] <- runif(N)
    ttt2[[i]] <- runif(N2)  # a2, ttt2 for epsilon(t)
    for(k in 1:N){
      xi <- xi + a[[i]][k] * exp(-gamma*(ttt[[i]][k]-tt)^2)
      xi2 <- xi2 + a2[[i]][k] * exp(-gamma*(ttt2[[i]][k]-tt)^2)
    }
    x[[i]] <- xi
    xx[[i]] <- xi2  # epsilon(t)
  }
  # X=\sum_k a[[i]][k] k(.-ttt[[i]][k])
  
  ###################################
  # Generate Y depending on model   #
  ###################################
  b1 <- rep(0, n) # b1[i] <- <beta_1 , X_i>
  b2 <- rep(0, n) # b2[i] <- <beta_2 , X_i>
  b1sq <- rep(0, n) # b1sq[i] <- <beta_1 , X_i^2>
  b2sq <- rep(0, n) # b2sq[i] <- <beta_2 , X_i^2>
  true.pred <- rep(0, n)
  for(i in 1:n){
    b1[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t1)^2)) #+ sum(a[[i]]*exp(-(ttt[[i]]-t2)^2)) + sum(a[[i]]*exp(-(ttt[[i]]-t3)^2))
    b2[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t11)^2)) #+ sum(a[[i]]*exp(-(ttt[[i]]-t22)^2)) + sum(a[[i]]*exp(-(ttt[[i]]-t33)^2))
    true.pred[i] <- b1[i] + b2[i]
    #b1sq[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t1)^2))^2 #+ sum(a[[i]]*exp(-(ttt[[i]]-t2)^2))^2 + sum(a[[i]]*exp(-(ttt[[i]]-t3)^2))^2
    #b2sq[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t11)^2))^2 #+ sum(a[[i]]*exp(-(ttt[[i]]-t22)^2))^2+ sum(a[[i]]*exp(-(ttt[[i]]-t33)^2))^2
    y[[i]] <-  true.pred[i] * ( exp(-gamma*(t1-tt)^2) - exp(-gamma*(t2-tt)^2) + exp(-gamma*(t3-tt)^2)) + sigma * xx[[i]]
  }
  tt = matrix(tt,n,nt,byrow=TRUE)
  result <- list(true.pred = true.pred, x=x, y=y, tt=tt)
  return(result)
}
model.xfyf1 <- function(n, ntx, nty, p, sd=1, ssq){
  result <- generate.fpca(n=n, nt=ntx, p=p, sd=sd)
  tt <- seq(0,1,len=nty)
  BM <- matrix(0, n, nty)
  a <- matrix(rnorm(n*100,sd=sd),n,100)
  for(j in 1:100) {
    a[,j] <- a[,j]/((j-1/2)*pi)
    for(i in 1:n){
      BM[i,] <- BM[i,] + a[i,j]* fpcawiener(tt,j)
    }    
  }
  nnt <- nty * 5
  tty <-seq(0,1,len=nnt)
  # unbalanced
  ty.ub <- matrix(0, n, nty)
  for(i in 1:n){
    ty.ub[i, ] <- sort(-sample(seq(0,1,len=nnt), nty))
  }
  #balanced
  ty.b <- round(seq(1, nnt, len=nty))
  ty.b <- matrix(ty.b, nrow=n, ncol=nty, byrow=T)/ nnt
  

  a <- result$a
  x.ub <- rep(0, n)
  x.b <- rep(0, n)
  b1 <- matrix(0, 2, ntx)  # row 1 : unbalnced, row 2 : balanced
  b2 <- matrix(0, 2, ntx)  # row 1 : unbalnced, row 2 : balanced
  tx.ub <- result$t.ub
  tx.b <- result$t.b
  true.pred <- rep(0,n)
  y.ub <- NULL
  y.b <- NULL
  for(i in 1:n){
    true.pred[i] <- a[i,1] + a[i,2]
    rhotb <- rep(0,nty)
    rhotub <- rep(0,nty)
    for(j in 1:5){
      rhotb <-  rhotb +  fpcawiener(ty.b[i,],j) 
      rhotub <-  rhotub + fpcawiener(ty.ub[i,],j)
    }
    y.ub[[i]] <- true.pred[i] * (rhotub) + ssq * BM[i,]
    y.b[[i]] <- true.pred[i] * (rhotb)+ ssq * BM[i,]
  }
  result$true.pred <- true.pred
  result$y.ub <- y.ub
  result$y.b <- y.b
  result$ty.ub <- ty.ub
  result$ty.b <- ty.b
  return(result)
}
model.xfyf2 <- function(n, ntx, nty, p, sd=1, ssq){
  result <- generate.fpca(n=n, nt=ntx, p=p, sd=sd)
  tt <- seq(0,1,len=nty)
  BM <- matrix(0, n, nty)
  a <- matrix(rnorm(n*100,sd=sd),n,100)
  for(j in 1:100) {
    a[,j] <- a[,j]/((j-1/2)*pi)
    for(i in 1:n){
      BM[i,] <- BM[i,] + a[i,j]* fpcawiener(tt,j)
    }    
  }
  nnt <- nty * 5
  tty <-seq(0,1,len=nnt)
  # unbalanced
  ty.ub <- matrix(0, n, nty)
  for(i in 1:n){
    ty.ub[i, ] <- sort(-sample(seq(0,1,len=nnt), nty))
  }
  #balanced
  ty.b <- round(seq(1, nnt, len=nty))
  ty.b <- matrix(ty.b, nrow=n, ncol=nty, byrow=T)/ nnt
  
  
  a <- result$a
  x.ub <- rep(0, n)
  x.b <- rep(0, n)
  b1 <- matrix(0, 2, ntx)  # row 1 : unbalnced, row 2 : balanced
  b2 <- matrix(0, 2, ntx)  # row 1 : unbalnced, row 2 : balanced
  tx.ub <- result$t.ub
  tx.b <- result$t.b
  true.pred <- rep(0,n)
  y.ub <- NULL
  y.b <- NULL
  for(i in 1:n){
    true.pred[i] <- a[i,1]/(1+exp(a[i,2])) + a[i,3]
    rhotb <- rep(0,nty)
    rhotub <- rep(0,nty)
    for(j in 1:5){
      rhotb <-  rhotb +  fpcawiener(ty.b[i,],j) 
      rhotub <-  rhotub + fpcawiener(ty.ub[i,],j)
    }
    y.ub[[i]] <- true.pred[i] * (rhotub) + ssq * BM[i,]
    y.b[[i]] <- true.pred[i] * (rhotb)+ ssq * BM[i,]
  }
  result$true.pred <- true.pred
  result$y.ub <- y.ub
  result$y.b <- y.b
  result$ty.ub <- ty.ub
  result$ty.b <- ty.b
  return(result)
}

model.xfyf2.new <- function(n, ntx, nty, N=5, N2=20,  sigma=.1){
  # beta1(t) = k( - t1)+1k(-t2)+ 1k(-t3)
  # beta2(t) = k(  - t11)+ 1k(-t22) + 1k(-t33)
  t1 <- 0.2;t2<-0.5;t3<-0.8
  t11 <- 0.9#;t22<-0.5;t33<-0.8
  tt <- (1:nt)/nt
  ###################################
  # Generate X
  ###################################
  x <- NULL   # observed x 
  ttt <- NULL # it is for the function x
  a <- NULL
  a2 <- NULL
  xx <- NULL
  ttt2<- NULL
  for(i in 1:n){
    xi <- rep(0, nt)
    xi2 <- rep(0,nt)
    a[[i]] <- rnorm(N)
    a2[[i]] <- rnorm(N2)
    ttt[[i]] <- runif(N)
    ttt2[[i]] <- runif(N2)  # a2, ttt2 for epsilon(t)
    for(k in 1:N){
      xi <- xi + a[[i]][k] * exp(-gamma*(ttt[[i]][k]-tt)^2)
      xi2 <- xi2 + a2[[i]][k] * exp(-gamma*(ttt2[[i]][k]-tt)^2)
    }
    x[[i]] <- xi
    xx[[i]] <- xi2  # epsilon(t)
  }
  # X=\sum_k a[[i]][k] k(.-ttt[[i]][k])
  
  ###################################
  # Generate Y depending on model   #
  ###################################
  b1 <- rep(0, n) # b1[i] <- <beta_1 , X_i>
  b2 <- rep(0, n) # b2[i] <- <beta_2 , X_i>
  b1sq <- rep(0, n) # b1sq[i] <- <beta_1 , X_i^2>
  b2sq <- rep(0, n) # b2sq[i] <- <beta_2 , X_i^2>
  true.pred <- rep(0, n)
  for(i in 1:n){
    b1[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t1)^2)) #+ sum(a[[i]]*exp(-(ttt[[i]]-t2)^2)) + sum(a[[i]]*exp(-(ttt[[i]]-t3)^2))
    b2[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t11)^2)) #+ sum(a[[i]]*exp(-(ttt[[i]]-t22)^2)) + sum(a[[i]]*exp(-(ttt[[i]]-t33)^2))
    true.pred[i] <- log(b1[i]^2 + b2[i]^2)
    #b1sq[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t1)^2))^2 #+ sum(a[[i]]*exp(-(ttt[[i]]-t2)^2))^2 + sum(a[[i]]*exp(-(ttt[[i]]-t3)^2))^2
    #b2sq[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t11)^2))^2 #+ sum(a[[i]]*exp(-(ttt[[i]]-t22)^2))^2+ sum(a[[i]]*exp(-(ttt[[i]]-t33)^2))^2
    y[[i]] <-  true.pred[i] * ( exp(-gamma*(t1-tt)^2) - exp(-gamma*(t2-tt)^2) + exp(-gamma*(t3-tt)^2)) + sigma * xx[[i]]
  }
  tt = matrix(tt,n,nt,byrow=TRUE)
  result <- list(true.pred = true.pred, x=x, y=y, tt=tt)
  return(result)
}
model.xfyf3 <- function(n, ntx, nty, p, sd=1, ssq){
  result <- generate.fpca(n=n, nt=ntx, p=p, sd=sd)
  tt <- seq(0,1,len=nty)
  BM <- matrix(0, n, nty)
  a <- matrix(rnorm(n*100,sd=sd),n,100)
  for(j in 1:100) {
    a[,j] <- a[,j]/((j-1/2)*pi)
    for(i in 1:n){
      BM[i,] <- BM[i,] + a[i,j]* fpcawiener(tt,j)
    }    
  }
  nnt <- nty * 5
  tty <-seq(0,1,len=nnt)
  # unbalanced
  ty.ub <- matrix(0, n, nty)
  for(i in 1:n){
    ty.ub[i, ] <- sort(-sample(seq(0,1,len=nnt), nty))
  }
  #balanced
  ty.b <- round(seq(1, nnt, len=nty))
  ty.b <- matrix(ty.b, nrow=n, ncol=nty, byrow=T)/ nnt
  
  
  a <- result$a
  x.ub <- rep(0, n)
  x.b <- rep(0, n)
  b1 <- matrix(0, 2, ntx)  # row 1 : unbalnced, row 2 : balanced
  b2 <- matrix(0, 2, ntx)  # row 1 : unbalnced, row 2 : balanced
  tx.ub <- result$t.ub
  tx.b <- result$t.b
  true.pred <- rep(0,n)
  y.ub <- NULL
  y.b <- NULL
  for(i in 1:n){
    true.pred[i] <- a[i,1] + a[i,2]
    rhotb <- rep(0,nty)
    rhotub <- rep(0,nty)
    
    y.ub[[i]] <- true.pred[i] * sin(true.pred[i] * tt * pi)  * BM[i,]
    y.b[[i]] <- true.pred[i] * sin(true.pred[i] * tt * pi) * BM[i,]
  }
  result$true.pred <- true.pred
  result$y.ub <- y.ub
  result$y.b <- y.b
  result$ty.ub <- ty.ub
  result$ty.b <- ty.b
  return(result)
}
model.xfyf3.new <- function(n, ntx, nty, N=5, N2=20,  sigma=.1){
  # beta1(t) = k( - t1)+1k(-t2)+ 1k(-t3)
  # beta2(t) = k(  - t11)+ 1k(-t22) + 1k(-t33)
  t1 <- 0.2;t2<-0.5;t3<-0.8
  t11 <- 0.9#;t22<-0.5;t33<-0.8
  tt <- (1:nt)/nt
  ###################################
  # Generate X
  ###################################
  x <- NULL   # observed x 
  ttt <- NULL # it is for the function x
  a <- NULL
  a2 <- NULL
  xx <- NULL
  ttt2<- NULL
  for(i in 1:n){
    xi <- rep(0, nt)
    xi2 <- rep(0,nt)
    a[[i]] <- rnorm(N)
    a2[[i]] <- rnorm(N2)
    ttt[[i]] <- runif(N)
    ttt2[[i]] <- runif(N2)  # a2, ttt2 for epsilon(t)
    for(k in 1:N){
      xi <- xi + a[[i]][k] * exp(-gamma*(ttt[[i]][k]-tt)^2)
      xi2 <- xi2 + a2[[i]][k] * exp(-gamma*(ttt2[[i]][k]-tt)^2)
    }
    x[[i]] <- xi
    xx[[i]] <- xi2  # epsilon(t)
  }
  # X=\sum_k a[[i]][k] k(.-ttt[[i]][k])
  
  ###################################
  # Generate Y depending on model   #
  ###################################
  b1 <- rep(0, n) # b1[i] <- <beta_1 , X_i>
  b2 <- rep(0, n) # b2[i] <- <beta_2 , X_i>
  b1sq <- rep(0, n) # b1sq[i] <- <beta_1 , X_i^2>
  b2sq <- rep(0, n) # b2sq[i] <- <beta_2 , X_i^2>
  true.pred <- rep(0, n)
  for(i in 1:n){
    b1[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t1)^2)) #+ sum(a[[i]]*exp(-(ttt[[i]]-t2)^2)) + sum(a[[i]]*exp(-(ttt[[i]]-t3)^2))
    b2[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t11)^2)) #+ sum(a[[i]]*exp(-(ttt[[i]]-t22)^2)) + sum(a[[i]]*exp(-(ttt[[i]]-t33)^2))
    true.pred[i] <- b1[i] + b2[i]
    #b1sq[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t1)^2))^2 #+ sum(a[[i]]*exp(-(ttt[[i]]-t2)^2))^2 + sum(a[[i]]*exp(-(ttt[[i]]-t3)^2))^2
    #b2sq[i] <- sum(a[[i]]*exp(-gamma*(ttt[[i]]-t11)^2))^2 #+ sum(a[[i]]*exp(-(ttt[[i]]-t22)^2))^2+ sum(a[[i]]*exp(-(ttt[[i]]-t33)^2))^2
    y[[i]] <-  x[[i]] + true.pred[i] *  sigma * xx[[i]]
  }
  tt = matrix(tt,n,nt,byrow=TRUE)
  result <- list(true.pred = true.pred, x=x, y=y, tt=tt)
  return(result)
}


xfyf1 <- function(n, nt, p=10, sd){
  x <- matrix(rnorm(n*p,sd=sd),n,p)
  tt <- (1:nt)/nt
  y <- NULL
  true.pred <- rep(0, n)
  for(i in 1:n){
    rhoi <- rep(0, nt)
    for(k in 1:4){
      rhoi <- rhoi + mu_gauss(k,tt)
    }
    true.pred[i] <- f(x[i,],1)
    y[[i]] <- mu_gauss(0,tt) + rhoi*true.pred[i]
  }
  tt <- matrix(rep(tt, n),n,nt,byrow=T)
  return(list(y=y, x=x, true.pred=true.pred,tt=tt))
}
