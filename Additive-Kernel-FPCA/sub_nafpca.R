#--------------------------------------------------------------#
# get.ftn : get data matrix -> ftn coefficient, inner products..
#--------------------------------------------------------------#
## Input
# balanced case :  X = n X p   or X = n X nt   // tt = nt-dim'l vector
# unbalanced case : X : list  x[[1]] - x[[n]],   tt: list tt[[1]] - tt[[n]]

get.ftn = function(x, tt=NULL, type="ftn", unbalanced=FALSE, p=1,
                   basisname = "bspline", nbasis=11,
                   lambdascale=2, ncv=10){
  require(fda)
  if(type=="multivariate"){  # x : n x p matrix
    x <- as.matrix(x)
    #if(center==TRUE) x <- center(x)
    n <- nrow(as.matrix(x))
    normX <- normmat(x)
    basis=NULL;xcoef=NULL  ## those are not needed in vectors
    ntx=NULL; xip <- tcrossprod(x)
    basis <- NULL ; xcoef=NULL; xcoef.org=NULL;
    binnerprod=NULL;lambdaout=NULL

    out = list(n=n,coef=xcoef, lambdaout=lambdaout, xip=xip, normx=normX, binnerprod=binnerprod, basis=basis)
  }else if(type=="ftn"){
    if(nbasis%%2==0) nbasis <- nbasis - 1 ## odd number of basis
    lambda_all<- exp(seq(log(10^(-10)),log(10^(-lambdascale)),length=ncv))

    if(!is.list(tt)){
      # x: n X nt matrix  x[1,] is the first curve
      n = nrow(x)
      tmin = min(tt) ; tmax<-max(tt)
      basis <- eval(as.name(paste("create.",basisname,".basis",sep='')))(c(tmin,tmax), nbasis=nbasis) # maxnbasis
      binnerprod <- eval(as.name(paste(basisname,"pen",sep='')))(basis,0)  # inner products of basis
      RR<-eval(as.name(paste(basisname,"pen",sep='')))(basis,2)  # inner prodcut of second derivatives .

      temps <-approx_mat(t(x),tt,basis,RR,lambdagrid=lambda_all)
      xcoef <- temps$mucoef; lambdaout <- temps$lambdaout
      xip <- t(xcoef) %*% binnerprod %*% xcoef
      normX = t(-2 * xip + diag(xip)) + diag(xip)

      out = list(n=n,coef=xcoef, lambdaout=lambdaout, xip=xip, normx=normX, binnerprod=binnerprod, basis=basis, p=p)
    }else if(is.list(tt)){  ## Still on the same domain
      # x: n-list.. x[[1]] is the first curve
      tunique = unique(unlist(tt))
      tmin = min(tunique); tmax = max(tunique)
      basis <- eval(as.name(paste("create.",basisname,".basis",sep='')))(c(tmin,tmax), nbasis=nbasis) # maxnbasis
      binnerprod <- eval(as.name(paste(basisname,"pen",sep='')))(basis,0)  # inner products of basis
      RR<-eval(as.name(paste(basisname,"pen",sep='')))(basis,2)  # inner prodcut of second derivatives .

      if(p==1){

        n = length(x)
        temps <-approx_list(x,tt,basis,RR,lambdagrid=lambda_all)
        xcoef <- temps$mucoef; lambdaout <- temps$lambdaout
        xip <- t(xcoef) %*% binnerprod %*% xcoef
        normX = t(-2 * xip + diag(xip)) + diag(xip)

        out = list(n=n,coef=xcoef, lambdaout=lambdaout, xip=xip, normx=normX, binnerprod=binnerprod, basis=basis, p=p)
      }else{
        # x$x1 , ... , x$xp   :
        # x$xi : list of n functions, i = 1. ... p
        # tt$x1, ..., tt$xp
        n = length(x[[1]])
        xcoef = array(0,c(nbasis,n,p))
        xip = array(0, c(n, n, p))
        normXs = array(0, c(n,n,p))
        normX = matrix(0,n,n)
        for(j in 1:p){
          temps <-approx_list(x[[j]],tt[[j]],basis,RR,lambdagrid=lambda_all)

          xcoef[,,j] <- temps$mucoef; lambdaout <- temps$lambdaout  #lambdaout is not really needed

          xip[,,j] <- t(xcoef[,,j]) %*% binnerprod %*% xcoef[,,j]
          normXs[,,j] = t(-2 * xip[,,j] + diag(xip[,,j])) + diag(xip[,,j])
          normX = normX + normXs[,,j]
        }

        out = list(n=n,coef=xcoef, lambdaout=lambdaout, xip=xip, normx=normX, normxs=normXs,binnerprod=binnerprod, basis=basis, p=p)
      }
    }
  }
  return(out)
}


plotpred <- function(x, n.tuning, vec, xlab="x", ylab="y", main="", ind.info = NULL){
  n <- nrow(x)
  if(!is.null(ind.info)){
    ind.info <- as.numeric(factor(ind.info))  # change indicator number
    nclass <- length(unique(ind.info))
    plot(x[which(ind.info==1),vec[1]], x[which(ind.info==1),vec[2]],
         xlim=c(min(x[, vec[1]]),max(x[, vec[1]])),
         ylim=c(min(x[, vec[2]]),max(x[, vec[2]])),
         xlab=xlab, ylab=ylab, main=main, pch=16)
    for(i in 2:(nclass)){
      points(x[which(ind.info==i),vec[1]], x[which(ind.info==i),vec[2]], col=i,pch=16)
    }
  }else {
    if(n %% n.tuning != 0 ) stop('n is not multiple of n.tuning')
    nclass <- n/n.tuning
    plot(x[1:n.tuning, vec[1]],x[1:n.tuning, vec[2]],
         xlim=c(min(x[, vec[1]]),max(x[, vec[1]])),
         ylim=c(min(x[, vec[2]]),max(x[, vec[2]])),
         xlab=xlab, ylab=ylab, main=main,pch=paste(1))
    for(i in 1:(nclass-1)){
      points(x[(i*n.tuning+1):((i+1)*n.tuning), vec[1]],x[(i*n.tuning+1):((i+1)*n.tuning), vec[2]], col=i+1, pch=paste(i+1))
    }
  }


}

######################################
# misrate : mis classification rate of total
######################################
misrate <- function(pred, class){
  n=length(class)
  return(sum(pred!=class)/n)
}
eval.ftn = function(the_ftn, tt){
  coef=the_ftn$coef
  basis=the_ftn$basis
  return(t(eval.basis(tt, basis)%*% coef))
}



######################################
# My data <-> funData / MFPCA package
######################################
# my data -> funData
# MFPCA does not work properly for unbalanced data!@!!@!@!@@!
convert.fd = function(fd, tt=NULL, p=1, n=100){
  class.fd = class(fd)

  if(class.fd=="funData"){
    fd.type = 1
    out = convert.fd.t1(fd)
  }else if(class.fd=="multiFunData"){
    fd.type = 2
    if(p>2) stop('need coding')
    out1 = convert.fd.t1(fd[[1]])
    out2 = convert.fd.t1(fd[[2]])
    xx = list(out1$x,out2$x)
    tx = list(out1$tt, out2$tt)
    out = list(x=xx,tt=tx)

  }else if(class.fd=="list"){
    fd.type = 3  # my data
    if(p==1){
      ttt = sort(unique(unlist(tt)))
      nt = length(ttt)
      o1 = matrix(NA, n, nt)
      for(i in 1:n){
        tmp = rep(NA, nt)
        tmp[ttt %in% tt[[i]] ] = x[[i]]
        o1[i,] = tmp
      }
      out = funData(ttt, o1)
    }
    #only when p=2


    if(p==2 & is.list(tt)){
      ttt = sort(unique(unlist(tt)))
      nt = length(ttt)
      o1 = matrix(NA, n, nt)
      o2 = matrix(NA, n, nt)
      for(i in 1:n){
        tmp1 = rep(NA, nt)
        tmp2 = rep(NA, nt)
        tmp1[ttt %in% tt[[1]][[i]] ] = x[[1]][[i]]
        tmp2[ttt %in% tt[[2]][[i]] ] = x[[2]][[i]]
        o1[i,] = tmp1
        o2[i,] = tmp2
      }
      fd1 = funData(ttt, o1)
      fd2 = funData(ttt, o2)
      out = multiFunData(list(fd1,fd2))
    }
    if(p==2 & !is.list(tt)){
      M = length(tt)
      o1 = matrix(NA, n, M)
      o2 = matrix(NA, n, M)
      for(i in 1:n){
        o1[i,] = x[[1]][[i]]
        o2[i,] = x[[2]][[i]]
        #fi=multiFunData(list(funData(tt[[1]][[i]], matrix(x[[1]][[i]],nrow=1)),funData(tt[[2]][[i]], matrix(x[[2]][[i]],nrow=1))))
      }
      fd1 = funData(tt, o1)
      fd2 = funData(tt, o2)
      out = multiFunData(list(fd1,fd2))
    }
    if(p>2) stop('need coding')
  }else if(class.fd=="matrix"){
    fd.type = 4 # my data
    out = funData(tt, fd)
  }
  return(out)
}

convert.fd.t1 = function(fd){
  tt = (fd@argvals)[[1]]
  x = fd@X
  n = dim(x)[1];
  # check the sparsity of the functional data
  if(length(which(is.na(fd@X), arr.ind=TRUE))==0){
    sparse=FALSE
    tx = tt
    xx = x
  }else{
    sparse = TRUE
    sparse.ind = which(is.na(x), arr.ind=TRUE)
    id.sparse = unique(sparse.ind[,1])
    xx = NULL
    tx = NULL
    for(i in (1:n)){
      xx[[i]] = as.vector(na.omit(x[i,]))
      tx[[i]] = tt[!is.na(x[i,])]
    }
  }
  out = list(tt=tx, x=xx)
  return(out)
}

#-----------------------------------------------------#
# Formula (based on the number of predictor, d)
#-----------------------------------------------------#
class.formula = function(d){
  fx <- "X1"
  for(i in 1:(d-1)) fx <- paste(fx,"+X",i+1,sep="")
  return(as.formula(paste("y~",fx)))
}
#-----------------------------------------------------#
#-----------------------------------------------------#
# Classify
# possible methods : c("lda", "qda", "svm", "ctree", "randomForest", "naiveBayes")
#-----------------------------------------------------#
classify = function(my.formula, train.data, test.data=NULL, method){
  require(MASS)  # lda, qda
  #library(party) # Decision Tree Conditional Inference Tree
  #require(randomForest)  # as it is
  library(e1071) # SVM
  trained = eval(as.symbol(method))(formula=my.formula, data=train.data)
  if(is.null(test.data)){
    if(method=="lda"|method=="qda") return(predict(trained, train.data)$class)
    if(method=="naiveBayes") return(predict(trained, train.data))
    return(predict(trained))
  }else{
    if(method=="lda"|method=="qda") return(predict(trained, test.data)$class)
    return(predict(trained, test.data))
  }
}

class.table = function(pred, true){
  n = length(unique(true))
  pred = factor(pred); true=factor(true)
  clss = sort(levels(true))
  out = matrix(0,n+2,n+2)

  for(i in 1:n){
    ind = which(true==true[i])
    out[i,1:n] = as.vector(table(true[pred==clss[i]]))
  }
  out[1:n,n+1] = apply(out[1:n,1:n], 1,sum)
  out[n+1,1:n] = apply(out[1:n,1:n], 2,sum)
  for(i in 1:n){
    out[i,n+2] = out[i,i]/out[i,n+1]*100
    out[n+2,i] = out[i,i]/out[n+1,i]*100
  }
  out[n+2,n+2] = sum(diag(out[1:n,1:n]))/sum(out[n+1,1:n])*100
  return(out)
}
