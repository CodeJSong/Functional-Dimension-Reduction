# Set working directory to the folder, functional dimension reduction 
setwd("~/CloudStation/Codes/GitHub/Functional-Dimension-Reduction")
source('sub_codes.R')
Rcpp::sourceCpp('cpp_code/nonlinear_subcodes.cpp')

##### For BIC-type dimension selection
#fun1 <- function(k) 1/log(k+1)
fun1 <- function(n, c=1) c * n^{-1/4}*log(n)
fun2 <- function(k) k
fdr.bic <- function(eigval, fun1=fun1, fun2=fun2, c=1){
  n <- length(eigval)
  out <- rep(0,n)
  for(k in 1:n){
    out[k] <- sum(eigval[1:k]) - eigval[1] * fun1(n, c) * fun2(k)
  }
  return(out)
}

################################
# Part 2. main code : fdr ------
################################
##################################################################
# fdr(x, y, xftn=FALSE, tequal=FALSE, tt, d=1)------
##################################################################
# INPUT VARIABLES
# xftn : indicator whether x is function or not
# tequal : indicator whether t is equally spaced or not
# tt : list of time points (since t() is a transpose function,
#                          named it tt)
#   example : tt[[1]], ..., tt[[n]], 
#           for each tt[[i]] is a vector of time points
# x : n x p matrix (when x is not a function)  or
#    a list of function, or dataframe
# y :  n x p matrix (when y is not a function)  or
#    a list of function, or dataframe
##################################################################
# OUTPUT
#  list(Mat = Q %*% GX.inv.half %*% A,  - it is for prediction
#     if g is our estimated function for central class
#      then g(x.new) = KX.new %*% Mat[,1] 
#         where KX.new=(k(x.new,x1), ..., k(x.new, xn))
#
#       t.unique = t.unique,  - unique time points
#       gt = gt, KT = KT,     - width parameter, gram matrix for T
#       gx = gx, KX = KX,     - width parameter, gram matrix for X
#       gy=gy, KY=KY,         - width parameter, gram matrix for Y
#       pred = pred[,1])      - n-vector (g(x1), ..., g(xn))
################################################################## 

fdr<- function(x=NULL, y=NULL, tx=NULL, ty=NULL, xdata=NULL,ydaya=NULL,xftn=FALSE, yftn=FALSE, xtequal=FALSE, ytequal=FALSE,
               option="fir", kernel="brownian", yclass = FALSE, c_brown=1, cvt=FALSE, cvx=FALSE,
               dimx=1,dimy=1, ncv=20,
               sht=2, shx=2, shy=2,ex=0.3, et=0.3, ey=0.3, ety=0.3, shty=2){
  ## Make KX depending on whether X is functional or not
  if(dimx>1){
    xdata <- x
    id <- unique(xdata$id)
    xt.unique <- sort(unique(xdata$tt))
    n <- length(id)
    x <- NULL
    tx <- NULL
    for(i in 1:n){
      temp <- xdata[xdata$id==id[i],]
      tx[[i]] <- temp$tt
      x[[i]] <- matrix(as.numeric(as.matrix(temp[,3:(2+dimx)])),ncol=dimx)
    } 
  }
  if(dimy>1){
    ydata <- y
    y <- NULL
    ty <- NULL
    id <- unique(ydata$id)
    yt.unique <- sort(unique(ydata$tt))
    for(i in 1:n){
      temp <- ydata[ydata$id==id[i],]
      ty[[i]] <- temp$tt
      y[[i]] <- matrix(as.numeric(as.matrix(temp[,3:(2+dimy)])),ncol=dimy)
    }  
  }
  
  if(!xftn){  # x : n x p matrix
    x <- as.matrix(x)
    n <- nrow(as.matrix(x))
    normX <- normmat(x)
    gx <- gam(normX)/shx 
    KX <- exp(-gx*normX)  
  }
  if(!yftn && !yclass){  # y : n x p matrix
    y <- as.matrix(y)
    n <- nrow(as.matrix(y))
    normY <- normmat(y)
    gy <- gam(normY)/shy
    KY <- exp(-gy*normY)
  }
  if(yclass){ # For categorical Y // 
    n <- length(y)
    yid <- unique(y)
    nclass <- length(yid)
    KY <- matrix(0, n, n)
    for(i in 1:nclass){
      KY[y==yid[i], y==yid[i]] <- 1
    }
  }
  if(xftn){   # x is list of functions
    # input variable is little different.
    # x : list of xi(tij)
    # y : n-dimension
    # tt : list of tij
    n <- length(tx)
    xt.unique <- as.matrix(sort(unique(unlist(tx))))
    ntx <- length(xt.unique)
    # Calculate KT to specify HT
    if(kernel=="gaussian"){
      normmat.t <- normmat(xt.unique)
      gt <- gam(normmat.t)/sht
      KTX <- exp(-gt*normmat.t)  # gram matrix for T  
    }
    if(kernel=="brownian"){
      nt <- length(xt.unique)
      KTX <- matrix(0, nt, nt)
      for(i in 1:nt){
        KTX[i,] <- c(xt.unique[1:i],rep(xt.unique[i],(nt-i)))*c_brown
      }
    }
    
    # Calculate KX for function on HT
    #xip <- x.ip2(x=x, tt=tx, xdata=xdata, t.unique=xt.unique, KT=KTX, tequal=xtequal, et=et, dim=dimx)  ## XX_ij = <Xi, Xj>_{HT}
    xip <- xiplist(n, x, tx, xt.unique, KTX, et, as.numeric(xtequal), dimx) 
    
    normX = t(-2 * xip + diag(xip)) + diag(xip)
    gx <- gam(normX)/shx
    KX <- exp(-gx*normX)    
  }
  if(yftn){   # y is list of functions
    n <- length(ty)
    yt.unique <- as.matrix(sort(unique(unlist(ty))))
    nty = length(yt.unique)
    if(kernel=="gaussian"){
      normmat.ty <- normmat(yt.unique)
      gt <- gam(normmat.ty)/shty
      KTY <- exp(-gt*normmat.ty)  # gram matrix for T  
    }
    if(kernel=="brownian"){
      nt <- length(yt.unique)
      KTY<- matrix(0, nty, nty)
      for(i in 1:nty){
        KTY[i,] <- c(yt.unique[1:i],rep(yt.unique[i],(nty-i)))*c_brown
      }
    }
    #yip <- x.ip2(x=y, tt=ty, xdata=ydata, t.unique=yt.unique, KT=KTY, tequal=ytequal, et=et, dim=dimy)  ## XX_ij = <Xi, Xj>_{HT}
    yip <- xiplist(n, y, ty, yt.unique, KTY, ety, as.numeric(ytequal), dimy) 
    normY =  t(-2 * yip + diag(yip)) + diag(yip)
    gy <- gam(normY)/shy
    KY <- exp(-gy*normY)    
    #ey <- ey * n /2 
  }
  ## CV value
  
  if(cvt){
    cvwint<- c(2,exp(seq(log(0.05),log(20),len=ncv-1)))
    
    if(xftn){
      cvet <-  c(ntx^(-1/4),exp(seq(log(ntx^(-1/4)/50),log(ntx^(-1/4)*50),len=ncv-1)))
      if(kernel=="brownian"){
        cvtout <- cvtsearch(length(x), as.numeric(kernel=="brownian"),
                            cvwint,cvet, KTX, x, tx, xt.unique,
                            as.numeric(xtequal), dimx)
        et = cvtout$et
      }
      
      if(kernel=="gaussian"){
        cvtout <- cvtsearch(length(x), as.numeric(kernel=="gaussian"),
                            cvwint,cvet, KTX, x, tx, xt.unique,
                            as.numeric(xtequal), dimx)
        sht = cvtout$sht 
        et = cvtout$et
        gt <- gam(normmat.t)/sht
        KTX <- exp(-gt*normmat.t)  # gram matrix for T    
      }
      
      xip <- xiplist(n, x, tx, xt.unique, KTX, et, as.numeric(xtequal), dimx) 
      normX = t(-2 * xip + diag(xip)) + diag(xip)
      gx <- gam(normX)/shx
      KX <- exp(-gx*normX)      
    }
    if(yftn){
      cvet <-  c(nty^(-1/4),exp(seq(log(nty^(-1/4)/50),log(nty^(-1/4)*50),len=ncv-1)))
      if(kernel=="brownian"){
        cvtout <- cvtsearch(length(y), as.numeric(kernel=="brownian"),
                            cvwint,cvet, KTY, y, ty, yt.unique,
                            as.numeric(ytequal), dimy)
        ety = cvtout$et
      }
      if(kernel=="gaussian"){
        cvtout <- cvtsearch(length(y), as.numeric(kernel=="gaussian"),
                            cvwint,cvet, KTY, y, ty, yt.unique,
                            as.numeric(ytequal), dimy)
        shty = cvtout$sht 
        ety = cvtout$et
        gt <- gam(normmat.ty)/shty
        KTY <- exp(-gt*normmat.ty)  # gram matrix for T    
      }
      yip <- xiplist(n, y, ty, yt.unique, KTY, ety, as.numeric(ytequal), dimy) 
      normY =  t(-2 * yip + diag(yip)) + diag(yip)
      gy <- gam(normY)/shy
      KY <- exp(-gy*normY)    
    }
  }
  if(cvx){
    cvwinx <- c(2,exp(seq(log(0.05),log(20),len=ncv-1)))
    cvwiny <- c(2,exp(seq(log(0.05),log(20),len=ncv-1)))
    cvex <-  c(n^(-1/4),exp(seq(log(n^(-1/4)/50),log(n^(-1/4)*50),len=ncv-1)))
    if(yclass) normY=KY
    cvxout <- cvxsearch(cvwinx, cvwiny, cvex, normX, normY, as.numeric(yclass), KY) 
    #cvyout <- cvxsearch(cvwinx, cvwiny, cvex, normY, normX)
    #ey <- cvyout$ex; shyy <- cvyout$shy
    shx = cvxout$shx; shy=cvxout$shy; ex=cvxout$ex; ey = cvxout$ey;
    if(!yclass){
      gy <- gam(normY)/shy
      KY <- exp(-gy*normY)      
    }
    gx <- gam(normX)/shx
    KX <- exp(-gx*normX)
  }
  
  if(ey==.3) ey=n^(-1/4)
  
  ## Estimate central class 
  M <- matrix(0, n, n)
  one <- matrix(1, n, 1)
  Q <- diag(n) - one %*% t(one) / n    
  GX <- Q %*% KX %*% Q    # centered gram matrix
  GY <- Q %*% KY %*% Q    # centered gram matrix
  if(option=="fir"){
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
    pred <- KX %*% Mat
  }
  if(option=="FAVEold"){
    Vfave <- matrix(0, n, n)
    GY.inv <- matpower2(GY,-1,ey)
    GX.inv <- matpower2(GX,-1,ex)
    GX.inv.half <- matpower2(GX,-1/2,ex)
    KX.inv <- matpower(KX,-1)
    KY.inv <- matpower(KY,-1)
    EYX <- GY.inv %*% GX - KY.inv %*% one %*% t(one) %*% (KY/n) %*% GY.inv %*% GX + KY.inv %*% one %*% t(one) %*% KX /n
    A <- Q %*% GY.inv 
    for(i in 1:n){ #Y_i
      VV <-  A %*% (KY[,i] - KY %*% one/n) + one/n
      Vny <- diag(as.vector(VV)) - VV %*% t(VV)
      Vfave <- Vfave + t(Q/n-KX.inv %*% KX %*% Vny)%*% KX %*% (Q/n-KX.inv %*% KX %*% Vny)
    }
    Vfave <- faveloop(GX,GY,KY,KX,Q,ex,ey)
    KX.half <- matpower(KX, 1/2)
    KX.inv.half <- matpower2(KX, -1/2, ex)
    M <- KX.half %*% Vfave %*% KX.half
    A <- eigen(M, sym=TRUE)
    eigval <- A$values
    A <- A$vec
    Mat <- Q%*% KX.inv.half %*% A 
    pred <- KX%*% Mat
  }
  if(option=="FAVE"){   
    GX.inv.half <- matpower2(GX,-1/2,ex)
    if(yclass) ey=0
    M <- faveloop(GX,GY,KY,KX,Q,ex,ey, as.numeric(yclass))
    A <- eigen(sym(M), sym=TRUE)
    eigval <- A$values
    A <- A$vec
    Mat <- GX.inv.half %*% A
    pred <- KX %*% Mat
    
  }
  out <- list(Mat = Mat, eval = eigval,
              M = M, ex=ex,
              gx = gx, KX = KX,
              KY = KY,
              pred = pred)
  if(kernel=="gaussian") out$gt=gt
  if(xftn) out$KTX <- KTX
  if(yftn) out$KTY <- KTY
  if(cvt) {
    out$et= et
    out$sht=sht
    out$shty=shty
    out$ety=ety
  }
  if(cvx){
    out$shx=shx
    out$shy=shy
    out$ex=ex
    out$ey=ey
  }
  return(out)
}


center <- function(x){
  return(t(t(x)-apply(x,2,mean)))}

##################################################################
# fdr.predict, given f_hat from fdr methods, with training set, estimate f_hat(x.new)
##################################################################

# prediction ------
fdr.predict <- function(x.train, x.test, tt.train,tt.test,Mat, xftn=TRUE, yftn=TRUE, xtequal=FALSE, ytequal=FALSE, 
                        dimx=1, dimy=1, sht=2, shx=2, shy=2, gx=1, gt=1, c_brown=1, et=2,
                        option="fir", kernel="brownian"){
  if(dimx==1){  
    ## Make KX depending on whether X is functional or not
    if(!xftn){  # x : n x p matrix
      n.train <- nrow(as.matrix(x.train))
      n.test <- nrow(as.matrix(x.test))
      n <- n.train + n.test
      x <- rbind(x.train,x.test)
      normmat.x <- normmat(x)
      KX <- exp(-gx*normmat.x)  
      normX <- normmat.x
      #ex <- ex * n
    }
    if(xftn){   # x is list of functions
      # input variable is little different.
      # x : list of xi(tij)
      # y : n-dimension
      # tt : list of tij
      n.train <- length(x.train)
      n.test <- length(x.test)
      n <- n.train + n.test
      t.unique <- unique(sort(c(unlist(tt.train),unlist(tt.test))))
      # Calculate KT to specify HT
      if(kernel=="gaussian"){
        normmat.t <- normmat(as.matrix(t.unique))
        KTX <- exp(-gt*normmat.t)  # gram matrix for T  
      }
      if(kernel=="brownian"){
        nt <- length(t.unique)
        KTX <- matrix(0, nt, nt)
        for(i in 1:nt){
          KTX[i,] <- c(t.unique[1:i],rep(t.unique[i],(nt-i)))*c_brown
        }
      }
      x <- NULL # aggregate train and test to calculate inner product
      tx <- NULL
      for(i in 1:n.train){
        
        x[[i]] <- x.train[[i]]
        tx[[i]] <- tt.train[[i]]
      }
      for(i in (n.train+1):n){
        x[[i]] <- x.test[[i-n.train]]
        tx[[i]] <- tt.test[[i-n.train]]
      }
      # Calculate KX for function on HT
      xip <- xiplist(n, x, tx, t.unique, KTX, et, as.numeric(xtequal), dimx) 
      #XXP <- x.ip(x=x, tt=tx, t.unique=t.unique, KT=KTX, et=et, tequal=xtequal)  ## XX_ij = <Xi, Xj>_{HT}
      normX = t(-2 * xip + diag(xip)) + diag(xip)
    }
  }
  if(dimx>1){
    x.total <- rbind(x.train,x.test)
    id.train <- unique(x.train$id)
    id.test <- unique(x.test$id)
    id <- c(id.train,id.test)
    n.train <- length(id.train)
    n.test <- length(id.test)
    xt.unique <- sort(unique(x.total$tt))
    n <- n.train+n.test
    x <- NULL
    tx <- NULL
    for(i in 1:n){
      temp <- x.total[x.total$id==id[i],]
      tx[[i]] <- as.matrix(temp$tt)
      x[[i]] <- matrix(as.numeric(as.matrix(temp[,3:(2+dimx)])),ncol=dimx)      
    }          
    if(kernel=="gaussian"){
      normmat.t <- normmat(as.matrix(xt.unique))
      KTX <- exp(-gt*normmat.t)  # gram matrix for T  
    }
    if(kernel=="brownian"){
      nt <- length(xt.unique)
      KTX <- matrix(0, nt, nt)
      for(i in 1:nt){
        KTX[i,] <- c(xt.unique[1:i],rep(xt.unique[i],(nt-i)))*c_brown
      }
    }
    xip <- xiplist(n, x, tx, xt.unique, KTX, et, as.numeric(xtequal), dimx) 
    #xip <- x.ip2(x=x, tt=tx, t.unique=xt.unique, KT=KTX, et=et, dim=dim)  ## XX_ij = <Xi, Xj>_{HT}
    normX = t(-2 * xip + diag(xip)) + diag(xip)
  }
  
  KX.new <- exp(-gx*normX)[(n.train+1):n,1:n.train]
  predicted <- KX.new %*% Mat
  return(predicted)
}


##################################################################
# dim.search2 : current veresion dimension selection
##################################################################
dim.search2 <- function(nfdr_pred, ncvl=10){
  cvlambda <- exp(seq(log(0.01),log(1),len=ncvl))
  cvwinx <- c(2,exp(seq(log(0.05),log(20),len=ncvl-1)))
  n = nrow(nfdr_pred$pred)
  cvex <-  c(n^(-1/4),exp(seq(log(n^(-1/4)/50),log(n^(-1/4)*50),len=ncvl-1)))
  temp_dim <- rep(0, ncvl);cvdim<- rep(0, ncvl);
  for(k in 1:ncvl){
    temp_dim[k] <- which.max(fdr.bic(nfdr_pred$eval, fun1=fun1,fun2=fun2, c=cvlambda[k]))
  }
  temp_dim2 <- sort(unique(temp_dim))
  cvdim <- rep(0, length(temp_dim2))
  for(k in 1:length(temp_dim2)){
    cvdim[k] <- cvd2search(as.matrix(nfdr_pred$pred[,1:temp_dim2[k]]),nfdr_pred$KY, cvwinx)
  }
  dim=temp_dim2[which.min(cvdim)]
  return(list(const=cvlambda[which(temp_dim==dim)[1]],dim=dim))
}


