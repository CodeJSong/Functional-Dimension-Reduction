require(RcppArmadillo)
Rcpp::sourceCpp("sub_codes.cpp")
Rcpp::sourceCpp('ftn_approx.cpp')
source('sub_nafpca.R')
#####################################################
# Functional Additive PCA
#####################################################
# p = 1
# x : n X (p*nt) matrix   vectorized version
# t : time point vector
# p : dimension of functions
#####################################################
# p = 2
# x : a list
#    x$x1, x$x2, ..., x$xp
#    tt$x1, ..., tt$xp
#####################################################

nafpca = function(x=NULL, tt=NULL, p=1, type="ftn", unbalanced=FALSE,
                 basisname="bspline", nbasis=11,
                  cvx=FALSE, randomsearch=TRUE, ftnstruct="fda",
                 m1=10,m2=10,
                  shx=.5,ex=0){
  require(fda)
  start.time = Sys.time()
  if(type=="multivariate"){  # x : n x p matrix
    n = nrow(x)
    normX = normmat(x)
    gx = gam(normX)/shx
    KX = exp(-gx*normX)
    one = matrix(1, n, 1)
    Q = diag(n) - one %*% t(one) / n
    lambdaout=NULL
    basis=NULL;xcoef=NULL  ## those are not needed in vectors
  }else if(type=="ftn"){
    if(type!="2dsup")
      temp = get.ftn(x, tt=tt, type=type, p=p, unbalanced=unbalanced,
                       basis = basisname, nbasis=nbasis, lambdascale=2)
    n = temp$n
    xcoef = temp$coef
    xip = temp$xip;
    normX = temp$normx;
    GB = temp$binnerprod;
    basis = temp$basis
    lambdaout = temp$lambdaout


    if(p==1){

      one = matrix(1, n, 1)
      Q = diag(n) - one %*% t(one) / n
      gx = gam(normX)/shx
      KX = exp(-gx*normX)
      GX = Q %*% KX %*% Q
      GX.inv.half = matpower2(GX,-0.5, 0)

      eigenM = eigen(GX, sym=TRUE)
      Mat = eigenM$vec
      pred = t(GX)%*% GX.inv.half %*% Mat
    }else if(p>1){

      one = matrix(1, n, 1)
      Q = diag(n) - one %*% t(one) / n

      normxs = temp$normxs
      KXs = matrix(0, n*p, n*p)
      gxs = rep(0,p)
      j = 1
      normX = normxs[,,j]
      gxs[j] = gam(normX)/shx
      KXj = exp(-gxs[j]*normX); KXs[((j-1)*n+1):(j*n), ((j-1)*n+1):(j*n)] = KXj
      GXj = Q %*% KXj %*% Q
      D.inv.half = matpower2(GXj,-0.5,0); GXs= GXj
      M.half = matpower2(GXj, 0.5,0)
      for(j in 2:p){
        normX = normxs[,,j]
        gxs[j] = gam(normX)/shx
        KXj = exp(-gxs[j]*normX); KXs[((j-1)*n+1):(j*n), ((j-1)*n+1):(j*n)] = KXj
        GXj = Q %*% KXj %*% Q
        D.inv.half = as.matrix(bdiag(D.inv.half, matpower2(GXj,-0.5,0)))
        GXs = rbind(GXs, GXj)
        M.half = rbind(M.half, matpower2(GXj, 0.5,0))
      }

      M.half = M.half %*% t(M.half)
      eigenM = eigen(M.half, sym=TRUE)
      Mat = eigenM$vec
      pred = t(GXs) %*% D.inv.half %*% Mat
      gx = gxs
      GX.inv.half = D.inv.half
    }
  }else if(type=="2dsup"){
    temp=NULL
    # input : vectorized   x = n * nt^2
    n = dim(x)[1]
    one = matrix(1, n, 1)
    Q = diag(n) - one %*% t(one) / n
    tmp=dim(x)[2]; n1 = n2 = sqrt(tmp[1])

    U = fourier.basis(1:n1/n1, m1)
    V = fourier.basis(1:n2/n2, m2)
    m1 = dim(U)[1]; m2= dim(V)[1]   # redefine m1, m2 just in case pre-defiend m1 and m2 are even numbers.

    UU.inv = matpower(U %*% t(U), -1)
    VV.inv = matpower(V %*% t(V), -1)

    coefmat =  kronecker(VV.inv %*% V, UU.inv %*% U)
    normX = normmat.tensor(x, U,V, UU.inv, VV.inv)

    gx = gam(normX)/shx
    KX = exp(-gx*normX)
    GX = Q %*% KX %*% Q
    GX.inv.half = matpower2(GX,-0.5, 0)

    eigenM = eigen(GX, sym=TRUE)
    Mat = eigenM$vec
    pred = t(GX)%*% GX.inv.half %*% Mat
  }
  comp.time = Sys.time()-start.time
  out = list(pred=pred, comp.time=comp.time, eval = eigenM$val, ftn=temp, gx=gx, GX.inv.half = GX.inv.half, Mat=Mat)
  return(out)
}


fpca = function(ftn){
  temp = ftn
  p = temp$p
  n = temp$n
  xcoef = temp$coef
  xip = temp$xip;
  normX = temp$normx;
  GB = temp$binnerprod;
  basis = temp$basis
  lambdaout = temp$lambdaout
  one = matrix(1, n, 1)
  Q = diag(n) - one %*% t(one) / n
  B.half = matpower2(GB,0.5,0)
  if(p==1){
    B.half = matpower2(GB,0.5,0)
    Sigma = B.half %*% xcoef %*% Q %*% t(xcoef) %*% B.half / n
    egn = eigen(Sigma, sym=TRUE)
    B.inv.half = matpower2(GB,-0.5,0)
    pred = Q %*% t(xcoef) %*% GB %*% B.inv.half %*% egn$vec
    out = list(pred=pred, eval=egn$val, mat = GB %*% B.inv.half %*% egn$vec)
  }else if(p>1){ # BX is the same for now...
    M.half = B.half %*% xcoef[,,1] %*% Q
    B.inv.half = matpower2(GB,-0.5,0)
    D.inv.half = B.inv.half
    BB = GB %*% xcoef[,,1] %*% Q
    for(j in 2:p){
      D.inv.half = as.matrix(bdiag(D.inv.half, B.inv.half))
      M.half = rbind(M.half, B.half %*% xcoef[,,j] %*% Q)
      BB = rbind(BB, GB %*% xcoef[,,j] %*% Q)
    }
    DD = M.half
    M.half = M.half %*% t(M.half)/n
    egn = eigen(M.half, sym=TRUE)
    pred = t(BB) %*% D.inv.half %*% egn$vec
    out = list(pred=pred, eval=egn$val)
  }
  return(out)
}

fpca.predict = function(fpcobj, ftn, x, tt){
  # p==1
  temp = ftn
  p = temp$p
  n = temp$n
  xcoef = temp$coef
  xip = temp$xip;
  normX = temp$normx;
  GB = temp$binnerprod;
  basis = temp$basis
  lambdaout = temp$lambdaout

  bname = basis$type
  nbasis = basis$nbasis

  mat = fpcobj$mat

  temp = get.ftn(x, tt, p=p, basisname=bname, nbasis=nbasis)
  test.coef = temp$coef
  n = ncol(test.coef)
  one = matrix(1, n, 1)
  Q = diag(n) - one %*% t(one) / n
  pred = Q %*% t(test.coef) %*% mat
  return(pred)
}

