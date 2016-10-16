require(RcppArmadillo)
require(Rcpp)
library(inline)
# Load Subcodes
##############################################
# 1.  some C codes those have for-loop
# 2. main function (fdr)
# 3. codes related to fdr like dimensino selection, prediction of test set, FSIR..
##############################################

# some c codes
#########################################
# favelopp : needed in f-gsave iteration
#########################################
faveloopsrc='
Function matpower2("matpower2");
Function matpower("matpower");
double ex=as<double>(epx), ey=as<double>(epy), yclas=as<double>(yclass);
arma::mat GX = as<arma::mat>(gx), GY = as<arma::mat>(gy), KY = as<arma::mat>(ky), KX = as<arma::mat>(kx);
arma::mat Q=as<arma::mat>(q), Ay(GX.n_rows,GX.n_cols), Vfave(GX.n_rows,GX.n_cols); Vfave.zeros();
arma::mat GXinv=as<arma::mat>(matpower2(GX, -1, ex));
arma::mat GXinvhalf=as<arma::mat>(matpower2(GX, -0.5, ex));
arma::mat KYinv;
if(yclas==0){
KYinv=as<arma::mat>(matpower2(KY, -1, ey));
for(int i=0;i<GX.n_rows;i++){
Ay = KX * (diagmat(KYinv * KY.col(i)) - KYinv * KY.col(i) * KY.col(i).t() * KYinv ) * KX;
Vfave = Vfave + (Q - GXinv * Ay * GXinv) * GX * (Q - GXinv * Ay * GXinv);
}
}
if(yclas==1){
KYinv=as<arma::mat>(matpower(KY, -1));
for(int i=0;i<GX.n_rows;i++){
Ay = KX * (diagmat(KYinv * KY.col(i)) - KYinv * KY.col(i) * KY.col(i).t() * KYinv ) * KX;
Vfave = Vfave + (Q - GXinv * Ay * GXinv) * GX * (Q - GXinv * Ay * GXinv);
}
}
Vfave = Q * GXinvhalf * Vfave * GXinvhalf * Q;
return(wrap(Vfave));'

faveloop <- cxxfunction(signature(gx="numeric",gy="numeric", ky="numeric", kx="numeric",q="numeric",epx="numeric", epy="numeric", yclass="numeric"),faveloopsrc, plugin="RcppArmadillo")


#########################################
# cvx : needed in finding ex gammax, ey, gammay
#########################################
cvxsrc <- '
Function gam("gam");
Function matpower2("matpower2");
Function whichmin("which.min");
int indexshx, indexshy, indexex, indexey;
arma::vec shxgrid = as<arma::vec>(cvwinx);
arma::vec shygrid = as<arma::vec>(cvwiny);
arma::vec exgrid = as<arma::vec>(cvex),  eygrid=exgrid;
arma::vec cv_shx(shxgrid.n_rows), cv_shy(shygrid.n_rows), cv_ex(exgrid.n_rows), cv_ey(eygrid.n_rows);
cv_shx.zeros();cv_shy.zeros();cv_ex.zeros();cv_ey.zeros();
arma::mat normx = as<arma::mat>(normX), normy = as<arma::mat>(normY);
arma::mat KX, KY, KXi, KYi, KXii, KYii, KXiinv;
double gx, gy, shxout, shyout, exout, eyout;


for(int o=0; o < shxgrid.n_rows; o++){
  gx = as<double>(gam(normx))/shxgrid(o); 
  gy = as<double>(gam(normy))/shygrid(0);
  KX = exp(-gx*normx); KY =exp(-gy*normy);
  for(int i=0;i<KX.n_rows;i++){
    KXi=KX; KYi=KY; KXii=KX; KYii=KY;
    KXii.shed_row(i); KXii.shed_col(i);
    KYii.shed_row(i); KYii.shed_col(i);
    KXi.shed_col(i); KYi.shed_col(i);
    KXiinv = as<arma::mat>(matpower2(KXii, -1, exgrid(0)));
    cv_shx(o)=cv_shx(o) +  as_scalar((KYi.row(i) - KXi.row(i) * KXiinv * KYii) * ( KYi.row(i) - KXi.row(i) * KXiinv * KYii    ).t());
  }
}
indexshx = as<int>(whichmin(cv_shx)) - 1;
shxout = shxgrid(indexshx);


gx = as<double>(gam(normx))/shxout; 
gy = as<double>(gam(normy))/shygrid(0);
KX = exp(-gx*normx); KY =exp(-gy*normy);
for(int o=0; o < exgrid.n_rows; o++){
  for(int i=0;i<KX.n_rows;i++){
    KXi=KX; KYi=KY; KXii=KX; KYii=KY;
    KXii.shed_row(i); KXii.shed_col(i);
    KYii.shed_row(i); KYii.shed_col(i);
    KXi.shed_col(i); KYi.shed_col(i);
    KXiinv = as<arma::mat>(matpower2(KXii, -1, exgrid(o)));
    cv_ex(o)=cv_ex(o) +  as_scalar((KYi.row(i) - KXi.row(i) * KXiinv * KYii) * ( KYi.row(i) - KXi.row(i) * KXiinv * KYii    ).t());
  }
}
indexex = as<int>(whichmin(cv_ex)) - 1;
exout = exgrid(indexex);

double yclas = as<double>(yclass);
// if y is not categorical
if(yclas==0){
  for(int o=0; o < shygrid.n_rows; o++){  
    gy = as<double>(gam(normy))/shygrid(o);
    KY =exp(-gy*normy);
    for(int i=0;i<KX.n_rows;i++){     // Exchange KX and KY
                                      KXi=KY; KYi=KX; KXii=KY; KYii=KX;
                                      KXii.shed_row(i); KXii.shed_col(i);
                                      KYii.shed_row(i); KYii.shed_col(i);
                                      KXi.shed_col(i); KYi.shed_col(i);
                                      KXiinv = as<arma::mat>(matpower2(KXii, -1, exout));
                                      cv_shy(o)=cv_shy(o) +  as_scalar((KYi.row(i) - KXi.row(i) * KXiinv * KYii) * ( KYi.row(i) - KXi.row(i) * KXiinv * KYii    ).t());
    }
  }
  indexshy = as<int>(whichmin(cv_shy)) - 1;
  shyout = shygrid(indexshy);
  gy = as<double>(gam(normy))/shyout; KY =exp(-gy*normy);
for(int o=0; o < eygrid.n_rows; o++){  
  for(int i=0;i<KX.n_rows;i++){     // Exchange KX and KY
    KXi=KY; KYi=KX; KXii=KY; KYii=KX;
    KXii.shed_row(i); KXii.shed_col(i);
    KYii.shed_row(i); KYii.shed_col(i);
    KXi.shed_col(i); KYi.shed_col(i);
    KXiinv = as<arma::mat>(matpower2(KXii, -1, eygrid(o)));
    cv_ey(o)=cv_ey(o) +  as_scalar((KYi.row(i) - KXi.row(i) * KXiinv * KYii) * ( KYi.row(i) - KXi.row(i) * KXiinv * KYii    ).t());
  }
}
indexey = as<int>(whichmin(cv_ey)) - 1;
eyout = eygrid(indexey);

}
if(yclas==1){  // if y is categorical
  KY = as<arma::mat>(ky);
  shyout = shygrid(0);
}



return(wrap(List::create(Named("shx",shxout),Named("shy",shyout),Named("ex",exout), Named("ey", eyout))));
'
cvxsearch <- cxxfunction(signature(cvwinx="numeric", cvwiny="numeric", cvex="numeric",
                                    normX="numeric", normY="numeric", yclass="numeric", ky="numeric"),
                         cvxsrc, plugin="RcppArmadillo")


#########################################
# cvt : needed in finding et gammat
#########################################
cvtsrc <- '
Function gam("gam");
Function matpower2("matpower2");
Function whichmin("which.min");
List x(xx), tt(ttt);
 // kernel = 1  : gaussian, kernel = brownian
int teq=as<int>(tequal), dim=as<int>(dimension), n=as<int>(num), indexsht, indexet, kern=as<int>(kernel);
arma::vec shtgrid = as<arma::vec>(cvwint);  // if kernel==brownian, skip it
arma::vec etgrid = as<arma::vec>(cvet);
arma::vec cv_sht(shtgrid.n_rows), cv_et(etgrid.n_rows); cv_sht.zeros(); cv_et.zeros();
arma::mat normt = as<arma::mat>(normmatt), KT, KTinv, H;   // if kernel==brownian, normt=KT
double gt, shtout, etout;
arma::vec tuni = as<arma::vec>(tuniq);
if(dim==1){
  arma::vec xi;
  double nt = normt.n_rows;
  if(kern==1){   // do it if kernel==gaussian
    for(int o=0;o<shtgrid.n_rows;o++){
      gt = as<double>(gam(normt))/shtgrid(o);
      KT = exp(-gt*normt);
      KTinv = as<arma::mat>(matpower2(KT, -1, etgrid(0)));
      H = KTinv * KT;
      
      for(int i=0;i<n;i++){
        if(teq==0){
          arma::vec ti = as<arma::vec>(tt[i]);
          int nti=ti.n_rows;
          arma::vec index(nti);
        
          for(int k=0;k<nti;k++){
            index(k) = as<int>(wrap(find(tuni==ti(k))));
          }
          arma::uvec indexi = arma::conv_to<arma::uvec>::from(index);
          arma::mat Vii = KT.submat((indexi),(indexi));
          arma::mat Viiinv = as<arma::mat>(matpower2(Vii, -1, etgrid(0)));
          xi=as<arma::vec>(x[i]);
          H = Viiinv * Vii;
          nt = nti;
        }
        
        if(teq==1){
          xi=as<arma::vec>(x[i]);
        }
        cv_sht(o) = cv_sht(o) + pow(norm(H * xi - xi, 2), 2) / nt / pow(1- trace(H)/nt  , 2 );
      }
    }
  }
  
  indexsht = as<int>(whichmin(cv_sht)) - 1;
  shtout = shtgrid(indexsht);
  for(int o=0;o<shtgrid.n_rows;o++){
    if(kern==1){  // if kernel==gaussian
      gt = as<double>(gam(normt))/shtout;
      KT = exp(-gt*normt);
    }
    if(kern==0){   // if kernel=brownian
      KT = normt;
    }
    KTinv = as<arma::mat>(matpower2(KT, -1, etgrid(o)));
    H = KTinv * KT;
    for(int i=0;i<n;i++){
      if(teq==0){
        arma::vec ti = as<arma::vec>(tt[i]);
        int nti=ti.n_rows;
        arma::vec index(nti);
      
        for(int k=0;k<nti;k++){
          index(k) = as<int>(wrap(find(tuni==ti(k))));
        }
        arma::uvec indexi = arma::conv_to<arma::uvec>::from(index);
        arma::mat Vii = KT.submat((indexi),(indexi));
        arma::mat Viiinv = as<arma::mat>(matpower2(Vii, -1, etgrid(o)));
        xi=as<arma::vec>(x[i]);
        H = Viiinv * Vii;
        nt = nti;
      }
      if(teq==1){
        xi=as<arma::vec>(x[i]);
      }

      cv_et(o) = cv_et(o) + pow(norm(H * xi - xi, 2), 2) / nt / pow(1- trace(H)/nt  , 2 );
    }
  }
  indexet = as<int>(whichmin(cv_sht)) - 1;
  etout = etgrid(indexet);
}
if(dim>1){
  arma::mat xi;
  double nt = normt.n_rows;
  if(kern==1){   // do it if kernel==gaussian
     for(int o=0;o<shtgrid.n_rows;o++){
       gt = as<double>(gam(normt))/shtgrid(o);
       KT = exp(-gt*normt);
       KTinv = as<arma::mat>(matpower2(KT, -1, etgrid(0)));
       H = KTinv * KT;
       
       for(int i=0;i<n;i++){
         if(teq==0){
           arma::vec ti = as<arma::vec>(tt[i]);
           int nti=ti.n_rows;
           arma::vec index(nti);
           
           for(int k=0;k<nti;k++){
             index(k) = as<int>(wrap(find(tuni==ti(k))));
           }
           arma::uvec indexi = arma::conv_to<arma::uvec>::from(index);
           arma::mat Vii = KT.submat((indexi),(indexi));
           arma::mat Viiinv = as<arma::mat>(matpower2(Vii, -1, etgrid(0)));
           xi=as<arma::mat>(x[i]);
           H = Viiinv * Vii;
           nt = nti;
         }
         
         if(teq==1){
           xi=as<arma::vec>(x[i]);
         }
         cv_sht(o) = cv_sht(o) + pow(norm(H * xi - xi, 2), 2) / nt / pow(1- trace(H)/nt  , 2 );
       }
     }
  }
  
  indexsht = as<int>(whichmin(cv_sht)) - 1;
  shtout = shtgrid(indexsht);
  for(int o=0;o<shtgrid.n_rows;o++){
    if(kern==1){  // if kernel==gaussian
                  gt = as<double>(gam(normt))/shtout;
                  KT = exp(-gt*normt);
    }
    if(kern==0){   // if kernel=brownian
                   KT = normt;
    }
    KTinv = as<arma::mat>(matpower2(KT, -1, etgrid(o)));
    H = KTinv * KT;
    for(int i=0;i<n;i++){
      if(teq==0){
        arma::vec ti = as<arma::vec>(tt[i]);
        int nti=ti.n_rows;
        arma::vec index(nti);
        
        for(int k=0;k<nti;k++){
          index(k) = as<int>(wrap(find(tuni==ti(k))));
        }
        arma::uvec indexi = arma::conv_to<arma::uvec>::from(index);
        arma::mat Vii = KT.submat((indexi),(indexi));
        arma::mat Viiinv = as<arma::mat>(matpower2(Vii, -1, etgrid(o)));
        xi=as<arma::mat>(x[i]);
        H = Viiinv * Vii;
        nt = nti;
      }
      if(teq==1){
        xi=as<arma::vec>(x[i]);
      }
      
      cv_et(o) = cv_et(o) + pow(norm(H * xi - xi, 2), 2) / nt / pow(1- trace(H)/nt  , 2 );
    }
  }
  indexet = as<int>(whichmin(cv_sht)) - 1;
  etout = etgrid(indexet);
}
return(wrap(List::create(Named("sht",shtout),Named("et",etout))));
'
cvtsearch <- cxxfunction(signature(num="numeric", kernel="numeric", cvwint="numeric",cvet="numeric", normmatt="numeric", xx="list", ttt="list", tuniq="numeric",
                                   tequal="numeric", dimension="numeric"),cvtsrc,plugin="RcppArmadillo")


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

#########################################
# cvdsearch : needed in finding the dimension (old version) . the current version is cvd2search
#########################################
cvdsrc <- '
Function gam("gam");
Function matpower2("matpower2");
Function normmat("normmat");
Function whichmin("which.min");
arma::mat xmat = as<arma::mat>(xxmat), KY = as<arma::mat>(ky);
arma::mat KX, KXi, KYi, KXii, KYii, KXiinv;
arma::mat normx=as<arma::mat>(normmat(xmat));
double cvd=0,gx, ex, shxout, exout;
arma::vec shxgrid = as<arma::vec>(cvwinx);
arma::vec exgrid = as<arma::vec>(cvex);
arma::vec cv_shx(shxgrid.n_rows),  cv_ex(exgrid.n_rows);
cv_shx.zeros();cv_ex.zeros();

int indexshx, indexex;

for(int o=0; o < shxgrid.n_rows; o++){
  gx = as<double>(gam(normx))/shxgrid(o);
KX = exp(-gx*normx); 
for(int i=0;i<KX.n_rows;i++){
KXi=KX; KYi=KY; KXii=KX; KYii=KY;
KXii.shed_row(i); KXii.shed_col(i);
KYii.shed_row(i); KYii.shed_col(i);
KXi.shed_col(i); KYi.shed_col(i);
KXiinv = as<arma::mat>(matpower2(KXii, -1, exgrid(0)));
cv_shx(o)=cv_shx(o) +  as_scalar((KYi.row(i) - KXi.row(i) * KXiinv * KYii) * ( KYi.row(i) - KXi.row(i) * KXiinv * KYii    ).t());
}
}
indexshx = as<int>(whichmin(cv_shx)) - 1;
shxout = shxgrid(indexshx);


gx = as<double>(gam(normx))/shxout; 
KX = exp(-gx*normx);
for(int o=0; o < exgrid.n_rows; o++){
for(int i=0;i<KX.n_rows;i++){
KXi=KX; KYi=KY; KXii=KX; KYii=KY;
KXii.shed_row(i); KXii.shed_col(i);
KYii.shed_row(i); KYii.shed_col(i);
KXi.shed_col(i); KYi.shed_col(i);
KXiinv = as<arma::mat>(matpower2(KXii, -1, exgrid(o)));
cv_ex(o)=cv_ex(o) +  as_scalar((KYi.row(i) - KXi.row(i) * KXiinv * KYii) * ( KYi.row(i) - KXi.row(i) * KXiinv * KYii    ).t());
}
}
indexex = as<int>(whichmin(cv_ex)) - 1;
exout = exgrid(indexex);


ex= exout;

KX = exp(-gx*normx); 
for(int i=0;i<KX.n_rows;i++){
KXi=KX; KYi=KY; KXii=KX; KYii=KY;
KXii.shed_row(i); KXii.shed_col(i);
KYii.shed_row(i); KYii.shed_col(i);
KXi.shed_col(i); KYi.shed_col(i);
KXiinv = as<arma::mat>(matpower2(KXii, -1, ex));
cvd=cvd +  as_scalar((KYi.row(i) - KXi.row(i) * KXiinv * KYii) * ( KYi.row(i) - KXi.row(i) * KXiinv * KYii    ).t());
}

return(wrap(cvd));
'
cvdsearch <- cxxfunction(signature(xxmat="numeric",  ky="numeric", cvex="numeric", cvwinx="numeric"),
                         cvdsrc, plugin="RcppArmadillo")


#########################################
# cvd2search : needed in finding the dimension (new vesion). old is not used 
#########################################
cvd2src <- '
Function gam("gam");
Function matpower("matpower");
Function normmat("normmat");
Function whichmin("which.min");
arma::mat xmat = as<arma::mat>(xxmat), KY = as<arma::mat>(ky);
arma::mat KX, KXi, KYi, KXii, KYii, KXiinv;
arma::mat normx=as<arma::mat>(normmat(xmat));
double cvd=0,gx, shxout;
arma::vec shxgrid = as<arma::vec>(cvwinx);
// arma::vec cv_shx(shxgrid.n_rows);
// cv_shx.zeros();
/*
int indexshx;

for(int o=0; o < shxgrid.n_rows; o++){
gx = as<double>(gam(normx))/shxgrid(o);
KX = exp(-gx*normx); 
for(int i=0;i<KX.n_rows;i++){
KXi=KX; KYi=KY; KXii=KX; KYii=KY;
KXii.shed_row(i); KXii.shed_col(i);
KYii.shed_row(i); KYii.shed_col(i);
KXi.shed_col(i); KYi.shed_col(i);
KXiinv = as<arma::mat>(matpower(KXii, -1));
cv_shx(o)=cv_shx(o) +  as_scalar((KYi.row(i) - KXi.row(i) * KXiinv * KYii) * ( KYi.row(i) - KXi.row(i) * KXiinv * KYii    ).t());
}
}
indexshx = as<int>(whichmin(cv_shx)) - 1;
shxout = shxgrid(indexshx);
*/

shxout=2;

gx = as<double>(gam(normx))/shxout; 
KX = exp(-gx*normx);
for(int i=0;i<KX.n_rows;i++){
KXi=KX; KYi=KY; KXii=KX; KYii=KY;
KXii.shed_row(i); KXii.shed_col(i);
KYii.shed_row(i); KYii.shed_col(i);
KXi.shed_col(i); KYi.shed_col(i);
KXiinv = as<arma::mat>(matpower(KXii, -1));
cvd=cvd +  as_scalar((KYi.row(i) - KXi.row(i) * KXiinv * KYii) * ( KYi.row(i) - KXi.row(i) * KXiinv * KYii    ).t());
}

return(wrap(cvd));
'
cvd2search <- cxxfunction(signature(xxmat="numeric",  ky="numeric", cvwinx="numeric"),
                         cvd2src, plugin="RcppArmadillo")
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
        cvtout <- cvtsearch(num=length(x), kernel=as.numeric(kernel=="gaussian"),
                            cvwint,cvet, normmatt=KTX, xx=x, ttt=tx, tuniq=xt.unique,
                            tequal=as.numeric(xtequal), dimension=dimx)
        et = cvtout$et
      }
      
      if(kernel=="gaussian"){
        
        cvtout <- cvtsearch(num=length(x), kernel=as.numeric(kernel=="gaussian"),
                            cvwint,cvet, normmatt=normmat.t, xx=x, ttt=tx, tuniq=xt.unique,
                            tequal=as.numeric(xtequal), dimension=dimx)
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
        cvtout <- cvtsearch(num=length(y), kernel=as.numeric(kernel=="gaussian"),
                            cvwint,cvet, normmatt=KTY, xx=y, ttt=ty, tuniq=yt.unique,
                            tequal=as.numeric(ytequal), dimension=dimy)
        ety = cvtout$et
      }
      if(kernel=="gaussian"){
        cvtout <- cvtsearch(num=length(y), kernel=as.numeric(kernel=="gaussian"),
                            cvwint,cvet, normmatt=normmat.ty, xx=y, ttt=ty, tuniq=yt.unique,
                            tequal=as.numeric(ytequal), dimension=dimy)
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
    #     Vfave <- matrix(0, n, n)
    #     GY.inv <- matpower2(GY,-1,ey)
    #     GX.inv <- matpower2(GX,-1,ex)
    #     
    #     KY.inv <- matpower2(KY,-1,ey)
    #     for(i in 1:n){ #Y_i
    #       by <- KY[,i]
    #       Ay <- KX %*% diag(as.vector(KY.inv %*% by)) %*% KX - KX %*% KY.inv %*% by %*% t(by) %*% KY.inv %*% KX
    #       Vfave <- Vfave + (Q- GX.inv %*% Ay %*% GX.inv) %*% GX %*%  (Q- GX.inv %*% Ay %*% GX.inv)
    #     }
    #   
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


# fkpca ------

fkpca <- function(xdata, xftn=TRUE,  xtequal=FALSE, dim=2,
                  option="fir", CV=FALSE, kernel="brownian", c_brown=1,
                  d=1, sht=2, shx=2, shy=2,ex=0.05, et=0.05, ey=0.05){
  id <- unique(xdata$id)
  xt.unique <- sort(unique(xdata$tt))
  n <- length(id)
  x <- NULL
  tx <- NULL
  for(i in 1:n){
    temp <- xdata[xdata$id==id[i],]
    tx[[i]] <- temp$tt
    x[[i]] <- temp[,2:3]
  }
  ## Make KX depending on whether X is functional or not
  
  if(xftn){   # x is list of functions
    # input variable is little different.
    # x : list of xi(tij)
    # y : n-dimension
    # tt : list of tij
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
    xip <- x.ip2(x=x, tt=tx, t.unique=xt.unique, KT=KTX, et=et, dim=dim)  ## XX_ij = <Xi, Xj>_{HT}
    normX = normmat(xip=xip, ftn=TRUE)
    gx <- gam(normX)/shx
    KX <- exp(-gx*normX)    
    #ex <- ex * n / 2
  }
  
  ## Estimate central class 
  M <- matrix(0, n, n)
  one <- matrix(1, n, 1)
  Q <- diag(n) - one %*% t(one) / n    
  GX <- Q %*% KX %*% Q
  GX.half.inv <- matpower2(GX, -1, ex)
  M <- GX.half.inv %*% GX %*% GX %*% GX.half.inv
  Mat <- eigen(M)$vec
  pred <- GX.half.inv %*% Mat
  return(pred)
}

################################
# GCV for bspline function estimation  (not used in our paper)
################################
gcvsrc <- 
  '
arma::mat X = as<arma::mat>(xmat), R = as<arma::mat>(rr);
arma::mat B(X.n_cols, X.n_cols), Binv(X.n_cols, X.n_cols), Iden_B = arma::eye<arma::mat>(X.n_cols,X.n_cols);
arma::mat H(X.n_rows, X.n_rows), Iden_H = arma::eye<arma::mat>(X.n_rows, X.n_rows);
arma::vec y = as<arma::vec>(ystar);
NumericVector lambdagrd = as<NumericVector>(lambdagrid);
arma::vec gcv(lambdagrd.size()), mucoef(X.n_cols);
int mincv;

Function whichmin("which.min");
for(int i=0; i <lambdagrd.size(); i++){
B = X.t() * X + lambdagrd(i) * R;
Binv = inv(B);
H = X * Binv * X.t();
gcv[i] = arma::as_scalar(y.t() * (Iden_H - H) *(Iden_H - H) * y )/ pow((1-trace(H)/X.n_rows),2.0);
}
mincv = as<int>(whichmin(gcv));
B = X.t() * X + lambdagrd(mincv-1) * R;
Binv = inv(B);
mucoef = Binv * X.t() * y;
return(List::create(Named("gcv",gcv), Named("mucoef",mucoef),Named("lambdas",lambdagrd),
Named("lambda",lambdagrd(mincv-1))));
'
bsplinegcv <- cxxfunction(signature(xmat="numeric",rr="numeric",ystar="numeric",lambdagrid="numeric"),gcvsrc, plugin="RcppArmadillo")
################################
# eval.basis   : c code for eval.basis function so that I can manipulate. (not used in our paper)
################################
evalbasis<-'
arma::vec tt = as<arma::vec>(t), knot=as<arma::vec>(knots);
arma::vec dervis = as<arma::vec>(derivs);
Function splinedesign("splineDesign");
return(splinedesign(knot,tt,4, derivs));
'
evalbasis<-cxxfunction(signature(t="numeric", knots="numeric", derivs="numeric"), evalbasis, plugin="RcppArmadillo")



##################################################################
# fsir : FSIR  - classic. not regularized. 
##################################################################

fsir <- function(x, y, H, k){
  n <- length(y)
  J <- ncol(x)
  z <- center(x)
  ### Make a slice  # of slice is H
  d <- (max(y)-min(y))/H   ## length of slice
  m <- min(y)
  yindex <- rep(0,n)         # yy[i]=h implies i'th observation is in slice h out of H
  yindex <- (y<=m+d)*1
  yindex <- yindex+(y>m+(H-1)*d)*H
  for ( i in 2:(H-1))
  {
    yindex <- yindex+(m+(i-1)*d<y)*(y<=m+i*d)*i
  }
  ####  Calculate ph, mh, V , mh[h,] is a vector
  ph <- rep(0, H) #ps
  hh <- rep(0, J) # hs
  R <- t(z)%*%z/n # RnJn
  AA<-matrix(0,J,J)
  V <- matrix(0, J, J)
  for ( i in 1:H)
  {
    ph[i] <- sum(yindex==i)/n
    for (j in 1:J)
      if(ph[i]!=0) hh[j] <- mean(z[yindex==i,j])
      V <- V + ph*hh%*%t(hh)
  }
  # k : reduced dimension
  
  p <- eigen(R)$vec
  Pk <- p[,1:k]%*%((solve(p%*%t(p))%*%t(p))[1:k,])
  Rk <- Pk %*% R %*% Pk
  Rk.inv.half <- matpower(Rk, -1/2)
  M <- Rk.inv.half %*% V %*% Rk.inv.half
  beta <- eigen(M,sym=TRUE)$vec[,1:k]
  xi <- x%*%beta
  return(list(beta=beta,pred=xi))
}


##################################################################
# fsir.reg : regularized fsir. there is a regularized parameter when we get function estimation
##################################################################

fsir.reg <- function(x, tt, y, H=4, alpha=10){
  require(fda)
  # x : matrix n x nt
  # tt : vector of length nt
  # y : scalar  // category
  # H : # of category = length(y) - 1
  n <- length(y)
  nbasis <- length(tt);nt <- length(tt)
  tmin <- min(tt); tmax <- max(tt)
  basis<-create.bspline.basis(c(tmin,tmax),nbasis,4)
  B <- t(eval.basis(tt,basis))
  R <- eval.basis(tt, basis,2)
  RR <- crossprod(R) / nt
  Me <- matrix(0, nbasis, nbasis)
  for(h in 1:H){
    indexy <- as.numeric(y==h)
    nh <- sum(indexy)
    Me = Me + nh/n * B %*% t(x) %*% indexy %*% t(indexy) %*% x %*% t(B)
  }
  Mxa <- B %*% t(x) %*% x %*% t(B) + alpha * RR
  Mxa.inv.half <- matpower(Mxa, -1/2)
  edcomp <- eigen(Mxa.inv.half , sym=TRUE)
  A <- Mxa.inv.half %*% edcomp$vec 
  pred <- x %*% t(B) %*% A
  return(list(pred=pred, Mat=t(B) %*% A))
}

##################################################################
# dim.search : old veresion dimension selection
##################################################################


dim.search <- function(nfdr_pred, ncvl=10){
  cvlambda <- exp(seq(log(0.01),log(1),len=ncvl))
  cvwinx <- c(2,exp(seq(log(0.05),log(20),len=ncvl-1)))
  n = nrow(nfdr_pred$pred)
  cvex <-  c(n^(-1/4),exp(seq(log(n^(-1/4)/50),log(n^(-1/4)*50),len=ncvl-1)))
  temp_dim <- rep(0, ncvl);cvdim<- rep(0, ncvl);
  for(k in 1:ncvl){
    temp_dim[k] <- which.max(fdr.bic(nfdr_pred$eval, fun1=fun1,fun2=fun2, c=cvlambda[k]))
  }
  temp_dim <- sort(unique(temp_dim))
  for(k in 1:length(temp_dim)){
    cvdim[k] <- cvdsearch(xxmat=as.matrix(nfdr_pred$pred[,1:temp_dim[k]]),ky=nfdr_pred$KY, cvex=cvex, cvwinx=cvwinx)
  }
  return(temp_dim[which.min(cvdim)])
  #return(list(temp_dim[which.min(cvdim)], temp_dim, cvdim))
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
    cvdim[k] <- cvd2search(xxmat=as.matrix(nfdr_pred$pred[,1:temp_dim2[k]]),ky=nfdr_pred$KY, cvwinx=cvwinx)
  }
  dim=temp_dim2[which.min(cvdim)]
  return(list(const=cvlambda[which(temp_dim==dim)[1]],dim=dim))
}

