// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// faveloop : implement for-loop with compiled C++
// [[Rcpp::export]]
mat faveloop(mat GX, mat GY, mat KY, mat KX, mat Q, double ex, double ey, double yclas){
  Function matpower2("matpower2");
  Function matpower("matpower");
  mat Ay(GX.n_rows,GX.n_cols), Vfave(GX.n_rows,GX.n_cols); Vfave.zeros();
  mat GXinv=as<mat>(matpower2(GX, -1, ex)), GXinvhalf=as<mat>(matpower2(GX, -0.5, ex));;
  mat KYinv;
  if(yclas==0){
    KYinv=as<mat>(matpower2(KY, -1, ey));
    for(int i=0;i<GX.n_rows;i++){
      Ay = KX * (diagmat(KYinv * KY.col(i)) - KYinv * KY.col(i) * KY.col(i).t() * KYinv ) * KX;
      Vfave = Vfave + (Q - GXinv * Ay * GXinv) * GX * (Q - GXinv * Ay * GXinv);
    }
  }
  if(yclas==1){   // if y is categorical
    KYinv=as<mat>(matpower(KY, -1));
    for(int i=0;i<GX.n_rows;i++){
      Ay = KX * (diagmat(KYinv * KY.col(i)) - KYinv * KY.col(i) * KY.col(i).t() * KYinv ) * KX;
      Vfave = Vfave + (Q - GXinv * Ay * GXinv) * GX * (Q - GXinv * Ay * GXinv);
    }
  }
  Vfave = Q * GXinvhalf * Vfave * GXinvhalf * Q;
  
  return(Vfave);
}


// tuning parameter selection about X
// [[Rcpp::export]]
List cvxsearch(vec shxgrid, vec shygrid, vec exgrid, mat normx, mat normy, double yclas, mat KY){
  Function gam("gam");
  Function matpower2("matpower2");
  Function whichmin("which.min");
  int indexshx, indexshy, indexex, indexey;
  vec eygrid=exgrid;
  vec cv_shx(shxgrid.n_rows), cv_shy(shygrid.n_rows), cv_ex(exgrid.n_rows), cv_ey(eygrid.n_rows);
  cv_shx.zeros();cv_shy.zeros();cv_ex.zeros();cv_ey.zeros();
  
  mat KX, KXi, KYi, KXii, KYii, KXiinv;
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
      KXiinv = as<mat>(matpower2(KXii, -1, exgrid(0)));
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
      KXiinv = as<mat>(matpower2(KXii, -1, exgrid(o)));
      cv_ex(o)=cv_ex(o) +  as_scalar((KYi.row(i) - KXi.row(i) * KXiinv * KYii) * ( KYi.row(i) - KXi.row(i) * KXiinv * KYii    ).t());
    }
  }
  indexex = as<int>(whichmin(cv_ex)) - 1;
  exout = exgrid(indexex);
  
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
                                        KXiinv = as<mat>(matpower2(KXii, -1, exout));
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
      KXiinv = as<mat>(matpower2(KXii, -1, eygrid(o)));
      cv_ey(o)=cv_ey(o) +  as_scalar((KYi.row(i) - KXi.row(i) * KXiinv * KYii) * ( KYi.row(i) - KXi.row(i) * KXiinv * KYii    ).t());
    }
  }
  indexey = as<int>(whichmin(cv_ey)) - 1;
  eyout = eygrid(indexey);
  
  }
  if(yclas==1){  // if y is categorical
    shyout = shygrid(0);
  }
  return(List::create(Named("shx",shxout),Named("shy",shyout),Named("ex",exout), Named("ey", eyout)));
}


// tuning parameter selection about time
// [[Rcpp::export]]
List cvtsearch(int n, int kern, vec shtgrid, vec etgrid, mat normt, List x, List tt, vec tuni, int teq, int dim){
  Function gam("gam");
  Function matpower2("matpower2");
  Function whichmin("which.min");
   // kernel = 1  : gaussian, kernel = brownian
  int indexsht, indexet;
  vec cv_sht(shtgrid.n_rows), cv_et(etgrid.n_rows); cv_sht.zeros(); cv_et.zeros();
  mat KT, KTinv, H;   // if kernel==brownian, normt=KT
  double gt, shtout, etout;
  if(dim==1){
    vec xi;
    double nt = normt.n_rows;
    if(kern==1){   // do it if kernel==gaussian
      for(int o=0;o<shtgrid.n_rows;o++){
        gt = as<double>(gam(normt))/shtgrid(o);
        KT = exp(-gt*normt);
        KTinv = as<mat>(matpower2(KT, -1, etgrid(0)));
        H = KTinv * KT;
        
        for(int i=0;i<n;i++){
          if(teq==0){
            vec ti = as<vec>(tt[i]);
            int nti=ti.n_rows;
            vec index(nti);
          
            for(int k=0;k<nti;k++){
              index(k) = as<int>(wrap(find(tuni==ti(k))));
            }
            uvec indexi = conv_to<uvec>::from(index);
            mat Vii = KT.submat((indexi),(indexi));
            mat Viiinv = as<mat>(matpower2(Vii, -1, etgrid(0)));
            xi=as<vec>(x[i]);
            H = Viiinv * Vii;
            nt = nti;
          }
          
          if(teq==1){
            xi=as<vec>(x[i]);
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
      KTinv = as<mat>(matpower2(KT, -1, etgrid(o)));
      H = KTinv * KT;
      for(int i=0;i<n;i++){
        if(teq==0){
          vec ti = as<vec>(tt[i]);
          int nti=ti.n_rows;
          vec index(nti);
        
          for(int k=0;k<nti;k++){
            index(k) = as<int>(wrap(find(tuni==ti(k))));
          }
          uvec indexi = conv_to<uvec>::from(index);
          mat Vii = KT.submat((indexi),(indexi));
          mat Viiinv = as<mat>(matpower2(Vii, -1, etgrid(o)));
          xi=as<vec>(x[i]);
          H = Viiinv * Vii;
          nt = nti;
        }
        if(teq==1){
          xi=as<vec>(x[i]);
        }
  
        cv_et(o) = cv_et(o) + pow(norm(H * xi - xi, 2), 2) / nt / pow(1- trace(H)/nt  , 2 );
      }
    }
    indexet = as<int>(whichmin(cv_sht)) - 1;
    etout = etgrid(indexet);
  }
  if(dim>1){
    mat xi;
    double nt = normt.n_rows;
    if(kern==1){   // do it if kernel==gaussian
       for(int o=0;o<shtgrid.n_rows;o++){
         gt = as<double>(gam(normt))/shtgrid(o);
         KT = exp(-gt*normt);
         KTinv = as<mat>(matpower2(KT, -1, etgrid(0)));
         H = KTinv * KT;
         
         for(int i=0;i<n;i++){
           if(teq==0){
             vec ti = as<vec>(tt[i]);
             int nti=ti.n_rows;
             vec index(nti);
             
             for(int k=0;k<nti;k++){
               index(k) = as<int>(wrap(find(tuni==ti(k))));
             }
             uvec indexi = conv_to<uvec>::from(index);
             mat Vii = KT.submat((indexi),(indexi));
             mat Viiinv = as<mat>(matpower2(Vii, -1, etgrid(0)));
             xi=as<mat>(x[i]);
             H = Viiinv * Vii;
             nt = nti;
           }
           
           if(teq==1){
             xi=as<vec>(x[i]);
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
      KTinv = as<mat>(matpower2(KT, -1, etgrid(o)));
      H = KTinv * KT;
      for(int i=0;i<n;i++){
        if(teq==0){
          vec ti = as<vec>(tt[i]);
          int nti=ti.n_rows;
          vec index(nti);
          
          for(int k=0;k<nti;k++){
            index(k) = as<int>(wrap(find(tuni==ti(k))));
          }
          uvec indexi = conv_to<uvec>::from(index);
          mat Vii = KT.submat((indexi),(indexi));
          mat Viiinv = as<mat>(matpower2(Vii, -1, etgrid(o)));
          xi=as<mat>(x[i]);
          H = Viiinv * Vii;
          nt = nti;
        }
        if(teq==1){
          xi=as<vec>(x[i]);
        }
        
        cv_et(o) = cv_et(o) + pow(norm(H * xi - xi, 2), 2) / nt / pow(1- trace(H)/nt  , 2 );
      }
    }
    indexet = as<int>(whichmin(cv_sht)) - 1;
    etout = etgrid(indexet);
  }
  return(List::create(Named("sht",shtout),Named("et",etout)));
}


// cvd2search : needed in finding the dimension (new vesion). old version(cvdserach) is not used 
// [[Rcpp::export]]
double cvd2search(mat xmat, mat KY, vec shxgrid){
  Function gam("gam");
  Function matpower("matpower");
  Function normmat("normmat");
  Function whichmin("which.min");
  
  mat KX, KXi, KYi, KXii, KYii, KXiinv;
  mat normx=as<arma::mat>(normmat(xmat));
  double cvd=0,gx, shxout;
  
  // vec cv_shx(shxgrid.n_rows);
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
  KXiinv = as<mat>(matpower(KXii, -1));
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
  KXiinv = as<mat>(matpower(KXii, -1));
  cvd=cvd +  as_scalar((KYi.row(i) - KXi.row(i) * KXiinv * KYii) * ( KYi.row(i) - KXi.row(i) * KXiinv * KYii    ).t());
  }
  
  return(cvd);  
}
