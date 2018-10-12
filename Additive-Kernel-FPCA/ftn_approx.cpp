// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
List approx_ftn(vec x, vec t, List basis, mat R, NumericVector lambdagrid){
  Function whichmin("which.min");
  Function getbasis("getbasismatrix");
  int nt=t.size(), nbasis = R.n_cols, ncv=lambdagrid.size();
  vec mucoef=zeros(nbasis), gcv=zeros(ncv);
  double lambdaout;
  List output;
  
  mat X=as<mat>(getbasis(t, basis)); // nt X nbasis matrix
  mat B(nbasis, nbasis), Binv(nbasis, nbasis), Iden_B = eye<mat>(nbasis,nbasis);
  mat H(nt, nt), Iden_H = eye<mat>(nt, nt);
  
  vec y = x;
  
  int mincv;
  for(int j=0; j <lambdagrid.size(); j++){
    B = X.t() * X + lambdagrid(j) * R;
    Binv = inv(B);
    H = X * Binv * X.t();
    gcv(j) =as_scalar(y.t() * (Iden_H - H) *(Iden_H - H) * y )/ pow((1-trace(H)/nt),2.0);
  }
  
  mincv = as<int>(whichmin(gcv));
  B = X.t() * X + lambdagrid(mincv-1) * R;
  Binv = inv(B);
  mucoef = Binv * X.t() * y;
  lambdaout=lambdagrid(mincv-1);
  
  output = List::create(Named("mucoef",mucoef),Named("lambdaout",lambdaout));
  return(output);
}

// [[Rcpp::export]]
List approx_mat(mat x, vec t, List basis, mat R, NumericVector lambdagrid){
  // input  : x= nt X x
  int n = x.n_cols;
  
  // output
  mat mucoef=zeros(R.n_rows,n);
  vec lambdaout(n);
  List output;
  
  for(int i=0;i<n; i++){
    List tmp = approx_ftn(x.col(i), t, basis, R, lambdagrid);
    mucoef.col(i) = as<vec>(tmp[0]);
    lambdaout(i) = as<double>(tmp[1]);
  }
  
  output = List::create(Named("mucoef",mucoef),Named("lambdaout",lambdaout));
  return(output);
}



// [[Rcpp::export]]
List approx_array(vec x, vec t, List basis, mat R, NumericVector lambdagrid){
  // vectorized X. (nt * p) X n matrix
  int n = x.n_cols, nt = t.size(), nbasis = R.n_cols;
  int p = x.n_rows / nt;
  mat mucoef(nbasis*p, n);
  mat lambdaout(n, p);
  for(int i=0; i<p;i++){
    List tmp = approx_mat(x.submat(i*nt,0,(i+1)*nt-1,n-1), t, basis, R, lambdagrid); //   n * nt matrix   . i-th column.
    mucoef.submat(i*nbasis,0,(i+1)*nbasis-1,n) = as<mat>(tmp[0]);
    lambdaout.col(i) = as<vec>(tmp[1]);
  }
  List output = List::create(Named("mucoef",mucoef),Named("lambdaout",lambdaout));
  return(output);
}


// [[Rcpp::export]]
List approx_list(List x, List t, List basis, mat R, NumericVector lambdagrid){
  int n = x.size();
  
  mat mucoef=zeros(R.n_rows,n);
  vec lambdaout(n);
  List output;
  
  for(int i=0;i<n; i++){
    
    List tmp = approx_ftn(x[i], t[i], basis, R, lambdagrid);
    mucoef.col(i) = as<vec>(tmp[0]);
    lambdaout(i) = as<double>(tmp[1]);
  }
  
  output = List::create(Named("mucoef",mucoef),Named("lambdaout",lambdaout));
  return(output);
}
// n* nt * p    n nt * p 
// [[Rcpp::export]]
mat testfun(cube x){
  x.reshape(x.n_cols*x.n_rows,x.n_slices,1);
  mat A=x.slice(0);
  return(A);
}

// [[Rcpp::export]]
vec testfun2(List x){
  return(x[0]);
}