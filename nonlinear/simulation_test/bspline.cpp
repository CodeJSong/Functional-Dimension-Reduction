// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]

SEXP evalbasis(NumericVector t, List basisobj){
  using namespace arma;
  Function splinedesign("splineDesign");
  List basis(basisobj);
  vec breaks=as<vec>(basis["params"]);
  CharacterVector type=as<CharacterVector>(basis["type"]);
  vec argval=as<vec>(basis["rangeval"]), tt=as<vec>(t);
  mat a1(1,1), a2(1,1);
  a1[0]=argval[0];
  a2[0]=argval[1];
  
    double nderiv=0,  m=tt.size(), norder=4;
    mat breakmat=mat(breaks);
    breakmat.insert_rows(breaks.size(),a2);
    breakmat.insert_rows(0,a1);
    breaks = vectorise(breakmat);
    
    NumericVector front, back;
    front=rep(breaks[0], norder-1);
    back=rep(breaks[breaks.size()-1],norder-1);
    
    breakmat.insert_rows(breaks.size(),mat(back));
    breakmat.insert_rows(0,mat(front));
    NumericVector derivs = rep(nderiv,m);
    vec deriv=as<vec>(derivs), knots=vectorise(breakmat);
  
  return(splinedesign(knots, tt, norder, deriv));
}

// [[Rcpp::export]]
List bspline_mat(mat x, vec t, List basis, mat R, NumericVector lambdagrid, int gcvall){
  Function whichmin("which.min");
  Function getbasis("getbasismatrix");
  // input
  int n = x.n_cols, nt = t.size();
  
  // output
  mat gcv=zeros(lambdagrid.size(),n), mucoef=zeros(R.n_rows,n);
  vec lambdaout(n);
  List output;
  
  
  mat X=as<mat>(getbasis(t, basis)); // nt X nbasis matrix
  for(int i=0;i<n; i++){
    mat B(X.n_cols, X.n_cols), Binv(X.n_cols, X.n_cols), Iden_B = eye<mat>(X.n_cols,X.n_cols);
    mat H(X.n_rows, X.n_rows), Iden_H = eye<mat>(X.n_rows, X.n_rows);
    
    vec y = x.col(i);
    
    int mincv, nt=X.n_rows;
    for(int j=0; j <lambdagrid.size(); j++){
      B = X.t() * X + lambdagrid(j) * R;
      Binv = inv(B);
      H = X * Binv * X.t();
      gcv(j,i) =as_scalar(y.t() * (Iden_H - H) *(Iden_H - H) * y )/ pow((1-trace(H)/nt),2.0);
    }
    
    vec gcvcols=gcv.col(i);
    mincv = as<int>(whichmin(gcvcols));
    B = X.t() * X + lambdagrid(mincv-1) * R;
    Binv = inv(B);
    mucoef.col(i) = Binv * X.t() * y;
    lambdaout[i]=lambdagrid(mincv-1);
  }
  
  output = List::create(Named("mucoef",mucoef),Named("lambdaout",lambdaout));
  return(output);
}
// [[Rcpp::export]]

List bspline_unbal(List x, List t, List basisobj, mat RR, NumericVector lambdagrid, int gcvall){
  Function whichmin("which.min");
  Function getbasis("getbasismatrix");
  List xx(x), tt(t), basis(basisobj), output;
  int n=x.size();
  mat X,R(RR);
  NumericVector lambdagrd = as<NumericVector>(lambdagrid);
  mat gcv=zeros<mat>(lambdagrd.size(),n), mucoef=zeros<mat>(R.n_rows,n);
  vec lambdaout(n);
  if(gcvall==0){
    for(int i=0; i<n; i++){
      //X=as<mat>(evalbasis(tt[i], basis));
      X=as<mat>(getbasis(tt[i], basis));
      mat B(X.n_cols, X.n_cols), Binv(X.n_cols, X.n_cols), Iden_B = eye<mat>(X.n_cols,X.n_cols);
      mat H(X.n_rows, X.n_rows), Iden_H = eye<mat>(X.n_rows, X.n_rows);
    
      vec y = as<vec>(xx[i]);
    
      int mincv, nt=X.n_rows;
      for(int j=0; j <lambdagrd.size(); j++){
        B = X.t() * X + lambdagrd(j) * R;
        Binv = inv(B);
        H = X * Binv * X.t();
        gcv(j,i) =as_scalar(y.t() * (Iden_H - H) *(Iden_H - H) * y )/ pow((1-trace(H)/nt),2.0);
      }
    
      vec gcvcols=gcv.col(i);
      mincv = as<int>(whichmin(gcvcols));
      B = X.t() * X + lambdagrd(mincv-1) * R;
      Binv = inv(B);
      mucoef.col(i) = Binv * X.t() * y;
      lambdaout[i]=lambdagrd(mincv-1);
    }
  }else if(gcvall==1){
    vec gcvvec(lambdagrd.size());
    for(int i=0; i<n; i++){
      //X=as<mat>(evalbasis(tt[i], basis));
      X=as<mat>(getbasis(tt[i], basis));
      mat B(X.n_cols, X.n_cols), Binv(X.n_cols, X.n_cols), Iden_B = eye<mat>(X.n_cols,X.n_cols);
      mat H(X.n_rows, X.n_rows), Iden_H = eye<mat>(X.n_rows, X.n_rows);
    
      vec y = as<vec>(xx[i]);
    
      int nt=X.n_rows;
      for(int j=0; j <lambdagrd.size(); j++){
        B = X.t() * X + lambdagrd(j) * R;
        Binv = inv(B);
        H = X * Binv * X.t();
        gcvvec(j) = gcvvec(j) + as_scalar(y.t() * (Iden_H - H) *(Iden_H - H) * y )/ pow((1-trace(H)/nt),2.0);
      }
    }
    int mincv = as<int>(whichmin(gcvvec));
    for(int i=0;i<n;i++){
      // X=as<mat>(evalbasis(tt[i], basis));
      X=as<mat>(getbasis(tt[i], basis));
      vec y = as<vec>(xx[i]);
      mat B = X.t() * X + lambdagrd(mincv-1) * R;
      mat Binv = inv(B);
      mucoef.col(i) = Binv * X.t() * y;
      lambdaout[i]=lambdagrd(mincv-1);
    }
  }
  
  output = List::create(Named("mucoef",mucoef),Named("lambdaout",lambdaout));
  return(output);
}


// [[Rcpp::export]]


List bsplinegcv(mat xmat, mat rr, vec ystar, NumericVector lambdagrid){
  mat X(xmat), R(rr);
  mat B(X.n_cols, X.n_cols), Binv(X.n_cols, X.n_cols), Iden_B = eye<mat>(X.n_cols,X.n_cols);
  mat H(X.n_rows, X.n_rows), Iden_H = eye<mat>(X.n_rows, X.n_rows);
  vec y(ystar);
  NumericVector lambdagrd(lambdagrid);
  vec gcv(lambdagrd.size()), mucoef(X.n_cols);
  int mincv;
  
  Function whichmin("which.min");
  for(int i=0; i <lambdagrd.size(); i++){
    B = X.t() * X + lambdagrd(i) * R;
    Binv = inv(B);
    H = X * Binv * X.t();
    gcv[i] = as_scalar(y.t() * (Iden_H - H) *(Iden_H - H) * y )/ pow((1-trace(H)/X.n_rows),2.0);
  }
  mincv = as<int>(whichmin(gcv));
  B = X.t() * X + lambdagrd(mincv-1) * R;
  Binv = inv(B);
  mucoef = Binv * X.t() * y;
  return(List::create(Named("gcv",gcv), Named("mucoef",mucoef),Named("lambdas",lambdagrd),
  Named("lambda",lambdagrd(mincv-1))));
}
  

// [[Rcpp::export]]

List bsplinegcv_mat(mat X, mat R, mat y, NumericVector lambdagrd){
  mat mucoef;
  mat B(X.n_cols, X.n_cols), Binv(X.n_cols, X.n_cols), Iden_B = eye<mat>(X.n_cols,X.n_cols);
  mat H(X.n_rows, X.n_rows), Iden_H = eye<mat>(X.n_rows, X.n_rows);
  vec gcv(lambdagrd.size());
  int mincv;
  
  Function whichmin("which.min");
  
  for(int i=0; i <lambdagrd.size(); i++){
    B = X.t() * X + lambdagrd(i) * R;
    Binv = inv(B);
    H = X * Binv * X.t();
    gcv[i] =trace(y.t() * (Iden_H - H) *(Iden_H - H) * y )/ pow((1-trace(H)/X.n_rows),2.0);
  }
  mincv = as<int>(whichmin(gcv));
  B = X.t() * X + lambdagrd(mincv-1) * R;
  Binv = inv(B);
  mucoef = Binv * X.t() * y;
  
  return(List::create(Named("gcv",gcv), Named("mucoef",mucoef),Named("lambdas",lambdagrd),
  Named("lambda",lambdagrd(mincv-1))));
}


// [[Rcpp::export]]

List bspline_pca(List xx, List tt, List basis, List pcaharmonics, mat R, NumericVector lambdagrid, int gcvall){
  // Using PCA basis 
  Function whichmin("which.min");
  Function getbasis("getbasismatrix");
  List output;
  int n=xx.size();
  mat X, pcacoef = as<mat>(pcaharmonics["coefs"]);
  NumericVector lambdagrd = as<NumericVector>(lambdagrid);
  mat gcv=zeros<mat>(lambdagrd.size(),n), mucoef=zeros<mat>(R.n_rows,n);
  vec lambdaout(n);
  if(gcvall==0){
    for(int i=0; i<n; i++){
      //X=as<mat>(evalbasis(tt[i], basis));
      X=as<mat>(getbasis(tt[i], basis)) * pcacoef;
      mat B(X.n_cols, X.n_cols), Binv(X.n_cols, X.n_cols), Iden_B = eye<mat>(X.n_cols,X.n_cols);
      mat H(X.n_rows, X.n_rows), Iden_H = eye<mat>(X.n_rows, X.n_rows);
      
      vec y = as<vec>(xx[i]);
      
      int mincv, nt=X.n_rows;
      for(int j=0; j <lambdagrd.size(); j++){
        B = X.t() * X + lambdagrd(j) * R;
        Binv = inv(B);
        H = X * Binv * X.t();
        gcv(j,i) =as_scalar(y.t() * (Iden_H - H) *(Iden_H - H) * y )/ pow((1-trace(H)/nt),2.0);
      }
      
      vec gcvcols=gcv.col(i);
      mincv = as<int>(whichmin(gcvcols));
      B = X.t() * X + lambdagrd(mincv-1) * R;
      Binv = inv(B);
      mucoef.col(i) = Binv * X.t() * y;
      lambdaout[i]=lambdagrd(mincv-1);
    }
  }else if(gcvall==1){
    vec gcvvec(lambdagrd.size());
    for(int i=0; i<n; i++){
      //X=as<mat>(evalbasis(tt[i], basis));
      X=as<mat>(getbasis(tt[i], basis)) * pcacoef;
      mat B(X.n_cols, X.n_cols), Binv(X.n_cols, X.n_cols), Iden_B = eye<mat>(X.n_cols,X.n_cols);
      mat H(X.n_rows, X.n_rows), Iden_H = eye<mat>(X.n_rows, X.n_rows);
      
      vec y = as<vec>(xx[i]);
      
      int nt=X.n_rows;
      for(int j=0; j <lambdagrd.size(); j++){
        B = X.t() * X + lambdagrd(j) * R;
        Binv = inv(B);
        H = X * Binv * X.t();
        gcvvec(j) = gcvvec(j) + as_scalar(y.t() * (Iden_H - H) *(Iden_H - H) * y )/ pow((1-trace(H)/nt),2.0);
      }
    }
    int mincv = as<int>(whichmin(gcvvec));
    for(int i=0;i<n;i++){
      // X=as<mat>(evalbasis(tt[i], basis));
      X=as<mat>(getbasis(tt[i], basis)) * pcacoef;
      vec y = as<vec>(xx[i]);
      mat B = X.t() * X + lambdagrd(mincv-1) * R;
      mat Binv = inv(B);
      mucoef.col(i) = Binv * X.t() * y;
      lambdaout[i]=lambdagrd(mincv-1);
    }
  }
  
  output = List::create(Named("mucoef",mucoef),Named("lambdaout",lambdaout));
  return(output);
}