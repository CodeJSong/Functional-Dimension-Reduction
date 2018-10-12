// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

/*
##################################################################
# largest_eval(A) : Find the largest eigenvalue of the symmetric matrix A --
##################################################################
# input : A : n X n symmetric matrix
# value : largest e_val of A
##################################################################
 */
// [[Rcpp::export]]
double largest_eval(mat M){
  vec eval;
  eig_sym(eval, M);
  double largest_eval = eval.max();
  return(largest_eval);
}


/*
##################################################################
# bDiag(A,d) : blog diagonal matrix --
##################################################################
# input : A : n X p matrix
# value : n*d X p*d block diagonal matrix filled with A
##################################################################
 */
// [[Rcpp::export]]
mat bDiag(mat A, int d){
  mat B(d*A.n_rows,d*A.n_cols,fill::zeros);
  for(int i=0; i<d;i++){
    B.submat(i*A.n_rows, i*A.n_cols, (i+1)*A.n_rows-1, (i+1)*A.n_cols-1) = A;
  }
  return(B);
}

/*
##################################################################
# normmat(x) ------
##################################################################
# Input - x is matrix (n x p)  
# Output - n x n matrix  out_ij = \| x_i - x_j \|^2
##################################################################
*/
// [[Rcpp::export]]
mat normmat(mat A){
  mat aa = A * A.t();
  mat daa(A.n_rows,A.n_rows); daa.each_col() = diagvec(aa);
  mat normA = (-2 * aa + daa).t() + daa;
  normA.elem(find(normA < pow(10,-16) )).zeros();
  return(normA);
}



/* 
##################################################################
# gam(normX) estimating width-parameter of kernel function -----
##################################################################
# input : normX : n x n matrix such that 
#       (normX)_ij = ||x_i - x_j||^2
# value : estimated gamma 
##################################################################
*/

// [[Rcpp::export]]
double gam(mat normx){
  int n = normx.n_rows;
  double avg;
  avg = accu(sqrt(normx));  // sum of all elements  = 2x sum_{i>j} ||x_i-x_j||^2
  avg = avg / (n * (n-1));  //  sum_{i>j} ||x_i-x_j||^2  / (n choose 2)
  return(1/(2*pow(avg,2)));
}



/*
##################################################################
# mppower(x,a) : power of matrix (moore penrose type)------
##################################################################
# input : x : matrix
# value : x^a ignoring the eigenvectors with eigenvalue < 10^(-15)
##################################################################
*/
// [[Rcpp::export]]
mat mppower(mat A, double alpha){
  double ignore=pow(10,-15);
  mat B=(A.t() + A) / 2;
  mat A_evec;
  vec A_eval;
  eig_sym(A_eval, A_evec, B);
  uvec indices = find((A_eval) > ignore);
  return(A_evec.cols(indices) * diagmat(pow(A_eval(indices),alpha)) * A_evec.cols(indices).t()); 
}

/*
##################################################################
# matpower2(x,a,epsilon) : power of matrix(regularized moore-penrose)------
##################################################################
# input : A (matrix), alpha(power), epsilon(sale)
# value : (A+ lambda_1 epsilon I)^+{alpha}
# lambda_1 : the largest eigenvalue of A
##################################################################
*/

// [[Rcpp::export]]

mat matpower2(mat A, double alpha, double epsilon){
  double ignore=pow(10,-15);
  mat B=(A.t() + A) / 2;
  mat A_evec, iden(A.n_rows,A.n_cols); iden.eye();
  vec A_eval, A_eval2;
  eig_sym(A_eval2, B);
  epsilon = epsilon * A_eval2(A_eval2.n_rows-1);   // largest eigenvalue of the matrix * epsilon
  eig_sym(A_eval, A_evec, B + epsilon *  iden);
  uvec indices = find((A_eval) > ignore);
  return(A_evec.cols(indices) * diagmat(pow(A_eval(indices),alpha)) * A_evec.cols(indices).t());
}


/*
##################################################################
# xiplist :  inner product matrix when x is vector-valued function------
##################################################################
# Inner product matrix of X
# output : n x n matrix,  (x.ip)ij= <xi,xj>_{H}
# x : list of x
# KT : gram matrix of T  which is used in the inner product
###############################################################
*/

// [[Rcpp::export]]

mat xiplist(int n, List x, List tt, vec tuni, mat KT, double et, double teq, int dim){
  Function makematrix("makematrix");
  Function matpower2("matpower2");
  if(teq==1){
    mat Vinv = as<mat>(matpower2(KT, -1, et));
    if(dim==1){ 
      mat xmat = as<mat>(makematrix(wrap(x),tuni.n_rows));
      return(xmat * Vinv * KT * Vinv * xmat.t());
    }
    if(dim>1){
      mat xxip(n,n); xxip.zeros();
      for(int d=0;d<dim;d++){
        mat xmat(n,KT.n_rows);
        for(int i=0;i<n;i++){
        xmat.row(i) = as<mat>(x[i]).col(d).t();
        }
        xxip = xxip + xmat * Vinv * KT * Vinv * xmat.t();
      }
    return(xxip);
    }
  }
  if(teq==0){
    Function which("which");
    mat xxip(n,n); xxip.zeros();
    if(dim > 1){
      for(int i=0;i<(n-1);i++){
        vec ti = as<vec>(tt[i]);
        int ni=ti.n_rows;
        vec index(ni);
        
        for(int k=0;k<ni;k++){
          index(k) = as<int>(wrap(find(tuni==ti(k))));
        }
        uvec indexi = conv_to<uvec>::from(index);
        mat Vii = KT.submat((indexi),(indexi));
        mat Viiinv = as<mat>(matpower2(Vii, -1, et));
        mat xi=as<mat>(x[i]);
        for(int d=0;d<dim;d++){
          xxip(i,i) = xxip(i,i) +  as_scalar(xi.col(d).t() * Viiinv * Vii * Viiinv * xi.col(d));
        }
        for(int j=(i+1);j<n;j++){
          vec tj = as<vec>(tt[j]);
          int nj = tj.n_rows;
          vec index(nj);
          for(int k=0;k<nj;k++){
            index(k) =  as<int>(wrap(find(tuni==tj(k))));
          }
          uvec indexj = conv_to<uvec>::from(index);
          mat Vjj = KT.submat(indexj,indexj);
          mat Vjjinv = as<mat>(matpower2(Vjj, -1, et));
          mat Vij = KT.submat(indexi,indexj);
          mat xj=as<mat>(x[j]);
          
          for(int d=0;d<dim;d++){
            xxip(i,j) = xxip(i,j) +  as_scalar(xi.col(d).t() * Viiinv * Vij * Vjjinv * xj.col(d));
            
          }
        }
      }
      int i=n-1;
      vec ti = as<vec>(tt[i]);
      int ni = ti.n_rows;
      vec index(ni);
      for(int k=0;k<ni;k++){
        index(k) =  as<int>(wrap(find(tuni==ti(k))));
      }
      uvec indexi = conv_to<uvec>::from(index);
      mat Vii = KT.submat(indexi,indexi);
      mat Viiinv = as<mat>(matpower2(Vii, -1, et));
      mat xi=as<mat>(x[i]);
      
      
      for(int d=0;d<dim;d++){
        xxip(i,i) = xxip(i,i) +  as_scalar(xi.col(d).t() * Viiinv * Vii * Viiinv * xi.col(d));
      }
      
    }
    if(dim==1){
      for(int i=0;i<(n-1);i++){
        vec ti = as<vec>(tt[i]);
        int ni=ti.n_rows;
        vec index(ni);
        
        for(int k=0;k<ni;k++){
          index(k) = as<int>(wrap(find(tuni==ti(k))));
        }
        uvec indexi = conv_to<uvec>::from(index);
        mat Vii = KT.submat((indexi),(indexi));
        mat Viiinv = as<mat>(matpower2(Vii, -1, et));
        vec xi=as<vec>(x[i]);
        xxip(i,i) = as_scalar(xi.t() * Viiinv * Vii * Viiinv * xi);
        
        for(int j=(i+1);j<n;j++){
          vec tj = as<vec>(tt[j]);
          int nj = tj.n_rows;
          vec index(nj);
          for(int k=0;k<nj;k++){
            index(k) =  as<int>(wrap(find(tuni==tj(k))));
          }
          uvec indexj = conv_to<uvec>::from(index);
          mat Vjj = KT.submat(indexj,indexj);
          mat Vjjinv = as<mat>(matpower2(Vjj, -1, et));
          mat Vij = KT.submat(indexi,indexj);
          vec xj=as<vec>(x[j]);
          xxip(i,j) =  as_scalar(xi.t() * Viiinv * Vij * Vjjinv * xj);
        }
      }
      int i=n-1;
      vec ti = as<vec>(tt[i]);
      int ni = ti.n_rows;
      vec index(ni);
      for(int k=0;k<ni;k++){
        index(k) =  as<int>(wrap(find(tuni==ti(k))));
      }
      uvec indexi = conv_to<uvec>::from(index);
      mat Vii = KT.submat(indexi,indexi);
      mat Viiinv = as<mat>(matpower2(Vii, -1, et));
      vec xi=as<vec>(x[i]);
      
      xxip(i,i) = as_scalar(xi.t() * Viiinv * Vii * Viiinv * xi);
      
      if(dim>1){
        for(int d=0;d<dim;d++){
          xxip(i,i) = xxip(i,i) +  as_scalar(xi.col(d).t() * Viiinv * Vii * Viiinv * xi.col(d));
        }
      }
    }
    return(xxip + xxip.t() - diagmat(xxip));  
  }
  return(0);
}