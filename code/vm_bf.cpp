#include <RcppArmadillo.h>
#define ARMA_NO_DEBUG

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::NumericMatrix vm_bf(int n_points, double mu, double sd) {
  //arma::mat vm_bf(int n_points, double mu, double sd) { // definition with arma return
  
  // NumericMatrix df(n_points, 4);
  // 
  // NumericMatrix::Column pvec = df(_ , 0);  // Reference to the first column
  // NumericMatrix::Column pdf = df(_ , 1);  // Reference to the second column
  // NumericMatrix::Column basis  = df(_ , 2);  // Reference to the third column
  // NumericMatrix::Column basis_norm  = df(_ , 3);  // Reference to the fourth column
  
  arma::mat df(n_points, 4);
  
  double kappa = 1/sd;
  
  // for reasons that remain mysterious to me, using the pointer to a column with linspace generates a functioning
  // vector (e.g., the calculations and Rcout below work), but the values are not actually preserved in the resulting matrix
  //arma::vec x(df.colptr(0), n_points, false); //in-memory pointer to first column
  arma::vec x(n_points);
  arma::vec pdf(df.colptr(1), n_points, false); //in-memory pointer to second column
  arma::vec basis(df.colptr(2), n_points, false);
  arma::vec basis_norm(df.colptr(3), n_points, false);
  
  //approach in std, not as easy to convert into arma
  //std::vector<double> x(n_points);
  //std::iota(x.begin(), x.end(), 2*M_PI/(n_points - 1)); // fill x 0--2*pi

  x = linspace<vec>(0, 2*M_PI, n_points); // generate positions between 0 and 2*pi
  //Rcout << x << endl;
  
  // von Mises probability density function (pdf)
  pdf = arma::exp(kappa * arma::cos(x - mu))/(2*M_PI*R::bessel_i(kappa, 0, 1));
  
  double m = max(pdf);
  double s = sum(pdf);
  
  // AUC 1.0: pdf / sum 
  basis = pdf/s;
  
  // max 1.0: pdf / max
  basis_norm = pdf/m;
  
  //debugging attempts
  //x = pdf/5;
  //df.col(0) = pdf/5;
  //arma::vec w(df.colptr(4), n_points, false);
  //w = pdf/5;
  //df.col(4) = linspace<vec>(0, 2*M_PI, n_points);
  
  df.col(0) = x;
  
  // It's a little slower to wrap as a numeric matrix and add column names, but worth it for convenience
  Rcpp::NumericMatrix rdf(wrap(df)); // convert
  colnames(rdf) = Rcpp::CharacterVector::create("pos", "pdf", "basis", "basis_norm");
  
  return(rdf);
}


/***
vm_bf(50, 1, 0.2)

vm_mf = function(x=NULL, mu, kappa, normalize=TRUE) {
  y <- exp(kappa * cos(x - mu))/(2*pi*besselI(kappa, nu = 0))
  if (isTRUE(normalize)) {
    y <- y/max(y)
  }
  return(y)
}

get_b <- function(n_points, mu, sd) {
  x <- seq(0, 2*pi, length.out=n_points)
  vm <- vm_mf(x=x, mu=mu, kappa=1/sd, normalize=FALSE) # probability density
  avm <- vm/sum(vm) # AUC 1.0 basis (used for update eligibility)
  nvm <- vm/max(vm) # AUC 1.0 basis (used for update eligibility)
  
  basis_df <- data.frame(pvec = x, pdf = vm, basis = avm, basis_norm=nvm)
  return(basis_df)
}

#microbenchmark::microbenchmark(
# cpp= vm_bf(50, 1, 0.2),
# r= get_b(50, 1, 0.2),
# times=5000
#)
*/
