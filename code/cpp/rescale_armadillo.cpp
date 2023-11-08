#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::rowvec rescale_armadillo(arma::rowvec x, double to_l, double to_h)
{
  if (abs(to_l - to_h) < 1e-5) stop("No range in to values");
  //NumericVector x2 = na_omit(x);
  //double from_l = min(x2);
  //double from_h = max(x2);
  double from_l = min(x);
  double from_h = max(x);
  //Rcout << "Min is: " << from_l << endl;
  //Rcout << "Max is: " << from_h << endl;
  x = (x - from_l)/(from_h - from_l) * (to_h - to_l) + to_l;
  
  //return(wrap(x));
  return(x);
}

   
/*** R
# comparison function
 v <- c(1:10000)
 
 b1 <- scales::rescale(v, to=c(5,15))
 b2 <- rescale_armadillo(v, 5, 15)
 identical(b1, b2)
 
 #Approx 3x speedup
 microbenchmark::microbenchmark(
  cpp = rescale_armadillo(v, 5, 15),
  r = scales::rescale(v, to=c(5,15)),
  times=50000
 )

*/
