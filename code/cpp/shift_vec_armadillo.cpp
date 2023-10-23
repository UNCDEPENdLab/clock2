#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// note that this is slower than the std::rotate approach

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec shift_vec_arma(arma::vec x, int delta)
{
  int N = x.size();
  if (abs(delta) >= N) delta = delta % N;
  
  return(shift(x, delta));
}

   
/*** R
# comparison function
 shift_vec_r = function(x, delta) {
   N <- length(x)
   delta <- delta %% N # wrap around if shift is larger than length
   if (delta == 0L) x else c(tail(x, -delta), head(x, delta))
 }
 
 b1 <- shift_vec_r(1:100, 202)
 print(b1)
 b2 <- shift_vec_arma(1:100, 202)
 print(b2)

 microbenchmark::microbenchmark(
  cpp = shift_vec_arma(1:100, -99),
  r = shift_vec_r(1:100, -99),
  times=50000
 )
    
*/
