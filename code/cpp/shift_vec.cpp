#include <Rcpp.h>

using namespace Rcpp;

// fast circular shift of vector. Approx 5x than pure R method
// note that this always returns a numeric result, even if integers are passed

// [[Rcpp::export]]

std::vector<double> shift_vec(std::vector<double> x, int delta)
{
  int N = x.size();
  
  // wrap around if shift is larger than length
  if (abs(delta) >= N) delta = delta % N;
  
  if (delta > 0)
    std::rotate(x.begin(), x.begin() + delta, x.end());
  else if (delta < 0)
    std::rotate(x.rbegin(), x.rbegin() + abs(delta), x.rend());
  
  return(x);
}

   
/***
# comparison pure R function
 shift_vec_r = function(x, delta) {
   N <- length(x)
   delta <- delta %% N # wrap around if shift is larger than length
   if (delta == 0L) x else c(tail(x, -delta), head(x, delta))
 }
 
 b1 <- shift_vec_r(1:100, 202)
 print(b1)
 b2 <- shift_vec(1:100, 202)
 print(b2)
 identical(b1, b2)
 microbenchmark::microbenchmark(
  cpp = shift_vec(1:100, -99),
  r = shift_vec_r(1:100, -99),
  times=50000
 )

 // Unit: microseconds
 //expr    min     lq      mean median     uq      max neval
 // cpp  1.633  2.210  3.150724  2.661  3.004 10611.47 50000
 //   r 10.572 12.899 15.676957 13.313 13.982 11536.16 50000
 
*/
