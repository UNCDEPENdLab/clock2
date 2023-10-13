#include <Rcpp.h>

using namespace Rcpp;

// somewhat faster version of scales::rescale (~2x)

// [[Rcpp::export]]

NumericVector rescale(NumericVector x, double to_l, double to_h)
{
  // NumericVector v  {10,20,30,40,50};
  // double x1 = v[0];
  // Annoying problems with trying to get scalar out of Nullable NumericVector
  // //double fl;
  // double fh;
  // //if (from_l.isNotNull())
  //   //double fl = Rcpp::as<double>(from_l[0]);
  //   double fl = from_l[0];
  // else
  //   fl = min(x);
  // 
  // if (from_h.isNotNull())
  //   fh = from_h[0];
  // else
  //   fh = max(x);
  // 
  // 
  
  if (abs(to_l - to_h) < 1e-5) stop("No range in to values");
  NumericVector x2 = na_omit(x);
  double from_l = min(x2);
  double from_h = max(x2);
  
  x = (x - from_l)/(from_h - from_l) * (to_h - to_l) + to_l;
  
  return(x);
}

   
/*** R
# comparison function
 v <- c(NA, 1:100000)
 
rescale(v, 5, 15)

 b1 <- scales::rescale(v, to=c(5,15))
 b2 <- rescale(v, 5, 15)
 identical(b1, b2)
 microbenchmark::microbenchmark(
  cpp = rescale(v, 5, 15),
  r = scales::rescale(v, to=c(5,15)),
  times=10000
 )

*/
