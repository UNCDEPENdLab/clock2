#include <Rcpp.h>

using namespace Rcpp;

// description rescale a vector to have a given 2*spread range around the mean and a designated mean
// note that this doesn't handle NAs elegantly!

// [[Rcpp::export]]

NumericVector v_rescale(NumericVector x, double spread = 20, double mean_val = 50, double force_min = 1)
{
  
  double to_l = mean_val - spread;
  double to_h = mean_val + spread;
  //if (is.na(l) || is.na(h)) warning("Missing rescale values")
  if (abs(to_l - to_h) < 1e-3) warning("Low and high rescale values do not differ");
  
  if (abs(to_l - to_h) < 1e-5) stop("No range in to values");
  //NumericVector x2 = na_omit(x);
  //double from_l = min(x2);
  //double from_h = max(x2);
  double from_l = min(x);
  double from_h = max(x);
  
  // doesn't handle NAs elegantly  
  x = (x - from_l)/(from_h - from_l) * (to_h - to_l) + to_l;
  
  double adjust = mean_val - mean(x);
  x = x + adjust;
  
  // if we want to force a given minimum to hold, apply it here (can undermine the mean)
  if (force_min > 0) x = x - (min(x) - force_min);
  return(x);
}

/*** 
  # comparison function
  v_rescale_r <- function(x, spread = 20, mean_val=50, force_min = 1) {
    checkmate::assert_number(force_min, lower=0)
    # mx <- mean(x)
    l <- mean_val - spread
    h <- mean_val + spread
    if (is.na(l) || is.na(h)) warning("Missing rescale values")
    if (abs(l - h) < 1e-3) warning("Low and high rescale values do not differ")
    y <- scales::rescale(x, to = c(l, h))
    adjust <- mean_val - mean(y)
    y <- y + adjust
    
    # if we want to force a given minimum to hold, apply it here (can undermine the mean)
    if (force_min > 0) y <- y - (min(y) - force_min)
    return(y)
  }
  
  v <- rnorm(100, mean=10, sd=2)
   
 b1 <- v_rescale_r(v, spread=50, mean_val=100, force_min=2)
 b2 <- v_rescale(v, spread=50, mean_val=100, force_min=2)
 identical(b1, b2)
 microbenchmark::microbenchmark(
   cpp = v_rescale(v, spread=50, mean_val=100, force_min=2),
   r = v_rescale_r(v, spread=50, mean_val=100, force_min=2),
   times=10000
 )

 # about 10x faster
 # Unit: microseconds
 # expr    min     lq      mean  median      uq      max neval
 # cpp  2.508  3.215  5.071137  3.6575  4.3210 9925.332 10000
 # r 30.192 37.021 41.865751 38.1495 40.6185 9314.148 10000

*/
