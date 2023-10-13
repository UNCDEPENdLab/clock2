#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]

double compute_overlap(NumericVector b1, NumericVector b2) {
  double s1 = sum(b1);
  if (std::abs(s1) - 1.0 > 1e-5) 
    stop("b1 does not sum to 1");
  
  double s2 = sum(b2);
  if (std::abs(s2) - 1.0 > 1e-5) 
    stop("b2 does not sum to 1");
  
  double res = sum(pmin(b1, b2));
  return(res);
}

/***
# b1 <- vm_cpp(100, 0.5, 0.2)[,"basis"]
# b2 <- vm_cpp(100, 1.5, 0.2)[,"basis"]

# function that computes the proportion overlap between two vm bfs
# compute_overlap_r <- function(b1, b2) {
#   stopifnot(abs(sum(b1) - 1) < 1e-5) # verify AUC = 1
#   stopifnot(abs(sum(b2) - 1) < 1e-5)
#   sum(pmin(b1, b2))
# }
# 
# compute_overlap(b1, b2)
# compute_overlap_r(b1, b2)
# 
# microbenchmark::microbenchmark(
#  cpp = compute_overlap(b1, b2),
#  r = compute_overlap_r(b1, b2),
#  times=50000
# )

#6x improvement
# Unit: microseconds
# expr   min     lq      mean median     uq       max neval cld
# cpp  1.568  2.050  2.429526  2.219  2.586  1238.723 50000  a 
# r    9.996 11.446 14.677729 13.257 14.036 18013.581 50000   b

*/
