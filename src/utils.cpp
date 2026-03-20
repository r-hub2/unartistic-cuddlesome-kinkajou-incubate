#include <cpp11.hpp>
#include <Rmath.h>

using namespace cpp11;

//' Difference in Log-Space
//'
//' difference of two values in log-space
//' @param lx log of first value
//' @param ly log of second value
//' @return difference in log-space: `log(exp(lx)-exp(ly))`
//' @export
[[cpp11::register]]
double logspace_sub_cpp(double lx, double ly) {
  // use logspace_sub from Rmath
  return logspace_sub(lx, ly);
}

//' Difference in log-space
//'
//' log-space difference of first two values in a vector
//' @param lxy vector with two values
//' @return difference in log-space: `log(exp(lxy[2nd]) - exp(lxy[1st]))`
//' @export
[[cpp11::register]]
double logspace_sub2_cpp(doubles lxy) {
  // use logspace_sub from Rmath
  return logspace_sub(lxy[1], lxy[0]);
}
