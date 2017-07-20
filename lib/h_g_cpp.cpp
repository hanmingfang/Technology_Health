#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;
using namespace arma;
// Enable C++11
// [[Rcpp::plugins(cpp11)]]
// Enable OpenMP (excludes macOS)
// [[Rcpp::plugins(openmp)]]

//Parameters
double phi = 1.0;
double xi = 0.455;
double util_min = 0.001;
double rate_g = 2.0;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

//This is just to try how to generate a random vector of Gamma dist in Armadillo
arma::vec h_g_rng(int n_elem, double a, double b) {
  return arma::randg<arma::vec>(n_elem, distr_param(a,b));
}

/*** R
h_g_rng(10,1.0,1.0)
*/
