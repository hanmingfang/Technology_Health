#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//Parameters
double phi = 1.0;
double xi = 0.455;
double util_min = 0.001;
double rate_g = 2.0;

// [[Rcpp::export]]
double h_g(double m) {
  return rate_g*exp(-rate_g*m);//exponential density
}

// [[Rcpp::export]]
double l0_s_cpp(double w0) {
  return pow(w0/phi,1/xi);
}

// [[Rcpp::export]]
double u0_cpp(double theta, double w0, double m) {
  double l0 = l0_s_cpp(w0);
  return (1/(1-theta))*(pow(std::max(w0*l0 - m - phi*(pow(l0,1+xi)/(1+xi)),util_min),1-theta)-1);
}

// [[Rcpp::export]]
NumericVector integrand_u0(double theta, double w0, NumericVector m){
  NumericVector aux(m.size());
  //draw vector of random numbers with exponential distribution
  int n = m.size();
  for(int i = 0; i < n; i++){
    aux[i] = u0_cpp(theta,w0,m[i])*h_g(m[i]);    
  }
  return aux;
}


/*** R
integrand_u0(theta = 0.5, w0 = 2.4, m = c(1,2,3,4,5))
*/
