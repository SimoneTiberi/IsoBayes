#include <cassert>
#include <ctime>
#include <iostream>
#include <random>
#include <vector>

// [[Rcpp::plugins(cpp11)]]

// Armadillo
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
//using namespace arma;


// [[Rcpp::export]]
List MCMC_Unique(Rcpp::NumericVector const& Y_unique,
          Rcpp::NumericVector const& delta_prior,
          unsigned int const& N,
          unsigned int const& K,
          unsigned int const& burn_in,
          unsigned int const& thin
) {
  Rcpp::NumericVector pi(N); // dirichlet sample.

  Rcpp::NumericMatrix PI_mat( (K-burn_in)/thin, N);

  unsigned int k, i;
  
  double dir_sample_sum, alpha;
  
  pi.fill(1);
  
  for (k=0 ; k<K ; k++) {
    // sample Y|pi
    // .copy or .clone

    // Dirichlter sampling, pi|Y, delta
    dir_sample_sum = 0.0;
    
    for (i=0 ; i < N; i++) {
      alpha = Y_unique[i] + delta_prior[i];
      pi[i] = as<double>(Rcpp::rgamma(1, alpha, 1));
      // test:
      // pi[i] = R::rgamma( alpha, 1);
      dir_sample_sum += pi[i];
    }
    
    for (i=0 ; i < N; i++) {
      pi[i] /= dir_sample_sum;
    }
    
    // store 1 value every "thin" iterations:
    if(k % thin == 0){
      // only keep values after burn_in:
      if(k >= burn_in){
        PI_mat( (k-burn_in)/thin, _) = pi;
      }
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("PI") = PI_mat);
}
