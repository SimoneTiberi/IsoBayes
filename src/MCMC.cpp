#include <cassert>
#include <ctime>
#include <iostream>
#include <random>
#include <vector>

// [[Rcpp::plugins(cpp17)]]

// Armadillo
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
//using namespace arma;


// [[Rcpp::export]]
List MCMC(Rcpp::ListOf<Rcpp::NumericVector> const& EC_numeric_multi_map,
          Rcpp::NumericVector const& Y_unique,
          Rcpp::NumericVector const& N_peptides_per_protein,
          Rcpp::NumericVector const& delta_prior,
          Rcpp::NumericVector const& PSM_multi_map,
          unsigned int const& N,
          unsigned int const& M,
          unsigned int const& K,
          unsigned int const& burn_in,
          unsigned int const& thin
) {
  Rcpp::NumericVector pi(N); // dirichlet sample.
  Rcpp::NumericVector y(N);
  
  Rcpp::NumericMatrix PI_mat( (K-burn_in)/thin, N);
  Rcpp::NumericMatrix Y_mat( (K-burn_in)/thin, N);
  
  unsigned int k, j, i, EC_len; //index, EC_len;
  
  double prob_tot, alpha;
  
  // initialize pi vector with the prior
  //for (i=0 ; i < N; i++) {
  //  pi[i] = PI_prior[i];
  //}
  // OR set al pi's to the same value:
  pi.fill(1);
  
  for (k=0 ; k<K ; k++) {
    // sample Y|pi
    // .copy or .clone
    for (i=0 ; i < N; i++) {
      y[i] = Y_unique[i];
    }
    
    // for every peptide: sample the allocation of peptides
    for (j=0 ; j<M ; j++) {
      EC_len = EC_numeric_multi_map[j].length();
      Rcpp::NumericVector pi_j(EC_len);
      Rcpp::IntegerVector y_tmp(EC_len); //, index(EC_len);
      
      for (i=0 ; i < EC_len; i++) {
        // index = EC_numeric_multi_map[j](i);
        pi_j[i] = pi[ EC_numeric_multi_map[j](i) - 1 ] / N_peptides_per_protein[ EC_numeric_multi_map[j](i) - 1 ];
      }
      prob_tot = std::accumulate(pi_j.begin(), pi_j.end(), 0.0);
      
      if( prob_tot > 0.0){
        for (i=0 ; i < EC_len; i++) {
          pi_j[i] /= prob_tot;
        }
        
        rmultinom( PSM_multi_map[j], pi_j.begin(), EC_len, y_tmp.begin());
        
        for (i=0 ; i < EC_len; i++) {
          //index = EC_numeric_multi_map[j](i);
          y[ EC_numeric_multi_map[j](i) -1 ] += y_tmp[i];
        }
      }
    }
    
    // Dirichlter sampling, pi|Y, delta
    for (i=0 ; i < N; i++) {
      alpha = y[i] + delta_prior[i];
      pi[i] = as<double>(Rcpp::rgamma(1, alpha, 1));
    }
    
    // store 1 value every "thin" iterations:
    if(k % thin == 0){
      // only keep values after burn_in:
      if(k >= burn_in){
        PI_mat( (k-burn_in)/thin, _) = pi;
        Y_mat( (k-burn_in)/thin, _) = y;
      }
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("PI") = PI_mat,
                            Rcpp::Named("Y") = Y_mat);
}

