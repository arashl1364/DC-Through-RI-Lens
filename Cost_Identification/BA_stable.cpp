// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
using namespace arma; // use the Armadillo library for matrix computations
using namespace Rcpp;

// [[Rcpp::export]]
List BA_upd(mat X, vec beta, double lambda, vec mu, mat states, int nobs, double betacmplx){
    // Blahut Arimoto algorithm (design: inside brands + outside good + price + bonus)
    // X here is only the design with simple attributes
    // states is the matrix of all the realization of states already multiplied by beta(complex)
    // mu is the prior
    // lambda is the cost of information acquisition
  int nstates = states.n_cols;
  double nalt = states.n_rows; 
    //int nvar = beta.n_elem;

  vec oldProbs = zeros<vec>(nalt);
  vec cProbs = oldProbs;
  mat oldCondProbs = zeros<mat>(nalt,nstates); 
  mat cCondProbs = oldCondProbs;
  mat post = zeros<mat>(nalt,nstates);
  
  // Calculate Omega (states) | beta
  ////////
  mat Xbeta = X*beta; //(span(0,nvar-1 - 1)); //last beta is the costly attribute (bonus)
  ////////
  mat Omega = zeros<mat>(nalt,nstates);
  for(int j=0;j<nstates;j++) Omega.col(j) = Xbeta;
  Omega = Omega + betacmplx*states; //beta(nvar-1) * states;  
  
  // Calculate P(v) | Omega, mu, lambda  using BA algorithm
  // max_iter: maximal number of iterations 
  // precision: threshold for the sum of squared differences between the choice 
  // probabilities of two subsequent iterations so that iterations stop
  int max_iter = 100000000;
  double precision = 10e-10; //10, needs to be a function of abs value
  
  // Transformation of the utility index matrix
  //mat Z = exp(Omega/lambda);
  
  
  mat G = Omega/lambda;
  
  rowvec maxl = max(Omega/lambda,0);
  mat Z = zeros<mat>(nalt,nstates);  
  for(int j=0;j<nstates;j++){
    Z.col(j)=exp(Omega.col(j)/lambda - maxl(j));
  }
  
  double delta = 1;
  int counter = 0;
  
  
  //Setting the starting value 
  double optimal_action_weight = 0.75;
  vec prior_choice_prob = Omega * mu; 
  oldProbs.fill((1-optimal_action_weight)/(nalt-1));
  oldProbs(prior_choice_prob.index_max()) = optimal_action_weight;
    
  // oldProbs.fill(1/nalt);
  // oldCondProbs.fill(1/nalt);
  
  while(delta>precision && counter<max_iter){
  
    // step 1
    //for(int j=0;j<nalt;j++){
    //  for(int k=0;k<nstates;k++){
    //    vec temp = trans(oldProbs) * Z.col(k);
    //    double D1 = temp(0);
    //    cCondProbs(j,k) = oldProbs(j) * Z(j,k) / D1; 
    //  }
    //} 
    for(int j=0;j<nalt;j++){ //numerically stable in logs
      for(int k=0;k<nstates;k++){
        vec temp = maxl(k) + log(trans(oldProbs) * Z.col(k)); //LSE
        double D1 = temp(0);
        cCondProbs(j,k) = log(oldProbs(j)) + G(j,k) - D1; //in logs
        cCondProbs(j,k) = exp(cCondProbs(j,k));
      }
    }
    
    
    // step 2
    cProbs = cCondProbs * mu;

    delta = - log(sum(sqrt(cProbs % oldProbs))); //sum(pow(cProbs - oldProbs,2)); 
    counter++;
    oldProbs = cProbs;
  }
  
  // Calculating the posterior beliefs
  for(int j=0;j<nalt;j++){
    if(cProbs(j) < 1e-10) cProbs(j) = 1e-10;
    post.row(j) = cCondProbs.row(j)/cProbs(j); //cProbs can be zero
    post.row(j) = post.row(j) % mu.t();
  }
  
  return List::create(Named("probs") = cProbs, Named("condprobs") = cCondProbs, Named("post") = post,
                      Named("Omega") = Omega);
}