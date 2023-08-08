#include "RcppArmadillo.h"
using namespace arma; // use the Armadillo library for matrix computations
using namespace Rcpp;
using namespace R; 

// [[Rcpp::depends(RcppArmadillo)]]


//timing functions
time_t itime;
char buf[100];

void startMcmcTimer() {
  itime = time(NULL);
  Rcout << " MCMC Iteration (est time to end - min) \n";
}

void infoMcmcTimer(int rep, int R) {
  time_t ctime = time(NULL);    
  char buf[32];
  
  double timetoend = difftime(ctime, itime) / 60.0 * (R - rep - 1) / (rep+1);
  sprintf(buf, " %d (%.1f)\n", rep+1, timetoend);
  Rcout <<  buf;
}

void endMcmcTimer() {
  time_t ctime = time(NULL);
  char buf[32];
  
  sprintf(buf, " Total Time Elapsed: %.2f \n", difftime(ctime, itime) / 60.0);     
  Rcout << buf;
  
  itime = 0;
}

//necessary functions for loop

//RI likelihood
// [[Rcpp::export]]
List RI_lik_single_looped(List Data, vec beta, double lambda, vec mu, mat states, int nobs, 
               uvec simp_idx, uvec comp_idx, mat optimal_action_weight2){  
  
  // Computes the joint likelihood of the RI choices of an individual 
  
  //// Arash Laghaie 2020
  //// Matteo Fina 2021
  int nstates = states.n_cols;
  double nalt = states.n_rows;   
  // int nvar = beta.n_elem; 
  // dim(optimal_action_weight)=c(nvar,nobs)
  vec lik(nobs);
  
  vec oldProbs = zeros<vec>(nalt);
  vec cProbs = oldProbs;
  mat oldCondProbs = zeros<mat>(nalt,nstates); 
  mat cCondProbs = oldCondProbs;
  mat post = zeros<mat>(nalt,nstates);
  mat outprobs(nalt,nobs);
  
  for(int i=0;i<nobs;i++){
    
    mat X = as<mat>(as<List>(Data[i])["X"]);
    int y = as<int>(as<List>(Data[i])["y"]);
    int v = as<int>(as<List>(Data[i])["v"]);
    
    // GENERAL VERSION
    // complex element(s) of beta and complex column(s) of X (indicated by complex_idx) are separated from the rest 
    // and then Omega is calculated accordingly
    mat X_simp(X.n_rows,simp_idx.n_elem);
    X_simp = X.cols(simp_idx);
    
    vec beta_simp(simp_idx.n_elem);
    beta_simp = beta.elem(simp_idx);
    
    // Calculate Omega (matrix of utility values at each state) | beta
    mat Xbeta = X.cols(simp_idx)*beta(simp_idx); 
    mat Omega = zeros<mat>(nalt,nstates);
    for(int j=0;j<nstates;j++) Omega.col(j) = Xbeta;
    //Omega.each_col()=Xbeta;
    //adding the utility of the complex attribute 
    Omega = Omega + conv_to<double>::from(beta(comp_idx)) * states;  //here beta(comp_idx) needs to converted to double for multiplication. change this when there are more than 1 simple attributes 
    
    // Calculate P(v) | Omega, mu, lambda  using BA algorithm
    // max_iter: maximal number of iterations 
    // precision: threshold for the sum of squared differences
    int max_iter = 100000000;
    double precision = 10e-10;
    
    // Transformation of the utility index matrix
    //mat Z = exp(Omega/lambda);
    
    //numerically stable
    mat G = Omega/lambda;
    
    rowvec maxl = max(Omega/lambda,0); //vectorwise maximum
    mat Z = zeros<mat>(nalt,nstates);  
    for(int j=0;j<nstates;j++){
      Z.col(j)=exp(Omega.col(j)/lambda - maxl(j));
    }
    
    double delta = 1;
    int counter = 0;
    
    //Starting value for BLAH 
    // double optimal_action_weight = 0.75; //input as vec.length()=nobs
    // vec prior_choice_prob = Omega * mu;
    // oldProbs.fill((1-optimal_action_weight)/(nalt-1));
    // oldProbs(prior_choice_prob.index_max()) = optimal_action_weight;
    
    oldProbs = optimal_action_weight2.col(i);
    
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
      
      delta = - log(sum(sqrt(cProbs % oldProbs))); // sum(pow(cProbs - oldProbs,2)); 
      counter++;
      oldProbs = cProbs;
    }
    outprobs.col(i)=cProbs;
    lik(i) = log(cCondProbs(y-1,v-1));
    if(lik(i) < -1e10) lik(i) = -1e10;
    //if(ISNAN(lik(i)) == TRUE) lik(i) = -1e10;
  }
  // lik.print();
  //return(sum(lik));
  
  return List::create(
    Named("lik") = sum(lik),
    Named("probs") = outprobs);
}

//function to draw from mv normal (only one k-dimensional MV draw)
//mat mvrnormArma(int n, vec mu, mat sigma) {
//  int ncols = sigma.n_cols;
//  mat Y = randn(n, ncols);
//  return repmat(mu, 1, n).t() + Y * chol(sigma);
//}

vec mvrnormArma(vec mu, mat sigma) {
  int ncols = sigma.n_cols;
  mat Y = randn(1, ncols);
  return vectorise(repmat(mu, 1,1).t() + Y * chol(sigma));
}

//function to compute normal density
double lndMvn(vec const& x, vec const& mu, mat const& rooti){
  
  //Wayne Taylor 9/7/2014
  
  // function to evaluate log of MV Normal density with  mean mu, var Sigma
  // Sigma=t(root)%*%root   (root is upper tri cholesky root)
  // Sigma^-1=rooti%*%t(rooti)   
  // rooti is in the inverse of upper triangular chol root of sigma
  //          note: this is the UL decomp of sigmai not LU!
  //                Sigma=root'root   root=inv(rooti)
  
  vec z = vectorise(trans(rooti)*(x-mu));
  
  return((-(x.size()/2.0)*log(2*M_PI) -.5*(trans(z)*z) + sum(log(diagvec(rooti))))[0]);
}



//Main loop

// [[Rcpp::export]]
List RI_MH_singlecomplex_cpp(List Data,int R,
                       double lambda,uvec simp_idx,uvec comp_idx,mat states, //RI primitives
                       mat stepsizes,vec beta_start, //stepsize and individual starting value
                       vec priomean, mat priovar) //priors
{
  
  //MH algorithm to draw from posterior of coefficients in a RI DCM
  
  //output:
  // betadraw: cube with dim=c(individuals, coefficient, MCMC draws)
  // lhdraw:   likelihood draws for each RI individual with dim=c(individuals, MCMC draws)
  // hiermean:
  // hiervar:
  
  
  //int nind =1;//Data.size();
  int nobs =Data.size(); //Data[0][0].length();
  int nalt = states.n_rows;
  int nvar = as<mat>(as<List>(Data[0])["X"]).n_cols;
  //int ncomplevels = 5;//rho.n_elem;
  int nstates = states.n_cols;//rho.size();//ncomplevels^(nalt-1);  // number of all possible states

  //simp_idx = simp_idx - 1;   //adopt cpp indexing
  //comp_idx = comp_idx - 1;   //adopt cpp indexing
  
  //creating storage for results
  //cube betadraw(nind,nvar,R),hiervar(nvar,nvar,R);
  //mat likdraw(nind,R),hiermean(nvar,R); 
  //vec oldlik(nind);
  
  mat betadraw(nvar,R);
  vec likdraw(R);
  
  // mu (prior over states) assuming a uniform prior
  vec mu(nstates);
  mu.fill(1.0/nstates);
  
  // initiliazing starting values for betadraws
  //mat betastar(nind, nvar);
  //for(int i=0;i<nind;i++){
  //  betastar.row(i) = beta_start;
  //}
  vec betastarc=beta_start;
  vec betastar=beta_start;

  //starting values for hierarchical prior     
  //vec  currmean=hier_beta_start;
  //mat  currvar=hier_var_start; 
  
  mat  rooti=inv(chol(priovar));
  //mat rooti=priovar;
  //variables needed for MCMC
  double alpha, unif, clik, clpostbeta, ldiff, oldlik, oldlpostbeta;
  //vec oldlpostbeta(nind), 
  vec alphaminv(2);
  List BAout;
  
  //initial evaluation of likelihood
  //  for(int i=0;i<nind;i++){
  //  oldlik(i) = RIDC2_lik(Data(i), trans(betastar.row(i)), lambda, mu, states, nobs,simp_idx, comp_idx) ; 
  //  oldlpostbeta(i) = oldlik(i) +  lndMvn(trans(betastar.row(i)), currmean,rooti); 
  //    }
  
  //cube BA_start(nalt,nobs,nind);
  mat BA_start(nalt,nobs);
  BA_start.fill(1.0/nalt); //1/nalt
  
  //double weight=1.0/nalt;
  
  //for(int i=0;i<nind;i++){
    
  //  BAout = RIDC2_lik(Data(i), trans(betastar.row(i)), lambda, mu, states, nobs,simp_idx, comp_idx,
  //                    BA_start.slice(i)) ;
  //  BA_start.slice(i)=as<mat>(BAout["probs"]);//*0.99;
  //  oldlik(i) = BAout["lik"];
  //  oldlpostbeta(i) = oldlik(i) +  lndMvn(trans(betastar.row(i)), currmean,rooti); 
  //}
  
    BAout = RI_lik_single_looped(Data, beta_start, lambda, mu, states, nobs,simp_idx, comp_idx,
                      BA_start);
    //BA_start=as<mat>(BAout["probs"])*0.5+(1/nalt)*0.5;//*0.99;
    //BA_start=as<mat>(BAout["probs"]);//*0.99;
    oldlik = BAout["lik"];
    oldlpostbeta= oldlik +  lndMvn(beta_start, priomean,rooti); 
  
  
  startMcmcTimer();
  
  //main iteration loop
  for(int reps=0;reps<R;reps++){
    
      //propose betastar 
      betastarc = mvrnormArma(betastar,stepsizes);
      
      //calculate likelihood + posterior  
      BAout = RI_lik_single_looped(Data, betastarc, lambda, mu, states, nobs,simp_idx, comp_idx,
                        BA_start) ;
      
      clik = BAout["lik"];
      clpostbeta =  clik +  lndMvn(betastarc, priomean,rooti); //rooti of currvar
      
      //acceptance ratio
      ldiff = clpostbeta - oldlpostbeta;
      alphaminv = {1.0, exp(ldiff)};
      alpha = min(alphaminv);
      
      if(alpha < 1.0) {
        unif = randu();
      }else{
        unif = 0.0;
      }
      if (unif <= alpha){ 
        betastar= betastarc;
        oldlik = clik;
        oldlpostbeta = clpostbeta;
        //update Blahut Arimoto starting value for next MCMC draw
        //BA_start=as<mat>(BAout["probs"])*0.5+(1/nalt)*0.5;//*0.99;
      }
      
     //i-break
    
    if ((reps+1)%100==0) infoMcmcTimer(reps, R);
    
    
    betadraw.col(reps) = betastar; 
    likdraw(reps)=oldlik;
  } //R break
  endMcmcTimer();
  
  return List::create(
    Named("betadraw") = betadraw,
    Named("lhdraw") = likdraw);
  
}
