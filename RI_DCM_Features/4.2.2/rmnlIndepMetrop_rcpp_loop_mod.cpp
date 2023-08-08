#include "RcppArmadillo.h" //Indicates shared library (aka package in R terms).
using namespace arma;      //Needed to avoid interference between libraries.
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]] 



//including necessary CPP(!) functions from bayesm package

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

//function to draw from mv student t distribution
//[[Rcpp::export]]
vec rmvst(double nu, vec const& mu, mat const& root){
  
  //  nu df, mean mu, Sigma=t(root)%*%root
  //  root is upper triangular cholesky root
  
  vec rnormd = rnorm(mu.size());
  vec nvec = trans(root)*rnormd;
  
  return(nvec/sqrt(rchisq(1,nu)[0]/nu) + mu); //rchisq returns a vectorized object, so using [0] allows for the conversion to double
}

//function to compute log likelihood of MNL
//[[Rcpp::export]]
double llmnl(vec const& beta, vec const& y, mat const& X){
  
  
  int n = y.size();
  int j = X.n_rows/n;
  mat Xbeta = X*beta;
  
  vec xby = zeros<vec>(n);
  vec denom = zeros<vec>(n);
  
  for(int i = 0; i<n;i++){      
    for(int p=0;p<j;p++) denom[i]=denom[i]+exp(Xbeta[i*j+p]);
    xby[i] = Xbeta[i*j+y[i]-1];
  }
  
  return(sum(xby - log(denom)));
}

//compute log density of mv normal distribution
//[[Rcpp::export]]
double lndMvn(vec const& x, vec const& mu, mat const& rooti){
  
  // mean mu, var Sigma
  // Sigma=t(root)%*%root   (root is upper tri cholesky root)
  // Sigma^-1=rooti%*%t(rooti)   
  // rooti is in the inverse of upper triangular chol root of sigma
  //          note: this is the UL decomp of sigmai not LU!
  //                Sigma=root'root   root=inv(rooti)
  
  vec z = vectorise(trans(rooti)*(x-mu));
  
  return((-(x.size()/2.0)*log(2*M_PI) -.5*(trans(z)*z) + sum(log(diagvec(rooti))))[0]);
}


//compute log density of mv student t distribution
//[[Rcpp::export]]
double lndMvst(vec const& x, double nu, vec const& mu, mat const& rooti, bool NORMC = false){
  
  
  // nu df, mean mu,
  // and with sigmai=rooti%*%t(rooti)   note: this is the UL decomp of sigmai not LU!
  // rooti is in the inverse of upper triangular chol root of sigma
  // or Sigma=root'root   root=inv(rooti)
  
  int dim = x.size();
  double constant;
  
  if(NORMC){
    constant = (nu/2)*log(nu)+lgamma((nu+dim)/2)-(dim/2.0)*log(M_PI)-lgamma(nu/2);
  } else {
    constant = 0.0;
  }
  
  vec z = vectorise(trans(rooti)*(x-mu));
  
  return((constant-((dim+nu)/2)*log(nu+trans(z)*z)+sum(log(diagvec(rooti))))[0]);
}


//Main loop function to be exported to R environment
//[[Rcpp::export]]
List rmnlIndepMetrop_rcpp_loop_mod(int R, int keep, double nu,
                                vec const& betastar, mat const& root,vec const& y,mat const& X,
                                vec const& betabar,mat const& rootpi,mat const& rooti,
                                double oldlimp,double oldlpost,int nprint) {

// Wayne Taylor 9/7/2014

  int mkeep = 0;
  int naccept = 0;    
  int ncolX = X.n_cols;
  
  mat betadraw(R/keep, ncolX);
  vec loglike(R/keep);
  vec betac = zeros<vec>(ncolX);
  rowvec beta = zeros<rowvec>(ncolX);
  double cloglike, clpost, climp, ldiff, alpha, unif, oldloglike;
  vec alphaminv;
  
  if(nprint>0) startMcmcTimer();
  
  // start main iteration loop
  for(int rep = 0; rep<R; rep++) {
    
    //propose new beta from multivariate student t, where betastar and root are
    //determined from tuning step (max log mnl likelihood)
    betac = rmvst(nu,betastar,root);
    
    //compute likelihood and posterior
    cloglike = llmnl(betac,y,X);
    clpost = cloglike+lndMvn(betac,betabar,rootpi);
    
    //compute acceptance probability
    climp = lndMvst(betac,nu,betastar,rooti,false);
    ldiff = clpost+oldlimp-oldlpost-climp;
    alphaminv << 1 << exp(ldiff); //intializes variables in the alphaminv vec: c(1,exp(ldiff))
    alpha = min(alphaminv);
  
    if(alpha < 1.0) {
        unif = runif(1)[0]; //rnorm returns a NumericVector, so using [0] allows for conversion to double
      }else{
        unif = 0.0;
      }
    if (unif <= alpha){ 
      beta = trans(betac);
      oldloglike = cloglike;
      oldlpost = clpost;
      oldlimp = climp;
      naccept++;
    }
          
    if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);
    
    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      betadraw(mkeep-1,span::all) = beta;
      loglike[mkeep-1] = oldloglike;
    }
  }
  
  if(nprint>0) endMcmcTimer();
      
  return List::create(
    Named("betadraw") = betadraw, 
    Named("loglike") = loglike, 
    Named("naccept") = naccept);
}
