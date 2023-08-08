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

//wishart draws for rmultireg

List rwishart(double nu, mat const& V){
  
  int m = V.n_rows;
  mat T = zeros(m,m);
  
  for(int i = 0; i < m; i++) {
    T(i,i) = sqrt(rchisq(1,nu-i)[0]); //rchisq returns a vectorized object, so using [0] allows for the conversion to double
  }
  
  for(int j = 0; j < m; j++) {  
    for(int i = j+1; i < m; i++) {    
      T(i,j) = rnorm(1)[0]; //rnorm returns a NumericVector, so using [0] allows for conversion to double
    }}
  
  mat C = trans(T)*chol(V);
  mat CI = solve(trimatu(C),eye(m,m)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  
  
  // C is the upper triangular root of Wishart therefore, W=C'C
  // this is the LU decomposition Inv(W) = CICI' Note: this is
  // the UL decomp not LU!
  
  // W is Wishart draw, IW is W^-1
  
  return List::create(
    Named("W") = trans(C) * C,
    Named("IW") = CI * trans(CI),
    Named("C") = C,
    Named("CI") = CI);
}

//rmultireg multiple regression
List rmultireg_mod(mat const& Y, mat const& X, mat const& Bbar, mat const& A, double nu, mat const& V) {
  //Note that we rename the function to avoid problems in shared libraries/R environment
  
  // Keunwoo Kim 09/09/2014
  
  // Purpose: draw from posterior for Multivariate Regression Model with natural conjugate prior
  
  // Arguments:
  //  Y is n x m matrix
  //  X is n x k
  //  Bbar is the prior mean of regression coefficients  (k x m)
  //  A is prior precision matrix
  //  nu, V are parameters for prior on Sigma
  
  // Output: list of B, Sigma draws of matrix of coefficients and Sigma matrix
  
  // Model: 
  //  Y=XB+U  cov(u_i) = Sigma
  //  B is k x m matrix of coefficients
  
  // Prior:  
  //  beta|Sigma  ~ N(betabar,Sigma (x) A^-1)
  //  betabar=vec(Bbar)
  //  beta = vec(B) 
  //  Sigma ~ IW(nu,V) or Sigma^-1 ~ W(nu, V^-1)
  
  int n = Y.n_rows;
  int m = Y.n_cols;
  int k = X.n_cols;
  
  //first draw Sigma
  mat RA = chol(A);
  mat W = join_cols(X, RA); //analogous to rbind() in R
  mat Z = join_cols(Y, RA*Bbar);
  // note:  Y,X,A,Bbar must be matrices!
  mat IR = solve(trimatu(chol(trans(W)*W)), eye(k,k)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  // W'W = R'R  &  (W'W)^-1 = IRIR'  -- this is the UL decomp!
  mat Btilde = (IR*trans(IR)) * (trans(W)*Z);
  // IRIR'(W'Z) = (X'X+A)^-1(X'Y + ABbar)
  mat E = Z-W*Btilde;
  mat S = trans(E)*E;
  // E'E
  
  // compute the inverse of V+S
  mat ucholinv = solve(trimatu(chol(V+S)), eye(m,m));
  mat VSinv = ucholinv*trans(ucholinv);
  
  List rwout = rwishart(nu+n, VSinv);
  
  // now draw B given Sigma
  //   note beta ~ N(vec(Btilde),Sigma (x) Covxxa)
  //       Cov=(X'X + A)^-1  = IR t(IR)  
  //       Sigma=CICI'    
  //       therefore, cov(beta)= Omega = CICI' (x) IR IR' = (CI (x) IR) (CI (x) IR)'
  //  so to draw beta we do beta= vec(Btilde) +(CI (x) IR)vec(Z_mk)  
  //  		Z_mk is m x k matrix of N(0,1)
  //	since vec(ABC) = (C' (x) A)vec(B), we have 
  //		B = Btilde + IR Z_mk CI'
  
  mat CI = rwout["CI"]; //there is no need to use as<mat>(rwout["CI"]) since CI is being initiated as a mat in the same line
  mat draw = mat(rnorm(k*m));
  draw.reshape(k,m);
  mat B = Btilde + IR*draw*trans(CI);
  
  return List::create(
    Named("B") = B, 
    Named("Sigma") = rwout["IW"]);
}

//RI likelihood
// [[Rcpp::export]]
List RIDC2_lik(List Data, vec beta, double lambda, vec mu, mat states, int nobs, 
                 uvec simp_idx, uvec comp_idx, mat optimal_action_weight){  
  
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
    Omega.each_col()=Xbeta;
    //adding the utility of the complex attribute 
    Omega = Omega + conv_to<double>::from(beta(comp_idx)) * states;  //here beta(comp_idx) needs to converted to double for multiplication. change this when there are more than 1 simple attributes 
    
    // Calculate P(v) | Omega, mu, lambda  using BA algorithm
    // max_iter: maximal number of iterations 
    // precision: threshold for the sum of squared differences
    int max_iter = 10000;
    double precision = 10e-7;
    
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
    //double optimal_action_weight = 0.75; //input as vec.length()=nobs
    //vec prior_choice_prob = Omega * mu;
    //oldProbs.fill((1-optimal_action_weight)/(nalt-1));
    //oldProbs(prior_choice_prob.index_max()) = optimal_action_weight;
    
    oldProbs = optimal_action_weight.col(i);
    
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
List RIDC2_MH_hier_cpp(List Data,int R,
                       vec rho,double lambda,uvec simp_idx,uvec comp_idx,mat states, //RI primitives
                       mat Bbar,mat Areg,int NUreg,mat Vreg, //reg priors
                       mat propvar,rowvec beta_start, vec hier_beta_start, //MCMC inputs
                       mat hier_var_start)
{
 int nind =Data.size();
 //List obs = Data[0];
 int nobs =20; //Data[0][0].length();
 int nalt = 2;//Data[0][0]["X"].nrow();
 int nvar = 2;//Data[0][0]["X"].ncol();
 //int ncomplevels = 5;//rho.n_elem;
 int nstates = 5;//21;//ncomplevels^(nalt-1);  // number of all possible states
 //IntegerVector simp_idx = {1,2}; //c(1:nvar)[-comp_idx] // index of simple attributes
  
  //simp_idx = simp_idx - 1;   //adopt cpp indexing
  //comp_idx = comp_idx - 1;   //adopt cpp indexing
  
//creating storage for results
  cube betadraw(nind,nvar,R),hiervar(nvar,nvar,R);
  mat likdraw(nind,R),hiermean(nvar,R); 
  vec oldlik(nind);

  vec testlik(nind);
    
// mu (prior over states) assuming a uniform prior
  vec mu(nstates);
  mu.fill(1.0/nstates);

// initiliazing starting values for betadraws
  mat betastar(nind, nvar);
  for(int i=0;i<nind;i++){
    betastar.row(i) = beta_start;
    }
  vec betastarc=trans(beta_start); //makes betastar rowvec?
  
  //starting values for hierarchical prior     
  vec  currmean=hier_beta_start;
  mat  currvar=hier_var_start; 
  mat  rooti=inv(chol(currvar));

  //variables needed for MCMC
  double alpha, unif, clik, clpostbeta, ldiff;
  vec oldlpostbeta(nind), alphaminv(2);
  List BAout;

//initial evaluation of likelihood
//  for(int i=0;i<nind;i++){
//  oldlik(i) = RIDC2_lik(Data(i), trans(betastar.row(i)), lambda, mu, states, nobs,simp_idx, comp_idx) ; 
//  oldlpostbeta(i) = oldlik(i) +  lndMvn(trans(betastar.row(i)), currmean,rooti); 
//    }
cube BA_start(nalt,nobs,nind);
BA_start.fill(1.0/nalt); //1/nalt
vec betalik(nvar+1);
//double weight=1.0/nalt;

for(int i=0;i<nind;i++){
  //BAout = RIDC2_lik(Data(i), trans(betastar.row(i)), lambda, mu, states, nobs,simp_idx, comp_idx,
  //                  BA_start.slice(i)) ;
  //for (int l=0;l<(nvar+1);l++){
  //  betalik(l)=betastar(i,l);
  //  if (l==nvar) {betalik(l)=-betastar(i,(l-1));}
  //}
  betalik(0)=betastar(i,0);
  betalik(1)=betastar(i,1);
  betalik(2)=-betastar(i,1);
  
  BAout = RIDC2_lik(Data(i), betalik, lambda, mu, states, nobs,simp_idx, comp_idx,
                    BA_start.slice(i)) ;
  BA_start.slice(i)=as<mat>(BAout["probs"])*0.8+0.1;
  oldlik(i) = BAout["lik"];
  oldlpostbeta(i) = oldlik(i) +  lndMvn(trans(betastar.row(i)), currmean,rooti); 
}

  startMcmcTimer();
  
  //main iteration loop
  for(int reps=0;reps<R;reps++){
    for(int i=0;i<nind;i++){

  // JUST FOR JOINT MODEL!
  //propose betastar per individual
  betastarc = mvrnormArma(trans(betastar.row(i)),propvar);
  //betastarc = vectorise(betastarc); //mvrnormArma returns a matrix, this makes it a colvec
  //for (int l=0;l<(nvar+1);l++){
  //  betalik(l)=betastar(i,l);
  //  if (l==nvar) {betalik(l)=-betastar(i,(l-1));}
  //}
  betalik(0)=betastarc(0);
  betalik(1)=betastarc(1);
  betalik(2)=-betastarc(1);
  //calculate likelihood + posterior  
  BAout = RIDC2_lik(Data(i), betalik, lambda, mu, states, nobs,simp_idx, comp_idx,
                    BA_start.slice(i)) ;
  //BAout = RIDC2_lik(Data(i), betastarc, lambda, mu, states, nobs,simp_idx, comp_idx,
  //                  BA_start.slice(i)) ;
  BA_start.slice(i)=as<mat>(BAout["probs"])*0.8+0.1;//randu();//0.25;//0.5*0.5;//randu();//weight;//(1.0/nalt);
  clik = BAout["lik"];

  rooti=inv(chol(currvar));
  clpostbeta =  clik +  lndMvn(betastarc, currmean,rooti); //rooti of currvar!

  //acceptance ratio
  ldiff = clpostbeta - oldlpostbeta(i);
  //alphaminv << 1.0 << exp(ldiff); //intializes variables in the alphaminv vec: c(1,exp(ldiff))
  alphaminv = {1.0, exp(ldiff)};
  alpha = min(alphaminv);
  
  if(alpha < 1.0) {
    unif = randu();//runif(1)[0]; //runif returns a NumericVector, so using [0] allows for conversion to double
  }else{
    unif = 0.0;
  }
  if (unif <= alpha){ 
    betastar.row(i) = trans(betastarc);
    oldlik(i) = clik;
    oldlpostbeta(i) = clpostbeta;
  }

  } //i-break
    
    if ((reps+1)%100==0) infoMcmcTimer(reps, R);
    
  //Hierarchical level
  
  //creating X matrix for rmultireg
    mat Xreg(nind,1,fill::ones);
    //please check input arguments for rmultireg, see lines 44ff
    List outreg = rmultireg_mod(betastar,Xreg,Bbar,Areg,NUreg,Vreg);
    
      hiermean.col(reps) = as<vec>(as<List>(outreg)["B"]);
      currmean=as<vec>(as<List>(outreg)["B"]); 
      hiervar.slice(reps) = as<mat>(as<List>(outreg)["Sigma"]);
      currvar=hiervar.slice(reps);
      
      betadraw.slice(reps) = betastar; //betadraw, vardraw are cubes/arrays with three dimensions!!!!
      likdraw.col(reps)=oldlik;
  } //R break
  
  endMcmcTimer();
  
return List::create(
  Named("betadraw") = betadraw,
  Named("lhdraw") = likdraw,
  Named("hiermean") = hiermean,
  Named("hiervar") = hiervar);
  
}
  