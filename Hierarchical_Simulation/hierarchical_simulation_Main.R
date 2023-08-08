#Simulate hierarchical data with fixed choice set:
rm(list=ls())
setwd("C:/Users/mfina/Downloads/hierarchical")
setwd("~/GitHub/RI-Mama/Hierarchical_Simulation")

library("RcppArmadillo")
library(bayesm)
source('RIDC2_SIM.R')
Rcpp::sourceCpp('BA.cpp')
Rcpp::sourceCpp('RI_DCM_Hier_Joint.cpp')


###############################################
############### Simulation  ###################
############################################### 
# Simulation from RI model for discrete choice 
# Version: Feb 2022
set.seed(66)
P = 1000#number of individuals
N = 20  #number of choices per individual
lambda = 0.25

#attribute levels
rho_simple=c(2.5,3,3.5,4,4.5) #simple attribute levels
rho_complex=c(0,0.5,1,1.5,2)  #complex attribute levels

#generating heterogenous preferences
#Note: coefficient for simple attribute = -(coefficient for complex attribute)
set.seed(66)
pricecoefs=rnorm(P,-1,0.2) #drawing price coefficients in population
beta_price=matrix(c(rnorm(P, 2.5, 0.5), #brand coefficients 
                    pricecoefs,
                    -pricecoefs), ncol= 3)

max(pricecoefs) #check if all pricecoefs are negative


#Simulating Data - Version Feb 2022
set.seed(66) #setting seed for replication
SIM_LIST<-NULL
for (i in 1:P){
  SIM_LIST[[i]]=RIDC2_Sim(beta_price[i,], lambda, rho_complex, rho_simple,N, complex="price")
}

#collecting Data in List form
Data = vector(mode="list",P)
for(i in 1:P){
  Data[[i]]=vector(mode="list",N)
}
  
for(j in 1:P){
  for(i in 1:N){
    Data[[j]][[i]]$X = SIM_LIST[[j]][[i]]$X
    Data[[j]][[i]]$y = SIM_LIST[[j]][[i]]$choice
    Data[[j]][[i]]$v = SIM_LIST[[j]][[i]]$curr_state
  } 
}


#creating state space
states=rbind(rho_complex,rep(0,length(rho_complex)))

set.seed(66)
out_joint=RIDC2_MH_hier_cpp(Data=Data, R=50000, rho=rho_complex, lambda=lambda,
                             simp_idx=c(0,1), comp_idx=2, states=states,
                             Bbar=matrix(rep(0,2),ncol=2),
                             Areg=matrix(0.01),NUreg=6,Vreg=diag(6,2), 
                             propvar=diag(0.01,2),
                             beta_start=c(2,-0.5),
                             hier_beta_start=c(2,-0.5),hier_var_start=diag(2))

summary.bayesm.mat(out_joint$hiermean[1,])
summary.bayesm.mat(out_joint$hiermean[2,])




#####################################################################################################################
#################################    MNL ESTIMATION     ################################### 
#####################################################################################################################

###################################
#bring data into bayesm format
###################################
lgtdata<-vector(mode="list", length(Data))
exes<-vector(mode="list", length(Data))
lgtxs<- vector(mode="list", length(Data))
lgtys<-matrix(rep(0,P*N), P,N)
for (i in 1:length(Data)){
  for (j in 1:length(Data[[1]])){
  lgtxs[[i]]$EXES=rbind(exes[[i]],Data[[i]][[j]]$X)
  exes[[i]]=lgtxs[[i]]$EXES
  lgtys[i,j]=Data[[i]][[j]]$y
  }
}

for (i in 1:length(Data)){
  lgtdata[[i]]$X=lgtxs[[i]]$EXES
  lgtdata[[i]]$y=lgtys[i,]
}

#####################################
#estimation
#####################################
data_mnl_sep=list(lgtdata=lgtdata,p=2)
p=2
ncomp=1
mcmc=list(R=50000,nprint=1000,keep=1) 
prior=list(ncomp=1)

out_mnl=rhierMnlRwMixture(Data=data_mnl_sep,Prior=prior, Mcmc=mcmc)

#######################################################
########### Joint Specification #######################
#######################################################
lgtdata<-vector(mode="list", length(Data))
exes<-vector(mode="list", length(Data))
lgtxs<- vector(mode="list", length(Data))
lgtys<-matrix(rep(0,P*N), P,N)
for (i in 1:length(Data)){
  for (j in 1:length(Data[[1]])){
    lgtxs[[i]]$EXES=rbind(exes[[i]],c(1,Data[[i]][[j]]$X[1,2]-Data[[i]][[j]]$X[1,3]), 
                          c(0,0))
    exes[[i]]=lgtxs[[i]]$EXES
    lgtys[i,j]=Data[[i]][[j]]$y
  }
}

for (i in 1:length(Data)){
  lgtdata[[i]]$X=lgtxs[[i]]$EXES
  lgtdata[[i]]$y=lgtys[i,]
}

#########################################
#### estimation of joint specification ##
#########################################
data_mnl_joint=list(lgtdata=lgtdata,p=2)

out_mnl_joint=rhierMnlRwMixture(Data=data_mnl_joint,Prior=prior, Mcmc=mcmc)

##################################################################################
###################################### RESULTS ###################################
##################################################################################


#################################
########### table 2 #############
#################################
#collecting hierarchical results for RU Logit:
hiermeans_mnl=matrix(nrow=mcmc$R,ncol=3)
for (i in 1:nrow(hiermeans_mnl)){
  hiermeans_mnl[i,]=out_mnl$nmix$compdraw[[i]][[1]]$mu
}

hiermeans_mnl_joint=matrix(nrow=mcmc$R,ncol=2)
for (i in 1:nrow(hiermeans_mnl_joint)){
  hiermeans_mnl_joint[i,]=out_mnl_joint$nmix$compdraw[[i]][[1]]$mu
}

#summary hierarchical means
#RI 
summary.bayesm.mat(out_joint$hiermean[1,])
summary.bayesm.mat(out_joint$hiermean[2,])
#separate RU MNL
summary.bayesm.mat(hiermeans_mnl[,1])
summary.bayesm.mat(hiermeans_mnl[,2])
summary.bayesm.mat(hiermeans_mnl[,3])
#joint RU MNL
summary.bayesm.mat(hiermeans_mnl_joint[,1])
summary.bayesm.mat(hiermeans_mnl_joint[,2])


#log marginal density
lh_RI_out=colSums(out_joint$lhdraw)
logMargDenNR(out_mnl$loglike[5000:50000]) #burnin 5k draws
logMargDenNR(out_mnl_joint$loglike[5000:50000])
logMargDenNR(lh_RI_out[5000:50000])

#################################
########### table 3 #############
#################################

#creating spaces
hiervars_RI=matrix(nrow=50000,ncol=2)
hiervars_mnl=matrix(nrow=50000,ncol=3) #rooti's, need to be transformed
hiervars_mnl_joint=matrix(nrow=50000,ncol=2) #rooti's, need to be transformed

hiervars_mnl_array=array(dim=c(3,3,50000))
hiervars_mnl_joint_array=array(dim=c(2,2,50000))

#Sigma = t(root)%*%root where root is upper triangular of chol root
#rooti = inv(root)
#Sigma^-1 = rooti%*%t(rooti)
for (i in 1:50000){
  hiervars_RI[i,]=diag(out_joint$hiervar[,,i])
  hiervars_mnl_array[,,i]=solve(out_mnl$nmix$compdraw[[i]][[1]]$rooti) #computing inverse
  hiervars_mnl_joint_array[,,i]=solve(out_mnl_joint$nmix$compdraw[[i]][[1]]$rooti) #computing inverse
}

for (i in 1:50000){
  hiervars_mnl[i,]=diag(t(hiervars_mnl_array[,,i])%*%hiervars_mnl_array[,,i])
  hiervars_mnl_joint[i,]=diag(t(hiervars_mnl_joint_array[,,i])%*%hiervars_mnl_joint_array[,,i])
  
}

summary.bayesm.mat(hiervars_RI)
summary.bayesm.mat(hiervars_mnl)
summary.bayesm.mat(hiervars_mnl_joint)

#summary.bayesm.nmix(out_mnl$nmix)
#summary.bayesm.nmix(out_mnl_joint$nmix)
 