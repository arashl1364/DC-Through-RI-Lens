rm(list = ls())
setwd("~/GitHub/RI-Mama/Cost_Identification")


library(bayesm)
library(Rcpp)
library(RcppArmadillo)

source('RI_DCM2_SIM_new.R')

Rcpp::sourceCpp('BA_stable.cpp')
Rcpp::sourceCpp('RI_MH_singlecomplex_cpp.cpp')
Rcpp::sourceCpp('RI_MH_singlecomplex_cpp_cat.cpp')

###############################################
############### Simulation  ###################
############################################### 
# Simulation from RI model for discrete choice : 
#           - # of alternatives: 1 inside goods, 1 outside good
#           - # of attributes: 1 price (simple), 1 discount (complex)
#           - realization of the complex attribute is independent of other attribute levels 
#           - rho are vectors of attribute levels
#
# Generating data
#number of choices
N = 2000
#cost of information processing
lambda001 = 0.001 #small lambda, perfect learning
lambda05  = 0.5   #intermediate lambda, partial learning
lambda5   = 5     #large lambda, no learning
#preferences
betas= c(-1,1) 

#attribute levels
rho_simp = c(2.25,2.75,5) #true simple
rho_comp =  c( 1,1.5,3.5) #true complex
#rho_comp =  c( 1,2.5,3.5) #true complex

#possible realizations of complex attributes for both alternatives
states_correct=rbind(rho_comp,rep(0,length(rho_comp)))
states_wrong=rbind(rho_simp,rep(0,length(rho_simp)))

#################################################################################
#################### CORRECT SPECIFICATION, Price=simple ########################
#################################################################################

#simulate choices 
set.seed(66)
sim.list001 = RI_DCM2_Sim_new(betas, lambda001, rho_comp, 
                              simplelvls = rho_simp, N, 2, brand=FALSE)
set.seed(66)
sim.list05  = RI_DCM2_Sim_new(betas, lambda05, rho_comp, 
                              simplelvls = rho_simp, N, 2, brand=FALSE)
set.seed(66)
sim.list5   = RI_DCM2_Sim_new(betas, lambda5, rho_comp, 
                              simplelvls = rho_simp, N, 2, brand=FALSE)
# create Data lists
Data001<-Data05<-Data5<- vector(mode = "list", N)
for(i in 1:N){
  Data001[[i]]$X = sim.list001[[i]]$X
  Data001[[i]]$y = sim.list001[[i]]$choice
  Data001[[i]]$v = sim.list001[[i]]$curr_state
  
  Data05[[i]]$X = sim.list05[[i]]$X
  Data05[[i]]$y = sim.list05[[i]]$choice
  Data05[[i]]$v = sim.list05[[i]]$curr_state
  
  Data5[[i]]$X = sim.list5[[i]]$X
  Data5[[i]]$y = sim.list5[[i]]$choice
  Data5[[i]]$v = sim.list5[[i]]$curr_state
} 

###############################################
############### Estimation  ###################
###############################################

####################################################
############## SMALL LAMBDA = 0.001 ################
####################################################
set.seed(66) #setting seed for replication
out001=RI_MH_singlecomplex_cpp(Data=Data001,R=120000,lambda=lambda001,
                            simp_idx=c(0),comp_idx=c(1), #cpp indexing
                            states=states_correct,
                            stepsizes=diag(0.1,2),beta_start=c(-2,2),
                            priomean=c(0,0),priovar=diag(100,2))

#results
windows()
par(mfrow=c(3,1))
plot(out001$betadraw[1,])
plot(out001$betadraw[2,])
plot(out001$lhdraw)
mean(out001$lhdraw[3000:120000])


####################################################
########### INTERMEDIATE LAMBDA = 0.5 ##############
####################################################

set.seed(66) #setting seed for replication
out05=RI_MH_singlecomplex_cpp(Data=Data05,R=120000,lambda=lambda05,
                               simp_idx=c(0),comp_idx=c(1), #cpp indexing
                               states=states_correct,
                               stepsizes=diag(0.001,2),beta_start=c(-2,2),
                               priomean=c(0,0),priovar=diag(100,2))

#results
windows()
par(mfrow=c(3,1))
plot(out05$betadraw[1,])
plot(out05$betadraw[2,])
plot(out05$lhdraw)
mean(out05$lhdraw[3000:120000])

####################################################
############### LARGE LAMBDA = 5 ###################
####################################################
set.seed(66) #setting seed for replication
out5=RI_MH_singlecomplex_cpp(Data=Data5,R=120000,lambda=lambda5,
                              simp_idx=c(0),comp_idx=c(1), #cpp indexing
                              states=states_correct,
                              stepsizes=diag(0.1,2),beta_start=c(-2,2),
                              priomean=c(0,0),priovar=diag(100,2))

#results
windows()
par(mfrow=c(3,1))
plot(out5$betadraw[1,])
plot(out5$betadraw[2,])
plot(out5$lhdraw)
mean(out5$lhdraw[3000:120000])



#################################################################################
#################### WRONG SPECIFICATION, Price=complex #########################
#################################################################################

# create Data lists
Data_wrg_spec001 <- Data001
Data_wrg_spec05  <- Data05
Data_wrg_spec5   <- Data5

row.names(states_wrong)=c("alt1","outside") #needed to find the state in wrong spec
nstates=ncol(states_wrong) #number of distinct states
comp_idx_wrg = 1   # index of the wrongly specified complex attribute

# finding the state number in the wrong specification
for(i in 1:N){
  for(j in 1:nstates){
    if(isTRUE(all.equal(Data_wrg_spec001[[i]]$X[,comp_idx_wrg],states_wrong[,j]))){Data_wrg_spec001[[i]]$v=j}  
    if(isTRUE(all.equal(Data_wrg_spec05[[i]]$X[,comp_idx_wrg],states_wrong[,j]))){Data_wrg_spec05[[i]]$v=j}  
    if(isTRUE(all.equal(Data_wrg_spec5[[i]]$X[,comp_idx_wrg],states_wrong[,j]))){Data_wrg_spec5[[i]]$v=j}  
    }
}
###############################################
############### Estimation  ###################
###############################################

####################################################
############## SMALL LAMBDA = 0.001 ################
####################################################
set.seed(66) #setting seed for replication
out_wrg_spec001 = RI_MH_singlecomplex_cpp(Data=Data_wrg_spec001,R=120000,
                                          lambda=lambda001,simp_idx=c(1),comp_idx=c(0),
                                          states=states_wrong,
                                          stepsizes=diag(0.1,2),beta_start=c(-2,2),
                                          priomean=c(0,0),priovar=diag(100,2))
#results
windows()
par(mfrow=c(3,1))
plot(out_wrg_spec001$betadraw[1,])
plot(out_wrg_spec001$betadraw[2,])
plot(out_wrg_spec001$lhdraw)

mean(out_wrg_spec001$betadraw[1,3000:120000])
mean(out_wrg_spec001$betadraw[2,3000:120000])
mean(out_wrg_spec001$lhdraw[3000:120000])


####################################################
########## INTERMEDIATE LAMBDA = 0.5 ###############
####################################################
set.seed(66) #setting seed for replication
out_wrg_spec05=RI_MH_singlecomplex_cpp(Data=Data_wrg_spec05,R=120000,
                                      lambda=lambda05,simp_idx=c(1),comp_idx=c(0),
                                      states=states_wrong,
                                      stepsizes=diag(0.0001,2),beta_start=c(-2,2),
                                      priomean=c(0,0),priovar=diag(100,2))

#results
windows()
par(mfrow=c(3,1))
plot(out_wrg_spec05$betadraw[1,])
plot(out_wrg_spec05$betadraw[2,])
plot(out_wrg_spec05$lhdraw)
mean(out_wrg_spec05$lhdraw[20000:120000])

####################################################
############### LARGE LAMBDA = 5 ###################
####################################################
set.seed(66) #setting seed for replication
out_wrg_spec5=RI_MH_singlecomplex_cpp(Data=Data_wrg_spec5,R=120000,
                                       lambda=lambda5,simp_idx=c(1),comp_idx=c(0),
                                       states=states_wrong,
                                       stepsizes=diag(0.1,2),beta_start=c(-2,2),
                                       priomean=c(0,0),priovar=diag(100,2))
#results
windows()
par(mfrow=c(3,1))
plot(out_wrg_spec5$betadraw[1,])
plot(out_wrg_spec5$betadraw[2,])
plot(out_wrg_spec5$lhdraw)
mean(out_wrg_spec5$lhdraw[3000:120000])


################################################################################
########################### RESULTS SECTION 3.3 ################################
################################################################################

#table 4
# summary.bayesm.mat(out001$lhdraw)
# summary.bayesm.mat(out_wrg_spec001$lhdraw)
# 
# summary.bayesm.mat(out05$lhdraw)
# summary.bayesm.mat(out_wrg_spec05$lhdraw)
# 
# summary.bayesm.mat(out5$lhdraw)
# summary.bayesm.mat(out_wrg_spec5$lhdraw)

quantile(out001$lhdraw[20001:120000]) #burnin 20k
quantile(out_wrg_spec001$lhdraw[20001:120000]) #burnin 20k

quantile(out05$lhdraw[20001:120000]) #burnin 20k
quantile(out_wrg_spec05$lhdraw[20001:120000]) #burnin 20k

quantile(out5$lhdraw[20001:120000]) #burnin 20k
quantile(out_wrg_spec5$lhdraw[20001:120000]) #burnin 20k


#Figure 1
library(ggplot2)

#ggplot inputs
betadraw001=out001$betadraw[2,20001:120000] #burnin 20k
betadraw05=out05$betadraw[2,20001:120000]
betadraw5=out5$betadraw[2,20001:120000]
type=c(rep("0.01",length(betadraw001)),rep("0.5",length(betadraw05)),rep("5",length(betadraw5)))
frame<-data.frame(list("betadraws"=as.vector(c(betadraw001,betadraw05,betadraw5)), "Type"=type))

figure1<-ggplot(frame, aes(x=betadraws, linetype=type, fill=type)) +
  geom_histogram(position="identity",alpha=0.4, size = 1, color="black", bins=80)+
  #geom_vline(mean(betadraw001), size=2, linetype="dashed")+
  #geom_vline(mean(betadraw05), size=2, linetype="dashed")+
  #geom_vline(mean(betadraw5), size=2, linetype="dashed")+
  geom_vline(xintercept = 1, size=0.8)+
  scale_linetype_manual(values=c("solid","dashed", "dotted"), labels = c("\u03BB = 0.01","\u03BB = 0.5","\u03BB = 5"))+
  scale_color_manual(values=c( "#F8766D","#00BFC4","#00BA38"), labels = c("\u03BB = 0.01","\u03BB = 0.5","\u03BB = 5"))+
  scale_fill_manual(values=c("#F8766D","#00BFC4","#00BA38"), labels = c("\u03BB = 0.01","\u03BB = 0.5","\u03BB = 5"))+
  
  labs(title="Posterior Draws Complex Coefficient",x="Complex Coefficient",y="Frequency")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.title = element_blank(),
        legend.position = "top",
        legend.key.size = unit(2, "lines"),
        legend.text = element_text(size=14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=14, face = "bold"))+
  coord_cartesian( ylim = c(0, 107000))


windows()
figure1


################################################################################
############################ APPENDIX RESULTS ##################################
################################################################################

#We estimate a categorical specification of the wrong model with intermediate lambda=0.5
#Baseline is rho_comp[1]

#appendix C, table 14
Data_cat=Data_wrg_spec05

for (i in 1:length(Data_wrg_spec05)){
  Data_cat[[i]]$X=matrix(rep(0,12),nrow=2,ncol=6)
  
  if(Data_wrg_spec05[[i]]$X[1,1]==rho_simp[1]){
    Data_cat[[i]]$X[1,1]=1
    
  }
  if(Data_wrg_spec05[[i]]$X[1,1]==rho_simp[2]){
    Data_cat[[i]]$X[1,2]=1
    
  }
  if(Data_wrg_spec05[[i]]$X[1,1]==rho_simp[3]){
    Data_cat[[i]]$X[1,3]=1
    
  }
  if(Data_wrg_spec05[[i]]$X[1,2]==rho_comp[2]){
    Data_cat[[i]]$X[1,5]=1
    
  }
  if(Data_wrg_spec05[[i]]$X[1,2]==rho_comp[3]){
    Data_cat[[i]]$X[1,6]=1
    
  }
  
}

set.seed(66)
# states2=rbind(c(-2,-3,-4),rep(0,3))
# beta_start2=c(-2,-3,-4,1,2,3)

states_cat=rbind(c(-0.5,-1,-3),rep(0,3)) #in categorical MH, this will be updated with every draw
beta_start_cat=c(-0.5,-1,-3,1,0.1,0.5)
set.seed(66) #for replication
out_wrg_cat=RI_MH_singlecomplex_cpp_cat(Data=Data_cat,R=15000,
                                        lambda=lambda05,simp_idx=c(3,4,5),comp_idx=c(0,1,2),
                                        states=states_cat,
                                        stepsizes=diag(0.001,6),beta_start=beta_start_cat,
                                        priomean=rep(0,6),priovar=diag(100,6))
windows()
par(mfrow=c(3,2))
plot(out_wrg_cat$betadraw[1,])
plot(out_wrg_cat$betadraw[2,])
plot(out_wrg_cat$betadraw[3,])
#plot(out_wrg_cat$betadraw[4,]) #baseline
plot(out_wrg_cat$betadraw[5,])
plot(out_wrg_cat$betadraw[6,])
plot(out_wrg_cat$lhdraw)


#likelihood mean comparison to true linear coding
mean(out05$lhdraw[20000:120000])
mean(out_wrg_cat$lhdraw[2000:15000])
#mean(out_wrg_cat2$lhdraw[2000:15000])

#Table 15, Appendix C
quantile(out05$lhdraw[20001:120000]) #burnin 20k
quantile(out_wrg_spec05$lhdraw[20001:120000]) #burnin 20k
quantile(out_wrg_cat$lhdraw[2001:15000]) #burnin 2k, only 15k draws for categorical

#Figure 15, Appendix C.1

windows()
par(mfrow=c(3,1))
plot(out001$lhdraw[110001:120000],type="l",ylim=c(-0.2,0.2),lwd=2, xlab="MCMC Draw", ylab="LLH",
     main=expression(paste(lambda," = 0.001")), cex.main=1.5)
lines(out_wrg_spec001$lhdraw[110001:120000],col="red", lwd=2)

plot(out05$lhdraw[110001:120000],type="l",ylim=c(-1100,-300), lwd=2, xlab="MCMC Draw", ylab="LLH", 
     main=expression(paste(lambda," = 0.5")), cex.main=1.5)
lines(out_wrg_spec05$lhdraw[110001:120000],col="red", lwd=2)

plot(out5$lhdraw[110001:120000],type="l",ylim=c(-0.5,0.5),lwd=2, xlab="MCMC Draw", ylab="LLH", 
     main=expression(paste(lambda," = 5")), cex.main=1.5)
lines(out_wrg_spec5$lhdraw[110001:120000],col="red", lwd=2)

