#rm(list = ls())
setwd("C:/Users/Matteo/Desktop/RI/Plots and Tables")

set.seed(69)

source("RICBC_choice_sets_and_incentives.R")
source("RICBC_RI_SC_Choice.R")
library(bayesm)
library(ggplot2)
library(gridExtra)
library(Rcpp)
library(RcppArmadillo)

source("rmnlIndepMetrop_rcpp_mod.R")
sourceCpp("rmnlIndepMetrop_rcpp_loop_mod.cpp")

PriceDraw<-function(pricelevels){
  drawprob<-rep(1/length(pricelevels),length(pricelevels))
  
  out<-sample(pricelevels,1,prob=drawprob)
}



#beta.vec2 <- c(4, -1)
#beta.vec2<-c(2,-1)
beta.vec2<-c(5,-1,1) #4 #3

lambda <- 0.5


#min and max price for uniform draws
prmin=8#7 #4.5
prmax=10#9 #5.5

#price levels for discrete draws
#plvls<-c(4,5.5,6)


#tt<-seq(-1,-6,by=-0.5)
#rhos<-matrix(rep(0,2*length(tt)),ncol=2)
#for (i in 1:length(tt)){
# rhos[i,]<-c(tt[i],0)
#}

#tt<-seq(0.2,2.25,by=0.02)
tt<-seq(0.1,4,by=0.02)

#tt<-c(tt,3)
rhos<-matrix(rep(0,2*length(tt)),ncol=2)
for (i in 1:length(tt)){
  rhos[i,]<-c((4-tt[i]),(4+tt[i])) #2.25
}
#rhos=rbind(c(-2.5,-2), c(-3,-1.5), c(-3.5,-1))

rho.prob <- c(0.5,0.5)


# Simulation --------------------------------------------------------------
N = 1000
L = nrow(rhos)
sim.list1<-array(list(),L)
for (i in 1:L){
  sim.list1[[L]] <-array(list(), N)
}


betas<-matrix(c(rep(0,(length(beta.vec2)+1)*L)),ncol=length(beta.vec2)+1)
standarderrors<-matrix(c(rep(0,(length(beta.vec2)+1)*L)),ncol=length(beta.vec2)+1)
logliks<-matrix(rep(0,L*10000),ncol=10000) #ncol=number of draws

set.seed(777)
#set.seed(61)
#set.seed(33) with seq by 0.05
for(j in 1:L){
  rho<-rhos[j,]
  t1 <- Sys.time()
  for(i in 1:N){
    
    X <- matrix(c(1,runif(1,prmin,prmax),
                  0,0), byrow = TRUE, ncol=2)
    
    #discrete price levels
    #X <- matrix(c(1,PriceDraw(plvls),
    #             0,0), byrow = TRUE, ncol=2)
    
    # Calculation of the states.
    states.and.prior <- CreateStatesAndPrior(X, beta.vec2, rho, rho.prob)
    
    
    # RI Optimization with Shannon Costs
    choice.probs.output <- CalcChoiceProbsUnderRIWithShannon(
      Omega = states.and.prior[[1]],
      mu = states.and.prior[[2]],
      lambda = lambda, 
      max.iter = 10^7,
      precision = 10^(-10))
    
    state.dependent.choice.probs<-choice.probs.output[[2]]
    uncond.choice.probs<-choice.probs.output[[1]]
    
    eqm.mutual.info <- MutualInfo(prior = states.and.prior[[2]], 
                                  signal.dist = choice.probs.output[[1]], 
                                  posterior = choice.probs.output[[3]])
    
    # Draw the state of the world.
    curr.state <- sample(1:ncol(states.and.prior[[1]]), 1, replace = TRUE,
                         prob = states.and.prior[[2]]) 
    #curr.state <- sample(1:ncol(states.and.prior[[1]]), 1, replace = TRUE,
    #prob = rep(0.25,4))   
    
    # Combine current state and conditional choice probs to draw the choice
    choice <- sample(1:nrow(X), 1, replace = TRUE,
                     prob = choice.probs.output[[2]][, curr.state])
    # Output
    choice.sim <- list(X, choice, curr.state, states.and.prior, 
                       eqm.mutual.info, rho, rho.prob, # Regulation and its distribution
                       beta.vec2, lambda,state.dependent.choice.probs,
                       uncond.choice.probs) 
    names(choice.sim) <- c("X", "choice", "curr.state", "states.and.prior", 
                           "eqm.mutual.info",  "rho", "rho.prob", "beta.vec", "lambda",
                           "state.dependent.choice.probs","uncond.choice.probs")
    
    
    
    
    
    #   -----------------------------------------------------------------------
    
    # print(paste0("simulation number ", i))
    sim.list1[[j]][[i]] <- choice.sim
    
  }
  
  
  t2 <- Sys.time()
  
  t2 - t1
  
  
  
  ############################## Estimation #################################
  ###########################################################################
  
  
  # obfuscated price bonuses at each state  
  num.alternatives = 2 #dim(sim.list1[[1]]$X)[1]
  reg.states = expand.grid(rep(list(rho), num.alternatives - 1 ))
  reg.states = rbind(t(reg.states),0)
  
  
  ## Model 2: MNL with both the observable price and the separate obfuscated price 
  ## (bonus) as attributes
  X = NULL
  y = rep(0,N)
  for(t in 1:N){
    
    #add the bonus to the design 
    adj.design = cbind(sim.list1[[j]][[t]]$X, reg.states[,sim.list1[[j]][[t]]$curr.state])
    colnames(adj.design)[3] = "bonus"
    
    #storing individual data in y and X
    X = rbind(X,adj.design)
    y[t] = sim.list1[[j]][[t]]$choice 
    
  }
  
  Data2 = list(X = X, y = y, p = num.alternatives)
  Mcmc2 = list(R = 10000)
  out_double_price1 = rmnlIndepMetrop_mod(Data=Data2, Mcmc=Mcmc2, posttune=TRUE)
  #out_double_price1 = rmnlIndepMetrop(Data=Data2, Mcmc=Mcmc2)
  
  
  
  betadraw <- colMeans(out_double_price1$betadraw)
  names(betadraw) <- c("brand", "price", "bonus")
  logliks[j,]<-out_double_price1$loglike
  #draws<-table(y)
  
  
  SE<-c(sqrt(var(out_double_price1$betadraw[,1])),
        sqrt(var(out_double_price1$betadraw[,2])),
        sqrt(var(out_double_price1$betadraw[,3])))
  
  
  
  #Bonus levels: -2,0
  betas[j,]<-round(betadraw,2)
  standarderrors[j,]<-round(SE,2)
  
}


plot(betas[,1])













g<-list("Discount"=rhos[,1]-rhos[,2],"Beta_b"=betas[,1],"Beta_p"=betas[,2],
        "Beta_d"=betas[,3])
df<-data.frame(g)

#df<-df[1:150,] #adjust length of discounts to exclude zeros


test1<-ggplot(data=df, aes(x=Discount,y=Beta_b))+
  geom_point( size=1) + 
  ggtitle("Brand")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, size=12, face = "bold"))+
  #theme(plot.title = element_text(hjust = 0.5))+
  #scale_linetype_manual(values=c("solid","dashed"))+
  #scale_color_manual(values=c( "#F8766D","#00BFC4"))+
  # geom_line(aes(y =diff2, colour="Larger Range"), size=1) + 
  labs(x="Discount Range",y="Estimated Part Worth")+
coord_cartesian( xlim = c(0, 8))

test2<-ggplot(data=df, aes(x=Discount,y=Beta_p))+
  geom_point( size=1) + 
  ggtitle("Price")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, size=12, face = "bold"))+
  #theme(plot.title = element_text(hjust = 0.5))+
  #scale_linetype_manual(values=c("solid","dashed"))+
  #scale_color_manual(values=c( "#F8766D","#00BFC4"))+
  # geom_line(aes(y =diff2, colour="Larger Range"), size=1) + 
  labs(x="Discount Range",y="Estimated Part Worth")+
coord_cartesian( xlim = c(0, 8))


test3<-ggplot(data=df, aes(x=Discount,y=Beta_d))+
  geom_point( size=1) + 
  ggtitle("Discount")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, size=12, face = "bold"))+
  #theme(plot.title = element_text(hjust = 0.5))+
  #scale_linetype_manual(values=c("solid","dashed"))+
  #scale_color_manual(values=c( "#F8766D","#00BFC4"))+
  # geom_line(aes(y =diff2, colour="Larger Range"), size=1) + 
  labs(x="Discount Range",y="Estimated Part Worth")+
  coord_cartesian( xlim = c(0, 8))
windows(height=3, width=9)
par(mfrow=c(1,3))
grid.arrange(test1,test2,test3,ncol=3)

lls<-rep(0,nrow(logliks))
for (i in 1:length(lls)){
  lls[i]<-mean(logliks[i,])
}

#calculating mcfadden r?
#no of alternatives
nalt<-2
#no of lambda sims
NL<-length(sim.list1)
#no of simulated choces per lambda
N<-length(sim.list1[[1]])
probs<-array(rep(0,length(sim.list1[[1]])*nalt*length(sim.list1)),
             dim = c(length(sim.list1),length(sim.list1[[1]]),nalt))

ys<-matrix(rep(0,length(sim.list1)*length(sim.list1[[1]])),
           ncol=length(sim.list1[[1]]))

zwischenprobs<-matrix(rep(0,length(sim.list1)*length(sim.list1[[1]])),
                      ncol=length(sim.list1[[1]]))

loglkmod<-rep(0,length(sim.list1))


nullchoiceprobs<-matrix(rep(0,NL*nalt),ncol=nalt)

for (j in 1:length(probs[,1,1])){
  
  for (i in 1:length(probs[1,,1])){
    
    #Logit choice probs
    probs[j,i,1]<-exp(betas[j,]%*%c(sim.list1[[j]][[i]]$X[1,],
                                    sim.list1[[j]][[i]]$rho[sim.list1[[j]][[i]]$curr.state]))/
      (exp(betas[j,]%*%c(sim.list1[[j]][[i]]$X[1,],
                         sim.list1[[j]][[i]]$rho[sim.list1[[j]][[i]]$curr.state]))+1)
    probs[j,i,2]<- 1-probs[j,i,1]
    
    #choice matrix, j rows is grid of lambdas and i columns are choice simulations
    #per lambda
    ys[j,i]<-sim.list1[[j]][[i]]$choice
    
    #helpmatrix
    zwischenprobs[j,i]<-log(probs[j,i,ys[j,i]])
  }
  
  loglkmod[j]<-sum(zwischenprobs[j,])
  
  nullchoiceprobs[j,]<-c(sum(ys[j,]==1)/N, 1-sum(ys[j,]==1)/N)
}

loglknull<-rep(0,NL)
zwischennull<-matrix(rep(0,N*NL),ncol=N)
Mcfaddenrsquare<-rep(0,NL)

for (j in 1:NL){
  for (i in 1:N){
    zwischennull[j,i]<-log(nullchoiceprobs[j,ys[j,i]])
  } 
  loglknull[j]<-sum(zwischennull[j,])
  Mcfaddenrsquare[j]<-1-(loglkmod[j]/loglknull[j])
}



#Mcfad<-Mcfaddenrsquare[1:150]
Mcfad<-Mcfaddenrsquare
g2<-list("Discount"=df[,1], "Mcfad"=Mcfad)
df2<-data.frame(g2)

test4<-ggplot(data=df2, aes(x=Discount,y=Mcfad))+
  geom_point( size=1) + 
  ggtitle("Model Fit")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, size=12, face = "bold"))+
  #theme(plot.title = element_text(hjust = 0.5))+
  #scale_linetype_manual(values=c("solid","dashed"))+
  #scale_color_manual(values=c( "#F8766D","#00BFC4"))+
  # geom_line(aes(y =diff2, colour="Larger Range"), size=1) + 
  labs(x="Discount Range",y="Mc Fadden R squared")+
  coord_cartesian( ylim = c(0.25, 1))

windows(height=3, width=5.5)
par(mfrow=c(1,1))
test4

windows(height=3, width=9)
par(mfrow=c(1,3))
grid.arrange(test1,test2,test3,ncol=3)
