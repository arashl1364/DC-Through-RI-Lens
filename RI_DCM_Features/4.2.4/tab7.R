
rm(list = ls())
setwd("C:/Users/Matteo/Desktop/RI/Plots and Tables")

#p in 2 to 3
#set.seed(66)
#set.seed(666)

#p in 1 to 3
#set.seed(66)
#set.seed(10)


source("RICBC_choice_sets_and_incentives.R")
source("RICBC_RI_SC_Choice.R")

set.seed(66)
#set.seed(777)

beta.vec <- c(3, -1,1) 
beta.vec2 <- c(3, 3, -1,1)



lambda <- 0.5 
rho <- c(0,2) #Possible Discounts
rho.prob <- rep(1, length(rho))/length(rho)

# Simulation --------------------------------------------------------------
N=1000
sim.list <- list()
length(sim.list) <- N
sim.list2 <- list()
length(sim.list2) <- N

for(i in 1:N){
  
  #fixing our design matrix
  X <- matrix(c(1,runif(1,2,4), #0,2
                0,0), byrow = TRUE, ncol = 2)
  #X <- matrix(c(1,rexp(1),
  #              0,0), byrow = TRUE, ncol = 2)
  #X <- CreateChoiceSet()
  # Calculation of the states.
  states.and.prior <- CreateStatesAndPrior(X, beta.vec, rho, rho.prob)
  
  # RI Optimization with Shannon Costs
  choice.probs.output <- CalcChoiceProbsUnderRIWithShannon(
    Omega = states.and.prior[[1]],
    mu = states.and.prior[[2]],
    lambda = lambda,
    max.iter = 10^7,
    precision = 10^(-10)     
  )      
  
  state.dependent.choice.probs<-choice.probs.output[[2]]
  uncond.choice.probs<-choice.probs.output[[1]]
  
  
  # Mutual Information as a measure of how much is learned 
  eqm.mutual.info <- MutualInfo2(prior = states.and.prior[[2]], 
                                 signal.dist = choice.probs.output[[1]], 
                                 posterior = choice.probs.output[[3]])
  
  
  
  # Simulating choice  ------------------------------------------------------
  # Draw the state of the world.
  curr.state <- sample(1:ncol(states.and.prior[[1]]), 1, replace = TRUE,
                       prob = states.and.prior[[2]])   
  
  # Combine current state and conditional choice probs to draw the choice
  choice <- sample(1:nrow(X), 1, replace = TRUE,
                   prob = choice.probs.output[[2]][, curr.state])
  # Output
  choice.sim <- list(X, choice, curr.state, states.and.prior, 
                     eqm.mutual.info,
                     rho, rho.prob, # Regulation and its distribution
                     beta.vec, lambda,state.dependent.choice.probs,
                     uncond.choice.probs)  # Preferences and learning costs
  names(choice.sim) <- c("X", "choice", "curr.state", "states.and.prior", 
                         "eqm.mutual.info",
                         "rho", "rho.prob", "beta.vec", "lambda",
                         "state.dependent.choice.probs","uncond.choice.probs")
  
  
  
  print(paste0("simulation number ", i))
  sim.list[[i]] <- choice.sim
  
}



###########################################################################
########################### Logit Estimation ##############################
###########################################################################
library(bayesm)
# obfuscated price bonuses at each state  
num.alternatives = dim(sim.list[[1]]$X)[1]
reg.states = expand.grid(rep(list(rho), num.alternatives - 1 ))
reg.states = rbind(t(reg.states),0)

## Model 2: MNL with both the observable price and the separate obfuscated price 
## (bonus) as attributes
X = NULL
y = rep(0,N)
for(i in 1:N){
  
  #add the bonus to the design 
  adj.design = cbind(sim.list[[i]]$X, reg.states[,sim.list[[i]]$curr.state])
  colnames(adj.design)[3] = "bonus"
  
  #storing individual data in y and X
  X = rbind(X,adj.design)
  y[i] = sim.list[[i]]$choice 
  
}

Data2 = list(X = X, y = y, p = num.alternatives)
Mcmc2 = list(R = 10000)
out_double_price = rmnlIndepMetrop(Data=Data2, Mcmc=Mcmc2)


betadraw1 <- colMeans(out_double_price$betadraw)
names(betadraw1) <- c("brand1", "price", "bonus")
draws1<-table(y)



sim.list1<-sim.list
nalt<-2
#no of simulated choces
N<-length(sim.list1)

probs<-matrix(rep(0,N*nalt),ncol=nalt)

ys<-c()
#ys<-rep(0,N)

zwischenprobs<-c()

#zwischenprobs<-matrix(rep(0,length(sim.list1)*length(sim.list1[[1]])),
#                      ncol=length(sim.list1[[1]]))

#loglkmod<-rep(0,length(sim.list1))

nullchoiceprobs<-c()
#nullchoiceprobs<-matrix(rep(0,NL*nalt),ncol=nalt)



for (i in 1:N){
  
  #Logit choice probs
  probs[i,1]<-exp(betadraw1%*%c(sim.list1[[i]]$X[1,],
                                sim.list1[[i]]$rho[sim.list1[[i]]$curr.state]))/
    (exp(betadraw1%*%c(sim.list1[[i]]$X[1,],
                       sim.list1[[i]]$rho[sim.list1[[i]]$curr.state]))+1)
  probs[i,2]<- 1-probs[i,1]
  
  #choice matrix, j rows is grid of lambdas and i columns are choice simulations
  #per lambda
  ys[i]<-sim.list1[[i]]$choice
  
  
  
  #helpmatrix
  zwischenprobs[i]<-log(probs[i,ys[i]])
}
#logit loglik
loglkmod<-sum(zwischenprobs)
nullchoiceprobs<-c(sum(ys==1)/N, 1-sum(ys==1)/N)



zwischennull<-c()

for (i in 1:N){
  zwischennull[i]<-log(nullchoiceprobs[ys[i]])
} 
loglknull<-sum(zwischennull)
Mcfaddenrsquare<-1-(loglkmod/loglknull)






###############################2 inside alternatives ##################################

#collect vector of previous prices to simulate choice set expansions:

prices<-rep(0,length(sim.list))
for (i in 1:length(prices)){
  prices[i]<-sim.list[[i]]$X[1,2]
}





t1 <- Sys.time()
for(i in 1:N){
  
  #X <- matrix(c(1,0,rexp(1),
  #              0,1,rexp(1),
  #              0,0,0), byrow = TRUE,ncol=3)
  
  X <- matrix(c(1,0,runif(1,2,4), #0,2
                0,1,runif(1,2,4),
                0,0,0), byrow = TRUE,ncol=3)
  #X<-CreateChoiceSet3()
  #X[1,3]<-prices[i]
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
  
  eqm.mutual.info <- MutualInfo2(prior = states.and.prior[[2]], 
                                 signal.dist = choice.probs.output[[1]], 
                                 posterior = choice.probs.output[[3]])
  
  # Draw the state of the world.
  curr.state <- sample(1:ncol(states.and.prior[[1]]), 1, replace = TRUE,
                       prob = states.and.prior[[2]])   
  #prob=rep(0.25,4) ) #manipulate prior 
  #
  # Combine current state and conditional choice probs to draw the choice
  choice <- sample(1:nrow(X), 1, replace = TRUE,
                   prob = choice.probs.output[[2]][, curr.state])
  # Output
  choice.sim <- list(X, choice, curr.state, states.and.prior, 
                     eqm.mutual.info,
                     rho, rho.prob, # Regulation and its distribution
                     beta.vec2, lambda,state.dependent.choice.probs,
                     uncond.choice.probs) 
  names(choice.sim) <- c("X", "choice", "curr.state", "states.and.prior", 
                         "eqm.mutual.info",
                         "rho", "rho.prob", "beta.vec", "lambda",
                         "state.dependent.choice.probs","uncond.choice.probs")
  
  
  
  
  
  #   -----------------------------------------------------------------------
  
  print(paste0("simulation number ", i))
  sim.list2[[i]] <- choice.sim
  
}


t2 <- Sys.time()

t2 - t1

prices2<-rep(0,length(sim.list2))
for (i in 1:length(prices2)){
  prices2[i]<-sim.list2[[i]]$X[2,3]
}

############################## Estimation #################################
###########################################################################


# obfuscated price bonuses at each state  
num.alternatives = dim(sim.list2[[1]]$X)[1]
reg.states = expand.grid(rep(list(rho), num.alternatives - 1 ))
reg.states = rbind(t(reg.states),0)


X = NULL
y = rep(0,N)
for(i in 1:N){
  
  #add the bonus to the design 
  adj.design = cbind(sim.list2[[i]]$X, reg.states[,sim.list2[[i]]$curr.state])
  colnames(adj.design)[4] = "bonus"
  
  #storing individual data in y and X
  X = rbind(X,adj.design)
  y[i] = sim.list2[[i]]$choice 
  
}

Data2 = list(X = X, y = y, p = num.alternatives)
Mcmc2 = list(R = 10000)
out_double_price2 = rmnlIndepMetrop(Data=Data2, Mcmc=Mcmc2)



betadraw2 <- colMeans(out_double_price2$betadraw)
names(betadraw2) <- c("brand1", "brand2", "price", "bonus")
draws2<-table(y)







SE1<-c(sqrt(var(out_double_price$betadraw[,1])), 
       sqrt(var(out_double_price$betadraw[,2])), 
       sqrt(var(out_double_price$betadraw[,3]))
)

SE2<-c(sqrt(var(out_double_price2$betadraw[,1])), 
       sqrt(var(out_double_price2$betadraw[,2])), 
       sqrt(var(out_double_price2$betadraw[,3])),
       sqrt(var(out_double_price2$betadraw[,4]))
)

#mean MCMC draws one inside alternative
round(betadraw1,2)
round(SE1,2)

#mean MCMC draws two inside alternatives
round(betadraw2,2)
round(SE2,2)


####Calculate MC Fadden R squared


sim.list1<-sim.list
nalt<-2
betas<-betadraw1
#no of simulated choces
N<-length(sim.list1)

probs<-matrix(rep(0,N*nalt),ncol=nalt)

ys<-c()
#ys<-rep(0,N)

zwischenprobs<-c()

#zwischenprobs<-matrix(rep(0,length(sim.list1)*length(sim.list1[[1]])),
#                      ncol=length(sim.list1[[1]]))

#loglkmod<-rep(0,length(sim.list1))

nullchoiceprobs<-c()
#nullchoiceprobs<-matrix(rep(0,NL*nalt),ncol=nalt)



for (i in 1:N){
  
  #Logit choice probs
  probs[i,1]<-exp(betas%*%c(sim.list1[[i]]$X[1,],
                            sim.list1[[i]]$rho[sim.list1[[i]]$curr.state]))/
    (exp(betas%*%c(sim.list1[[i]]$X[1,],
                   sim.list1[[i]]$rho[sim.list1[[i]]$curr.state]))+1)
  probs[i,2]<- 1-probs[i,1]
  
  #choice matrix, j rows is grid of lambdas and i columns are choice simulations
  #per lambda
  ys[i]<-sim.list1[[i]]$choice
  
  
  
  #helpmatrix
  zwischenprobs[i]<-log(probs[i,ys[i]])
}
#logit loglik
loglkmod<-sum(zwischenprobs)
nullchoiceprobs<-c(sum(ys==1)/N, 1-sum(ys==1)/N)



zwischennull<-c()

for (i in 1:N){
  zwischennull[i]<-log(nullchoiceprobs[ys[i]])
} 
loglknull<-sum(zwischennull)
Mcfaddenrsquare<-1-(loglkmod/loglknull)



sim.list1<-sim.list2
nalt<-3
betas<-betadraw2
#no of simulated choces
N<-length(sim.list1)

probs<-matrix(rep(0,N*nalt),ncol=nalt)

ys<-c()
#ys<-rep(0,N)

zwischenprobs<-c()

#zwischenprobs<-matrix(rep(0,length(sim.list1)*length(sim.list1[[1]])),
#                      ncol=length(sim.list1[[1]]))

#loglkmod<-rep(0,length(sim.list1))

nullchoiceprobs<-c()
#nullchoiceprobs<-matrix(rep(0,NL*nalt),ncol=nalt)

#states with two inside goods have 4 realizations, we need to map the state
#to a bonus level (1 or 0):
a<-c() #A1
b<-c() #A2
c<-0 #outside

for (i in 1:N){
  if (sim.list1[[i]]$curr.state == 1 || sim.list1[[i]]$curr.state == 3)
  {a[i] = 0
  } else{
    a[i] = 2}
  
  if (sim.list1[[i]]$curr.state == 1 || sim.list1[[i]]$curr.state == 2)
  {b[i] = 0
  } else{
    b[i] = 2}
  
}

for (i in 1:N){
  
  #Logit choice probs
  probs[i,1]<-exp(betas%*%c(sim.list1[[i]]$X[1,],
                            a[i]))/
    ( exp(betas%*%c(sim.list1[[i]]$X[1,],
                    a[i]))+
        exp(betas%*%c(sim.list1[[i]]$X[2,],
                      b[i]))+
        exp(betas%*%c(sim.list1[[i]]$X[3,],
                      c))
    )
  
  probs[i,2]<- exp(betas%*%c(sim.list1[[i]]$X[2,],
                             b[i]))/
    ( exp(betas%*%c(sim.list1[[i]]$X[1,],
                    a[i]))+
        exp(betas%*%c(sim.list1[[i]]$X[2,],
                      b[i]))+
        exp(betas%*%c(sim.list1[[i]]$X[3,],
                      c))
    )
  probs[i,3]<-exp(betas%*%c(sim.list1[[i]]$X[3,],
                            c))/
    ( exp(betas%*%c(sim.list1[[i]]$X[1,],
                    a[i]))+
        exp(betas%*%c(sim.list1[[i]]$X[2,],
                      b[i]))+
        exp(betas%*%c(sim.list1[[i]]$X[3,],
                      c))
    )
  #1-probs[i,1]-probs[i,2]
  
  #choice matrix, j rows is grid of lambdas and i columns are choice simulations
  #per lambda
  ys[i]<-sim.list1[[i]]$choice
  
  
  
  #helpmatrix
  zwischenprobs[i]<-log(probs[i,ys[i]])
}
#logit loglik
loglkmod<-sum(zwischenprobs)

nullchoiceprobs<-c(sum(ys==1)/N, sum(ys==2)/N, 1-sum(ys==1)/N-sum(ys==2)/N)



zwischennull<-c()

for (i in 1:N){
  zwischennull[i]<-log(nullchoiceprobs[ys[i]])
} 
loglknull<-sum(zwischennull)
Mcfaddenrsquare2<-1-(loglkmod/loglknull)

#####Results
#mean MCMC draws one inside alternative
round(betadraw1,2)
round(SE1,2)
round(Mcfaddenrsquare,2)

#mean MCMC draws two inside alternatives
round(betadraw2,2)
round(SE2,2)
round(Mcfaddenrsquare2,2)

