rm(list = ls())
setwd("C:/Users/Matteo/Desktop/RI/Plots and Tables")

#seed to replicate exact numbers from paper
set.seed(66)



source("RICBC_choice_sets_and_incentives.R")
source("RICBC_RI_SC_Choice.R")


beta.vec <- c(2, -1, 1)

lambda <- .5 #information processing cost
plow=2 #lower price limit
pup=4  #upper price limit
rho <- c(0,2) #complex attribute levels
rho.prob <- rep(1, length(rho))/length(rho) #distribution of complex attribute levels

# Simulation --------------------------------------------------------------
N=5000 #number of simulated choices
sim.list <- list()
length(sim.list) <- N

for(i in 1:N){
  #functions can be found in "RICBC_choice_sets_and_incentives.R"
  X <- CreateChoiceSet(plow,pup) 
  
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
  eqm.mutual.info <- MutualInfo(prior = states.and.prior[[2]], 
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


####################################################################################
############################## ESTIMATION ########################################
####################################################################################

library(bayesm)

num.alternatives = dim(sim.list[[1]]$X)[1]
reg.states = expand.grid(rep(list(rho), num.alternatives - 1 ))
reg.states = rbind(t(reg.states),0)


####################################################################################
############### MODEL WITHOUT INTERACTIONS #########################################

X = NULL
y = rep(0,N)
for(i in 1:N){
  
  # add the complex discount and the interaction terms to the design
  # 
  adj.design = cbind(sim.list[[i]]$X,
                     reg.states[,sim.list[[i]]$curr.state]
  )
  colnames(adj.design)[3] = "complex discount"
  
  
  
  
  #storing individual data in y and X
  X = rbind(X,adj.design)
  y[i] = sim.list[[i]]$choice 
  
}

#Data  contains data for model without interaction terms
Data = list(X = X, y = y, p = num.alternatives)
Mcmc = list(R = 20000)
out_interaction = rmnlIndepMetrop(Data=Data, Mcmc=Mcmc)

beta.interaction <- colMeans(out_interaction$betadraw)
#names(beta.interaction) <- c("brand", "price", "bonus", "interaction discount x price",
# "interaction discount x price²")


####################################################################################
############### MODEL WITHOUT INTERACTIONS #########################################
X = NULL
y = rep(0,N)
for(i in 1:N){
  
  # add the complex discount and the interaction terms to the design
  # 
  adj.design = cbind(sim.list[[i]]$X,
                     reg.states[,sim.list[[i]]$curr.state]
                     
                     ,c(reg.states[1,sim.list[[i]]$curr.state]*sim.list[[i]]$X[1,2],
                      0)
                     
                     ,c((reg.states[1,sim.list[[i]]$curr.state])*(sim.list[[i]]$X[1,2])^2,
                       0)
  )
  colnames(adj.design)[3] = "complex discount"
  colnames(adj.design)[4] = "complex discount x price"
  colnames(adj.design)[5] = "complex discount x price²"
  
  
  
  #storing individual data in y and X
  X = rbind(X,adj.design)
  y[i] = sim.list[[i]]$choice 
  
}

#Data2  contains data for model with interaction terms
Data2 = list(X = X, y = y, p = num.alternatives)
Mcmc2 = list(R = 20000)
out_interaction2 = rmnlIndepMetrop(Data=Data2, Mcmc=Mcmc2)

beta.interaction2 <- colMeans(out_interaction2$betadraw)
#names(beta.interaction) <- c("brand", "price", "bonus", "interaction discount x price",
# "interaction discount x price²")


################################# RESULTS #############################################

summary.bayesm.mat(out_interaction$betadraw)
colMeans(out_interaction$betadraw[2001:20000,])
summary.bayesm.mat(out_interaction2$betadraw)
colMeans(out_interaction2$betadraw[2001:20000,])



#traceplots MCMC draws
#windows()
#par(mfrow=c(3,2))
#plot(out_interaction$betadraw[,1], type="l")
#plot(out_interaction$betadraw[,2], type="l")
#plot(out_interaction$betadraw[,3], type="l")
#plot(out_interaction$betadraw[,4], type="l")
#plot(out_interaction$betadraw[,5], type="l")
