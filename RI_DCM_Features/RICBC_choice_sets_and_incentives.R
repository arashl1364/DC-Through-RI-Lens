# In this file, we 
# 1) create choice sets with a single inside good, a single outside good, and 
# a consumer who has preferences for the particular brand and some costs 
# of learning
# 2) Calculation of states based on a simple bonus scheme with three states. 

#  Change Log 
#  01.09.20: CreateStatesAndPrior() can be used with different price coef for 
#  the bonus
#  
#   

library("entropy")

########################################################################### 
# 
# Choice Sets -------------------------------------------------------------
# 
###########################################################################

# Simple choice set with single inside and single outside price exponential
CreateChoiceSetexp <- function(){
  # row = alternative, last row = outside option
  out <- matrix(c(1, rexp(1), 0, 0 ), byrow = TRUE, ncol = 2)
  colnames(out) <- c("brand", "price")
  rownames(out) <- c("alt 1", "outside")
  return(out)
} 

CreateChoiceSet <- function(plowbound,pupbound){
  # row = alternative, last row = outside option
  out <- matrix(c(1, runif(1,plowbound,pupbound), 
                  0, 0 ), byrow = TRUE, ncol = 2)
  colnames(out) <- c("brand", "price")
  rownames(out) <- c("alt 1", "outside")
  return(out)
} 
# 

# Choice set with two inside goods with individual brands / intercepts
# and exponentially distributed prices 

CreateChoiceSet2 <- function(){
  # row = alternative, last row = outside option
  out <- matrix(c(1, 0,  rexp(1),
                  0, 1,  rexp(1),
                  0, 0, 0), byrow = TRUE, ncol = 3)
  colnames(out) <- c(paste0("brand ", 1:(ncol(out) - 1)), 
                     "price")
  rownames(out) <- c(paste0("alt. ", 1:(ncol(out) - 1)), 
                     "outside")
  return(out)
} 

#uniform prices
CreateChoiceSet3 <- function(){
  # row = alternative, last row = outside option
  out <- matrix(c(1, 0,  runif(1,0.01,2),
                  0, 1,  runif(1,0.01,2),
                  0, 0, 0), byrow = TRUE, ncol = 3)
  colnames(out) <- c(paste0("brand ", 1:(ncol(out) - 1)), 
                     "price")
  rownames(out) <- c(paste0("alt. ", 1:(ncol(out) - 1)), 
                     "outside")
  return(out)
} 

CreateChoiceSet4 <- function(){
  # row = alternative, last row = outside option
  out <- matrix(c(1, 0, 0, runif(1,0.01,2),
                  0, 1, 0, runif(1,0.01,2),
                  0, 0, 1, runif(1,0.01,2),
                  0, 0, 0, 0), byrow = TRUE, ncol = 4)
  colnames(out) <- c(paste0("brand ", 1:(ncol(out) - 1)), 
                     "price")
  rownames(out) <- c(paste0("alt. ", 1:(ncol(out) - 1)), 
                     "outside")
  return(out)
} 

###########################################################################
# 
# Auxiliary functions -----------------------------------------------------
# 
###########################################################################
 
# Calculation of Learning Incentives --------------------------------------

IncentivesForFullInfo <- function(Omega, mu){
  value.prior.info <- max(Omega %*% mu)
  value.full.info <-  apply(Omega, 2, max) %*% mu
  incentives <- value.full.info - value.prior.info
  return(incentives)
}






# Calculation of state dependent payoffs ----------------------------------
# For the following, suppose that the regulation impact is independent across
# inside alternatives and suppose they are symmetric.


CreateStatesAndPrior <- function(X, beta.vec, 
                                 rho, rho.prob) {
  # This function creates a set of payoff vectors ("states") and their 
  # respective (prior) probabilities given a regulation (rho, rho.prob). 
  
  
  # Inputs # 
  # X: simple design
  # beta: preferences where last entry is complex attribute
  # rho, rho.prob: Regulation realizations and probs 
  
  # Outputs # 
  # Omega: Set of payoff vectors 
  # mu: "prior" probabilities
  
  util.index <- X %*% beta.vec[-length(beta.vec)]
  num.alternatives <- length(util.index)
  
  Omega <- matrix(util.index, 
                  ncol = length(rho)^(num.alternatives - 1), 
                  nrow = num.alternatives, 
                  byrow = FALSE)
  
  # Calculating the impact of regulation -----------------------------------------
  
  # output matrix
  # expand.grid creates all possible permutations 
  reg.states <- expand.grid(rep(list(rho), num.alternatives - 1 ))
  reg.states <- t(as.matrix(reg.states))
  reg.states <- rbind(reg.states, 0)
  Omega <- Omega + reg.states * beta.vec[length(beta.vec)]
  
  mu <- expand.grid(rep(list(rho.prob), num.alternatives - 1 ))
  mu <- t(as.matrix(mu))
  mu <- apply(mu, 2, prod)
  
  out <- list(Omega, mu)
  return(out)
}

CreateStatesAndPrior2 <- function(X, beta.vec, 
                                 rho, rho.prob) {
  # This function creates a set of payoff vectors ("states") and their 
  # respective (prior) probabilities given a regulation (rho, rho.prob). 
  
  
  # Inputs # 
  # X: Choice sets
  # beta: preferences
  # rho, rho.prob: Regulation realizations and probs 
  
  # Outputs # 
  # Omega: Set of payoff vectors 
  # mu: "prior" probabilities
  
  util.index <- X %*% beta.vec
  num.alternatives <- length(util.index)
  
  Omega <- matrix(util.index, 
                  ncol = length(rho)^(num.alternatives - 1), 
                  nrow = num.alternatives, 
                  byrow = FALSE)
  
  # Calculating the impact of regulation -----------------------------------------
  
  # output matrix
  # expand.grid creates all possible permutations 
  reg.states <- expand.grid(rep(list(rho), num.alternatives - 1 ))
  reg.states <- t(as.matrix(reg.states))
  reg.states <- rbind(reg.states, 0)
  Omega <- Omega + reg.states *(beta.vec[length(beta.vec)])
  
  mu <- expand.grid(rep(list(rho.prob), num.alternatives - 1 ))
  mu <- t(as.matrix(mu))
  mu <- apply(mu, 2, prod)
  
  out <- list(Omega, mu)
  return(out)
}


#no outside option
CreateStatesAndPriorNO <- function(X, beta.vec, 
                                 rho, rho.prob) {
  # This function creates a set of payoff vectors ("states") and their 
  # respective (prior) probabilities given a regulation (rho, rho.prob). 
  
  
  # Inputs # 
  # X: Choice sets
  # beta: preferences
  # rho, rho.prob: Regulation realizations and probs 
  
  # Outputs # 
  # Omega: Set of payoff vectors 
  # mu: "prior" probabilities
  
  util.index <- X %*% beta.vec[-length(beta.vec)]
  num.alternatives <- length(util.index)
  
  Omega <- matrix(util.index, 
                  ncol = length(rho)^(num.alternatives), 
                  nrow = num.alternatives, 
                  byrow = FALSE)
  
  # Calculating the impact of regulation -----------------------------------------
  
  # output matrix
  # expand.grid creates all possible permutations 
  reg.states <- expand.grid(rep(list(rho), num.alternatives))
  reg.states <- t(as.matrix(reg.states))
  #reg.states <- rbind(reg.states, 0)
  Omega <- Omega + reg.states * beta.vec[length(beta.vec)]
  
  mu <- expand.grid(rep(list(rho.prob), num.alternatives ))
  mu <- t(as.matrix(mu))
  mu <- apply(mu, 2, prod)
  
  out <- list(Omega, mu)
  return(out)
}


CreateStatesAndPriorCORR <- function(X, beta.vec, 
                                   rho, rho.prob, theta) {
  # This function creates a set of payoff vectors ("states") and their 
  # respective (prior) probabilities given a regulation (rho, rho.prob). 
  
  
  # Inputs # 
  # X: Choice sets
  # beta: preferences
  # rho, rho.prob: Regulation realizations and probs 
  
  # Outputs # 
  # Omega: Set of payoff vectors 
  # mu: "prior" probabilities
  
  util.index <- X %*% beta.vec
  num.alternatives <- length(util.index)
  
  Omega <- matrix(util.index, 
                  ncol = length(rho)^(num.alternatives), 
                  nrow = num.alternatives, 
                  byrow = FALSE)
  
  # Calculating the impact of regulation -----------------------------------------
  
  # output matrix
  # expand.grid creates all possible permutations 
  reg.states <- expand.grid(rep(list(rho), num.alternatives))
  reg.states <- t(as.matrix(reg.states))
  #reg.states <- rbind(reg.states, 0)
  Omega <- Omega + reg.states * beta.vec[length(beta.vec)]
  
  mu <- expand.grid(rep(list(rho.prob), num.alternatives ))
  mu <- t(as.matrix(mu))
  mu <- apply(mu, 2, prod)
  
  out <- list(Omega, mu)
  return(out)
}

#   -----------------------------------------------------------------------
# The following function CAN be used in order to save some calculation time 
# by precalculating once the impact of regulation 
CreateRegulationImpact <- function(rho, rho.prob, N){
  # rho, rho.prob: Bonus and its distribution for a single alternative 
  # N: number of alternatives in the choice set that can be independently (!!!)
  # affected by the regulation
  # There is an outside option
  
  # output matrix
  # expand.grid creates all possible permutations 
  
  reg.states <- expand.grid(rep(list(rho),N ))
  reg.states <- t(as.matrix(reg.states))
  reg.states <- rbind(reg.states, 0)
  
  reg.mu <- expand.grid(rep(list(rho.prob),N ))
  reg.mu <- t(as.matrix(reg.mu))
  reg.mu <- apply(reg.mu, 2, prod)
  
  out <- list(reg.states, reg.mu)
  
  names(out) <- c("reg.states", "reg.mu")
  
  
  return(out)
}



###########################################################################
# 
# Measures for the Impact of Learning -------------------------------------
# 
###########################################################################

# Difference between the unconditional and conditional probs --------------
DiffBetweenConditionalUnconditional <- function(uncond.probs, cond.probs){
  # Calculates, given optimal conditional and unconditional choice probs from 
  # RI choice, the "difference" between those probs as simple average. 
  # 
  # uncond.probs: vector of choice probs, lenght = number of actions 
  # cond.probs: choice probs 
  out <- sweep(cond.probs, 1, uncond.probs, "-" )
  out <- sum(abs(out)) / dim(cond.probs)[2]
  return(out)
}


# Difference between Prior and Posterior Probs ----------------------------

DiffPriorDecisionAndPosteriorDecision <- function(Omega, mu, cond.probs){
  # This calculates the changes in choice probs due to learning 
  # Prior probs are based on the prior information without any learning
  # These are NOT the unconditional choice probs!
  optimal.prior.choice <- which.max(Omega %*% mu)
  prior.choice.probs <- rep(0, dim(Omega)[1]) 
  prior.choice.probs[optimal.prior.choice] <- 1
  
  out <- sweep(cond.probs, 1, prior.choice.probs, "-" )
  out <- sum(abs(out)) / dim(cond.probs)[2]
  return(out)
}


# Equilibrium Mutual Information ------------------------------------------

MutualInfo <- function(prior, signal.dist, posterior){
  # prior: Prior distribution over states
  # signal.dist: distribution of observable states. Here: Signal = Choice Action 
  # posterior: Matrix of n.actions x n.states with posteriors cond. on action
  
  post.entropies <- apply(posterior, 1, entropy)
  expected.entropy <- signal.dist %*% post.entropies
  out <- entropy(prior) - expected.entropy 
  return(out)
}

MutualInfo2 <- function(prior, signal.dist, posterior){
  # prior: Prior distribution over states
  # signal.dist: distribution of observable states. Here: Signal = Choice Action 
  # posterior: Matrix of n.actions x n.states with posteriors cond. on action
  
  post.entropies <- apply(posterior, 1, entropy)
  expected.entropy <- signal.dist %*% post.entropies
  out <- entropy(prior) - expected.entropy 
  return(list(out,entropy(prior),expected.entropy))
}














