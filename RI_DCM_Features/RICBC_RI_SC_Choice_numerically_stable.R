# In this file, I implement a function that solves numerically for the optimal 
# choice probabilities with the Blahut-Arimoto algorithm. 
# These can be used as an input to determine choices of a RI SC DM.
# See Caplin, Dean, and Leahy (2019) "Rational inattention, 
# optimal consideration sets, and stochastic choice", p. 1067ff. 
# for details and formulas. 
# Author: Sergey Turlo, 12 July 2020
#         Matteo Fina, 2021
# CHANGES
# 
# 17 July 20: Use of Bhattacharyya distance as the convergence criterion. 
# 22 July 20: Change of the starting values (for the iteration)

# TO-Do: 
# 
# 1) Implement a check that looks at the conditions as suggested in Caplin et al. 
# and ensures that the solution actually satisfies them. 

BhatDistance <- function(p, q){
  # Bhattacharyya distance between two discrete probabilities on the same state
  # space. 
  out <- - log(sum(sqrt(p * q))) 
  return(out)
}





CalcChoiceProbsUnderRIWithShannon_numstable <- function(Omega, mu, lambda, 
                                              max.iter = 10^8, 
                                              precision = 10^(-16))
{
  
  # Description -------------------------------------------------------------
  
  # Omega is the matrix that describes the payoffs for each action in each state 
  # Columns = States, Rows = Actions 
  
  # mu: Prior distribution over States
  
  # lambda: Flow cost of learning, must be strictly positive! 
  # If too large(>20), adjust max.iter and precision! 
  
  # max.iter: maximal number of iterations 
  
  # precision: threshold for the sum of squared differences between the choice 
  # probabilities of two subsequent iterations so that iterations stop
  
  
  # Pre-Calculations --------------------------------------------------------
  nalt=nrow(Omega)
  nstates=ncol(Omega)
  maxl=c()
  Z_mat=matrix(nrow=nalt,ncol=nstates)
  G=Omega/lambda
  for (i in 1:ncol(G)){
    maxl[i]=max(G[,i])
    Z_mat[,i] <- exp(G[,i]-maxl[i])
    
  }

  # Transformation of the utility index matrix
  #Z_mat <- exp(Omega/lambda)
  
  N_actions <- dim(Omega)[1]
  N_states <- dim(Omega)[2]
  # CREATE A CHECK THAT DIMENSIONS ALL FIT!!!!!!!
  
  
  # Starting Guess ----------------------------------------------------------
  # Pick as a starting guess for the unconditional choice probs a distribution 
  # where the ex-ante optimal action is chosen with probability  
  # "optimal.action.weight". The objective is a better performance with high 
  # lambda values. 
  # 


  optimal.action.weight <- 0.75

  optimal.prior.choice <- which.max(Omega %*% mu)
  Start_Probs <- rep(0, N_actions)
  Start_Probs[optimal.prior.choice] <- optimal.action.weight
  Start_Probs[- optimal.prior.choice] <-
    (1 - optimal.action.weight) / ( N_actions -1)

  Start_CondProbs <- matrix(rep(Start_Probs, N_states), ncol = N_states)
  
  
  
  # Start_CondProbs <- matrix(1/N_actions, nrow = N_actions, ncol = N_states)
  # Start_Probs <- rep(1/N_actions, N_actions)
  
  
  
  # Auxiliary functions used for the iterations -----------------------------
  # Step 1: Calculation of conditional probabilities
  .step1 <- function(Current_Probs, Current_CondProbs){
    new_CondProbs <- Current_CondProbs  
    for(i in 1:N_actions){
      for(j in 1:N_states){
        D1 <- maxl[j]+ log(Current_Probs) %*% Z_mat[, j]
        new_CondProbs[i, j] <- log(Current_Probs[i]) + G[i,j] - D1 
        new_CondProbs[i, j] <-exp(new_CondProbs[i, j])
      }
    }
    return(new_CondProbs)
  }
  
  # Step 2: Calculation of actual (prior) choice probabilities
  .step2 <- function(Current_CondProbs){
    new_Probs <-Current_CondProbs %*% mu
    return(c(new_Probs))
  }
  
  # .iter: Takes a list that contains the current guess of 1) the choice probs 
  # and 2) the choice probs conditional on the state and calculates a new guess  
  # as its output. 
  
  .iter <- function(prob_list){
    old_Probs <- prob_list[[1]]
    old_CondProbs <- prob_list[[2]]
    
    new_CondProbs <-
      .step1(Current_Probs = old_Probs, Current_CondProbs = old_CondProbs)
    new_Probs <- .step2(Current_CondProbs = new_CondProbs)
    
    out <- list(new_Probs, new_CondProbs)
    return(out)
  }
  
  # Actual Iterations -------------------------------------------------------
  
  delta <-  1
  counter <- 0
  prob_list <- list(Start_Probs, Start_CondProbs)
  
  while(delta > precision && counter < max.iter){
    old_prob_list <- prob_list
    prob_list <- .iter(prob_list)
    # delta <- sum((old_prob_list[[1]] - prob_list[[1]])^2)
    # Bhattacharyya distance between old_prob_list[[1]] and prob_list[[1]]
    delta <- BhatDistance(old_prob_list[[1]], prob_list[[1]])    
    
    counter <- counter + 1
    
  }
  
  # print(counter)
  # print(delta)
  
  # print(paste0("Iterations until convergence: ", counter, "."))
  if(counter >= max.iter){
    print("NO CONVERGENCE!")
  }
  
  
  # Revealed Posterior Beliefs ----------------------------------------------
  post_beliefs <- sweep(prob_list[[2]], 1, prob_list[[1]], FUN = '/')
  post_beliefs <- sweep(post_beliefs, 2, mu, FUN = '*')
  
  prob_list[[3]] <- post_beliefs
  
  names(prob_list) <- 
    c("Choice Probabilities", "State Dependent Choice Probabilities", 
      "Revealed Posterior Beliefs")
  
  names(prob_list[[1]]) <- paste0(rep("Action ", N_actions), 1:N_actions)
  rownames(prob_list[[2]]) <- paste0(rep("Action ", N_actions), 1:N_actions, "|")
  colnames(prob_list[[2]]) <- paste0(rep("State ", N_states), 1:N_states)
  
  
  rownames(prob_list[[3]]) <- paste0(rep("Action ", N_actions), 1:N_actions, "|")
  colnames(prob_list[[3]]) <- paste0(rep("State ", N_states), 1:N_states)
  
  
  # Rounding of results (IS THIS OKAY? ESPECIALLY AS WE ARE INTERESTED 
  # IN Pr(a) = 0)! -------------------------------------
  prob_list[[1]] <- round(prob_list[[1]], 5)
  prob_list[[2]] <- round(prob_list[[2]], 5)
  prob_list[[3]] <- round(prob_list[[3]], 5)
  
  return(prob_list) 
}



















