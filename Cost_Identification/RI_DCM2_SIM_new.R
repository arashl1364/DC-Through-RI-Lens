RI_DCM2_Sim_new = function(beta, lambda, rho, simplelvls, N,comp_idx, brand=FALSE){
  #Arash Laghaie 2020
  #Matteo Fina 2022
  
  # Simulate choices from RI DCM with 1 inside 1 outside and 1 complex attribute
  #Rcpp::sourceCpp('BA_updated.cpp')
  
  nvar = length(beta) 
  nalt = 2
  ncomplevels = length(rho)  # number of complex attribute levels 
  rho.prob = rep(1,ncomplevels)/ncomplevels   
  sim.list <- list()
  length(sim.list) <- N
  nstates = ncomplevels^(nalt-1)  # number of all possible states
  
  # Calculating mu (prior over states) assuming a uniform prior
  mu = rep(1,nstates)/(nstates)
  
  # Calculating all possible realization of bonus schemes 
  states = expand.grid(rep(list(rho), nalt - 1 ))
  states = t(as.matrix(states))
  states = rbind(states, 0)
  row.names(states) = c("alt1","outside")
  
  for(i in 1:N){
    #
    # Creating the design matrix
    #
    # continuous simple price
    # X = rbind(cbind(diag(nalt-1),rexp(nalt-1)), rep(0,nalt)) # last row of X is for the outside good with all its elements zero
    ##### here needs to be generalized ####
    # simple attributes
     # if price is simple 
     # discount = simplelvls
    if (brand==TRUE){
    Xsimp = rbind(cbind(diag(nalt-1),sample(simplelvls,nalt-1,T)), 0) # last row of X is for the outside good with all its elements zero
      
    }
    else (Xsimp = rbind(cbind(sample(simplelvls,nalt-1,T)), 0) )
    
    
    
    #
    # Calculating cond. and uncond. choice probs using BA algorithm
    #
    out = BA_upd(Xsimp, beta[-comp_idx], lambda, mu, states, N, beta[comp_idx])
    
    
    # adding the realized complex attribute to the design 
    # Draw the state of the world.
    curr_state <- sample(1:length(mu), 1, prob = mu)
    X = cbind(Xsimp, states[,curr_state])
    
    
    #
    # Combine current state and conditional choice probs to draw the choice
    #
    choice = sample(1:nrow(X), 1, replace = TRUE,
                    prob = out$condprobs[, curr_state])
    
    # Output
    choice.sim <- list(X, choice, curr_state, out$Omega, mu, 
                       rho, rho.prob, # Regulation and its distribution
                       beta, lambda, out)  # Preferences and learning costs
    names(choice.sim) <- c("X", "choice", "curr_state", "Omega", "mu", 
                           "rho", "rho.prob", "beta.vec", "lambda", "BAoutput")
    
    #   -----------------------------------------------------------------------
    
    sim.list[[i]] <- choice.sim
    
  }  
  return(sim.list)
}

