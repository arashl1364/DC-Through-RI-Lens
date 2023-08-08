RIDC2_Sim = function(beta, lambda, rho,simplelevels, N,complex){

  nvar = length(beta)   
  nalt = 2 #nalt fixed to 2
  ncomplevels = length(rho)  # number of complex attribute levels 
  rho.prob = rep(1,ncomplevels)/ncomplevels  #distribution of attr levels 
  sim.list <- list()
  length(sim.list) <- N
  nstates = ncomplevels^(nalt-1)  # number of all possible states
  
  # Calculating mu (prior over states) assuming a uniform prior
  mu = rep(1,nstates)/(nstates)
  
  # Calculating all possible realization of complex attributes (States)
  states = expand.grid(rep(list(rho), nalt - 1 ))
  states = t(as.matrix(states))
  states = rbind(states, 0)
  row.names(states) = c("alt1","outside")
  
  for(i in 1:N){
    #
    # Creating the design matrix
    #
    # simple attributes
    if(complex=="bonus"){   # if price is simple 
      price = 0:3
      X = rbind(cbind(diag(nalt-1),sample(price,nalt-1,T)), rep(0,nalt)) # last row of X is for the outside good with all its elements zero
    }else if(complex=="price"){   # if bonus is simple 
      # discretized complex price
      bonus = simplelevels#c(2.5,3,3.5,4,4.5)#(2,2.2,2.5,2.8,3)#c(-2,-1.2,-1,-0.7,-0.5)##c(-2,-1.2,-1,-0.8,-0.5)#c(-5,-2)#
      X = rbind(cbind(diag(nalt-1),sample(bonus,nalt-1,T)), rep(0,nalt)) 
    }
    
    #
    # Calculating cond. and uncond. choice probs using Blahut Arimoto algorithm
    #
    out = BA(X, beta, lambda, mu, states, N)
    
    
    # adding the realized complex attribute to the design 
    # Draw the state of the world.
    curr_state <- sample(1:length(mu), 1, prob = mu)
    X = cbind(X, states[,curr_state])
    
    
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

