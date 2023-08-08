rm(list = ls())

setwd("C:/Users/Matteo/Desktop/RI/Plots and Tables")
source("RICBC_choice_sets_and_incentives.R")
source("RICBC_RI_SC_Choice.R")

set.seed(66)
library(ggplot2)
library(gridExtra)


#choice set
X<- matrix(c(2,    2,  #first column brand, second column price
             1.75,  2,
             1.5, 2,
             1,2), byrow = TRUE, ncol = 2) 
#Note that wlog outside is chosen s.t. it is weakly dominated in all states
colnames(X) <- c("brand ", "price")
nvar<-nrow(X)


#Primitives
beta.vec <- c(1, -1,1) #brand, price, discount
lambda <- 0.5
rho <- c(2, 0)
rho.prob <- rep(1, length(rho))/length(rho)
nstates<-(length(rho))^(nvar-1)


#calculating utilities
#states.and.prior <- CreateStatesAndPrior(X, beta.vec, rho, rho.prob)

#NO outside alternative:

#functions can be found in "RICBC_choice_sets_and_incentives.R"
states.and.prior <- CreateStatesAndPriorNO(X, beta.vec, rho, rho.prob)



#calculating choice probs-------------------------------------------------------------------
choice.probs.output <- CalcChoiceProbsUnderRIWithShannon(
  Omega = states.and.prior[[1]],
  mu = states.and.prior[[2]],
  lambda = lambda,
  max.iter = 10^7,
  precision = 10^(-10)     
)  


#collecting and summarizing choice probs for plotting-----------------------------------------
plot.mat.uncond<-rep(0,nvar) #unconditional choiceprobs
plot.mat.post<-rep(0,nvar) #revealed posteriors
plot.mat.cond<-rep(0,nvar) #conditional choiceprobs

#alternative a4, in the case with four alternatives, yields highest utility
#in state 8. to see type "states.and.prior[[1]]" which displays utilities in all
#possible states (columns).
plot.mat.uncond<-choice.probs.output$`Choice Probabilities`

#choice probs conditional on A4 being the best
plot.mat.cond<-choice.probs.output$`State Dependent Choice Probabilities`[,8]

#plotting results

#c("A1","A2","A3","A4")
xlabels<-expression(a[1],a[2],a[3],a[4])
g<-list("Alternatives"=c("A1","A2","A3","A4"),
        "Conditional_Choice_Probabilities"=plot.mat.cond,
        "Unconditional_Choice_Probabilities"=plot.mat.uncond)
df<-data.frame(g)


test1<-ggplot(data=df, aes(x=Alternatives,y=Conditional_Choice_Probabilities))+
  #ggtitle("Four Alternatives")+
  geom_bar(stat="identity", width=0.5, color="black", fill="#F8766D")+
  labs(x="Alternatives",y="Cond. Choice Prob.")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  coord_cartesian( ylim = c(0, 1))+
  scale_x_discrete(labels= xlabels)



test2<-ggplot(data=df, aes(x=Alternatives,y=Unconditional_Choice_Probabilities))+
  geom_bar(stat="identity", width=0.5, color="black", fill="#00BFC4")+
  labs(x="Alternatives",y="Uncond. Choice Prob.")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        plot.title = element_text(hjust = 0.5, size = 14))+
  coord_cartesian( ylim = c(0, 1))+
  scale_x_discrete(labels= xlabels)
  

#windows(height=3, width=6)
#par(mfrow=c(1,2))
#grid.arrange(test2,test1,ncol=2)


################################################################################
################################################################################
############################# removing alternative 1  ##########################
################################################################################
################################################################################

# smaller choice set
X<-X[-1,] #dropping the best (first) alternative

colnames(X) <- c("brand ", "price")
nvar<-nrow(X)

states.and.prior2 <- CreateStatesAndPriorNO(X, beta.vec, rho, rho.prob)



#calculating choice probs-------------------------------------------------------------------
choice.probs.output <- CalcChoiceProbsUnderRIWithShannon(
  Omega = states.and.prior2[[1]],
  mu = states.and.prior2[[2]],
  lambda = lambda,
  max.iter = 10^7,
  precision = 10^(-10)     
)  



#alternative a4, in the case with only three alternatives, yields highest utility
#in state 4. to see type "states.and.prior2[[1]]" which displays utilities in all
#possible states.
#collecting and summarizing choice probs for plotting-----------------------------------------
plot.mat.uncond2<-choice.probs.output$`Choice Probabilities`
plot.mat.cond2<-choice.probs.output$`State Dependent Choice Probabilities`[,4]

xlabels<-expression(a[2],a[3],a[4])
g<-list("Alternatives"=c("A2","A3","A4"),
        "Unconditional_Choice_Probabilities"=plot.mat.uncond2,
        "Conditional_Choice_Probabilities"=plot.mat.cond2)
df<-data.frame(g)


test3<-ggplot(data=df, aes(x=Alternatives,y=Unconditional_Choice_Probabilities))+
  #ggtitle("Three Alternatives")+
  geom_bar(stat="identity", width=0.5, color="black", fill="#00BFC4")+
  labs(x="Alternatives",y="Uncond. Choice Prob.")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        plot.title = element_text(hjust = 0.5, size = 14, face ="bold"))+
  coord_cartesian( ylim = c(0, 1))+
  scale_x_discrete(labels= xlabels)


test4<-ggplot(data=df, aes(x=Alternatives,y=Conditional_Choice_Probabilities))+
  #ggtitle("Three Alternatives")+
  geom_bar(stat="identity", width=0.5, color="black", fill="#F8766D")+
  labs(x="Alternatives",y="Cond. Choice Prob.")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        plot.title = element_text(hjust = 0.5, size = 14, face ="bold"))+
  coord_cartesian( ylim = c(0, 1))+
  scale_x_discrete(labels= xlabels)


#plotting-----------------------------------------------------------------------
#windows(height=3.5, width=9)
#par(mfrow=c(1,3))
#grid.arrange(test2,test1,test3,ncol=3)

#windows(height=3.5, width=9)
#par(mfrow=c(1,2))
#grid.arrange(test1,test4,ncol=2)

windows(height=3.5, width=9)
par(mfrow=c(2,2))
grid.arrange(test2,test1,test3, test4, ncol=2)
