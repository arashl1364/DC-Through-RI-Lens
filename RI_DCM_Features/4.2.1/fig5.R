

rm(list = ls())
setwd("C:/Users/Matteo/Desktop/RI/Plots and Tables")
#setwd("C:/Users/Matteo/Desktop/RI/R-codes/R_Files")

source("RICBC_choice_sets_and_incentives.R")
source("RICBC_RI_SC_Choice.R")

set.seed(66)
library(ggplot2)
library(gridExtra)

#In this file we compute choice probabilities, mutual information, and 
# the impact of a change in the complex attribute for an inside good that is the
#only inside option in a simple choice set. Simple attributes are price and brand,
#complex attribute is discount.

X<- matrix(c(1,    2,  
             0,0 ), byrow = TRUE, ncol = 2)
colnames(X) <- c("brand ", "price")

#Initializing-------------------------------------------------------------------------------
beta.vec <- c(1.2, -1,1) #brand, price, discount preferences
lambda <-seq(0.01,4.01, by=0.02) #sequence of information processing costs

rho <- c(0,2) #complex attribute levels
rho.prob <- rep(1, length(rho))/length(rho)

#calculating utilities-----------------------------------------------------------------------
states.and.prior <- CreateStatesAndPrior(X, beta.vec, rho, rho.prob)

#creating space for objects
probs<-matrix(rep(0,2*length(lambda)),ncol=2)
mutual.info<-rep(0,length(lambda))
stateprobs<-array(rep(0,length(lambda)*2*2),dim=c(length(lambda),2,2))
impact.discount<-rep(0,length(lambda))


for (i in 1:length(lambda)){
  choice.probs.output <- CalcChoiceProbsUnderRIWithShannon(
    Omega = states.and.prior[[1]],
    mu = states.and.prior[[2]],
    lambda = lambda[i],
    max.iter = 10^7,
    precision = 10^(-10)     
  )  
  
  mutual.info[i] <- MutualInfo(prior = states.and.prior[[2]], 
                               signal.dist = choice.probs.output[[1]], 
                               posterior = choice.probs.output[[3]])
  
  stateprobs[i,,]=choice.probs.output$`State Dependent Choice Probabilities`
  probs[i,]=choice.probs.output$`Choice Probabilities`
  impact.discount[i]=abs(stateprobs[i,1,1]-stateprobs[i,1,2])
}


#collect choice probs of inside good
#probabilities for bad complex attribute (=0)
state.probs<-rep(0,length(lambda))
for (i in 1:length(lambda)){
  state.probs[i]=stateprobs[i,1,2]
}

#probabilities for good complex attribute (=2)
state.probs2<-rep(0,length(lambda))
for (i in 1:length(lambda)){
  state.probs2[i]=stateprobs[i,1,1]
}

g<-list("Lambda"=lambda,
        "State_Probs"=state.probs,
        "State_Probs2"=state.probs2,
        "Mutual_Info"=mutual.info,
        "Impact_Discount"=impact.discount)
df<-data.frame(g)

test1<-ggplot(data=df, aes(x=Lambda))+
  ggtitle(bquote(bold(paste("Cond. Choice Prob. d = 2"))))+
  geom_line(aes(y =state.probs), size=1, color="#F8766D") + 
  labs(x="Information Processing Cost",y="Cond. Choice Prob.")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        plot.title = element_text(hjust = 0.5, size = 14))+
  geom_segment( x=0, y =0,aes(xend=0.01, yend=0), size=1,color="#F8766D")+ #little trick until simulator is rewritten to include lambda=0
  geom_segment( x=2.42, y =0,aes(xend=2.42, yend=1.5), size=0.6, color="black", linetype="solid")+
  annotate("text", x=2.55, y=0.1, label= expression(paste(lambda,"''")), size=6) +
  geom_segment( x=0, y =0,aes(xend=0, yend=1.5), size=0.6, color="black", linetype="solid")+
  annotate("text", x=0.08, y=0.1, label= expression(paste(lambda,"'")), size=6) +
  
  coord_cartesian( ylim = c(0, 1))

test12<-ggplot(data=df, aes(x=Lambda))+
  ggtitle(bquote(bold(paste("Cond. Choice Prob. d = 0"))))+
  geom_line(aes(y =state.probs2), size=1, color="#F8766D") + 
  labs(x="Information Processing Cost",y="Cond. Choice Prob.")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        plot.title = element_text(hjust = 0.5, size = 14))+
  #geom_segment( x=0, y =1,aes(xend=0.01, yend=1), size=1,color="#F8766D")+ #little trick until simulator is rewritten to include lambda=0
  geom_segment( x=2.42, y =0,aes(xend=2.42, yend=1.5), size=0.6, color="black", linetype="solid")+
  annotate("text", x=2.55, y=0.1, label= expression(paste(lambda,"''")), size=6) +
  geom_segment( x=0, y =0,aes(xend=0, yend=1.5), size=0.6, color="black", linetype="solid")+
  annotate("text", x=0.08, y=0.1, label= expression(paste(lambda,"'")), size=6) +
  
  coord_cartesian( ylim = c(0, 1))

test2<-ggplot(data=df, aes(x=Lambda))+
  ggtitle(bquote(bold(paste("Mutual Information"))))+
  geom_line(aes(y =mutual.info), size=1, color="#F8766D") + 
  labs(x="Information Processing Cost",y="Mutual Information")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        plot.title = element_text(hjust = 0.5, size = 14))+
 
  geom_segment( x=2.42, y =0,aes(xend=2.42, yend=1.5), size=0.6, color="black", linetype="solid")+
  annotate("text", x=2.55, y=0.1, label= expression(paste(lambda,"''")), size=6) +
  geom_segment( x=0, y =0,aes(xend=0, yend=1.5), size=0.6, color="black", linetype="solid")+
  annotate("text", x=0.08, y=0.1, label= expression(paste(lambda,"'")), size=6) +
  
   coord_cartesian( ylim = c(0, 0.8))

test3<-ggplot(data=df, aes(x=Lambda))+
  ggtitle(bquote(bold(paste("Choice Prob. Difference"))))+
  geom_line(aes(y =impact.discount), size=1, color="#F8766D") + 
  labs(x="Information Processing Cost",y="Choice Prob. Difference")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        plot.title = element_text(hjust = 0.5, size = 14))+
  geom_segment( x=2.42, y =0,aes(xend=2.42, yend=1.5), size=0.6, color="black", linetype="solid")+
  annotate("text", x=2.55, y=0.1, label= expression(paste(lambda,"''")), size=6) +
  geom_segment( x=0, y =0,aes(xend=0, yend=1.5), size=0.6, color="black", linetype="solid")+
  annotate("text", x=0.08, y=0.1, label= expression(paste(lambda,"'")), size=6) +
  
   coord_cartesian( ylim = c(0, 1))

#windows(height=3.3, width=6)
#par(mfrow=c(1,2))
#grid.arrange(test2,test3,ncol=2)

#windows(height=3.3, width=6)
#par(mfrow=c(1,2))
#grid.arrange(test1,test12,ncol=2)

windows(height=4, width=6)
par(mfrow=c(2,2))
grid.arrange(test2,test3,test12,test1,ncol=2)