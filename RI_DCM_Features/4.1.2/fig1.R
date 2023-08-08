rm(list = ls())
setwd("C:/Users/Matteo/Desktop/RI/Plots and Tables")
#setwd("C:/Users/Matteo/Desktop/RI/R-codes/R_Files")

source("RICBC_choice_sets_and_incentives.R")
source("RICBC_RI_SC_Choice.R")

set.seed(66)
library(ggplot2)
library(gridExtra)


beta.vec <- c(2, -1,1) #preference vector, last entry is complex
lambda <- 0.5 #information cost
rho <- c(0,2) #complex attribute levels
#rho <- seq(0,-2, by=-0.25)
rho.prob <- rep(1, length(rho))/length(rho)
p<-seq(2,4, by=0.02)#simple price levels

#creating necessary objects
N<-length(p)
choice.probs.diff<-rep(0,N)

#calculate choice probs for different price levels 
for(i in 1:N){
  
  X <- matrix(c(1,0,p[i],0), ncol=length(beta.vec)-1) #fixed design matrix
  
  states.and.prior <- CreateStatesAndPrior(X, beta.vec, rho, rho.prob)
  
  choice.probs.output <- CalcChoiceProbsUnderRIWithShannon(
    Omega = states.and.prior[[1]],
    mu = states.and.prior[[2]],
    lambda = lambda,
    max.iter = 10^7,
    precision = 10^(-10)     
  )
  choice.probs.diff[i]=abs(choice.probs.output$`State Dependent Choice Probabilities`[1,1]-
    choice.probs.output$`State Dependent Choice Probabilities`[1,2])
}

#plot probability difference
g<-list("Price"=p,"Diff"=choice.probs.diff)
frame<-data.frame(g)


test1<-ggplot(data=frame, aes(x=Price,y=Diff))+
  geom_line(size=2,color ="#F8766D")+
  ggtitle("Choice Prob. Difference")+
  labs(x="Price",y="Choice Probability Difference")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.title = element_blank(),
        legend.position = "top",
        legend.key.size = unit(1, "cm"),
        plot.title = element_text(hjust = 0.5, size=14, face = "bold"))+
  #scale_linetype_manual(values=c("dashed","solid"))+
  #scale_color_manual(values=c( "#00BFC4","#F8766D"))+
  coord_cartesian( ylim = c(0, 1))

#windows(height=3, width=5.5)
#test1

#plot choice probs for both states d=0 and d=1
g<-list("Price"=p,"Diff"=choice.probs.diff)
frame<-data.frame(g)


choice.probs=matrix(rep(0,2*length(p)),ncol=2)
for(i in 1:N){
  
  X <- matrix(c(1,0,p[i],0), ncol=length(beta.vec)-1) #fixed design matrix
  
  states.and.prior <- CreateStatesAndPrior(X, beta.vec, rho, rho.prob)
  
  choice.probs.output <- CalcChoiceProbsUnderRIWithShannon(
    Omega = states.and.prior[[1]],
    mu = states.and.prior[[2]],
    lambda = lambda,
    max.iter = 10^7,
    precision = 10^(-10)     
  )
  choice.probs[i,1]=choice.probs.output$`State Dependent Choice Probabilities`[1,1]
  choice.probs[i,2]=choice.probs.output$`State Dependent Choice Probabilities`[1,2]
}

g<-list("Price"=c(p,p),"probs"=c(choice.probs[,1],c(choice.probs[,2])), 
                                 "Type"=c(rep("d = 0",length(p)),rep("d = 2",length(p)))
                                 )
frame<-data.frame(g)

test2<-ggplot(data=frame, aes(x=Price,y=probs, linetype=Type, color=Type))+
  geom_line(size=2)+
  ggtitle("Choice Prob. of Inside Alternative")+
  labs(x="Price",y="Choice Probability")+
  scale_linetype_manual(values=c("dotted","solid"))+
  scale_color_manual(values=c( "#00BFC4","#F8766D"))+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        legend.key.size = unit(2, "lines"),
        
        #legend.key.width,
       # legend.position = "top",
        plot.title = element_text(hjust = 0.5, size=14, face = "bold"))+
  #scale_linetype_manual(values=c("dashed","solid"))+
  #scale_color_manual(values=c( "#00BFC4","#F8766D"))+
  coord_cartesian( ylim = c(0, 1))

#windows(height=3, width=5.5)
#test2

#arrange both plots
windows(height=4.5, width=12)
par(mfrow=c(1,2))
grid.arrange(test2,test1,ncol=2)