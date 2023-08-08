rm(list = ls())
#setwd("C:/Users/Matteo/Desktop/RI/R-codes/R_Files")
setwd("~/GitHub/RI-MAMA/RI_DCM_Features")

source("RICBC_choice_sets_and_incentives.R")
source("RICBC_RI_SC_Choice.R")

set.seed(66)
library(ggplot2)
library(gridExtra)
library(latex2exp)

#preferences
beta.vec <- c(1, -1,1)

#processing cost
lambda <- 0.5

#Discount Scheme
rho <- c(2,0)
rho.prob <- rep(1, length(rho))/length(rho)

#correlation parameter
theta<-seq(-1,1,by=0.02)

#design matrix
X <- matrix(c(1,2,
              0,0), byrow = TRUE, ncol=2)

#creating space
condprobinside<-matrix(rep(0,length(theta)*length(rho)), ncol=length(rho))
condprobout<-matrix(rep(0,length(theta)*length(rho)), ncol=length(rho))


for (i in 1:length(theta)){

# Calculation of the states.
states.and.prior <- CreateStatesAndPrior(X, beta.vec, rho, rho.prob)
states.and.prior[[2]]<-c(1/2*(1+theta[i]),1/2*(1-theta[i]))

# RI Optimization with Shannon Costs
choice.probs.output <- CalcChoiceProbsUnderRIWithShannon(
  Omega = states.and.prior[[1]],
  mu = states.and.prior[[2]],
  lambda = lambda, 
  max.iter = 10^7,
  precision = 10^(-10))

condprobinside[i,]<-choice.probs.output$`State Dependent Choice Probabilities`[1,]
condprobout[i,]<-choice.probs.output$`State Dependent Choice Probabilities`[2,]

}

#Correlation<-rep(theta,2)
#Probs1<-as.vector(cbind(condprobinside[,1],condprobout[,1]))
#Probs2<-as.vector(cbind(condprobinside[,2],condprobout[,2]))
#type=c(rep("Inside",length(theta)),rep("Outside",length(theta)))
#g<-list("Correlation"=Correlation,"Probs1"=Probs1,"Probs2"=Probs2, "Type"=type)
#frame<-data.frame(g)


#test1<-ggplot(data=frame, aes(x=Correlation,y=Probs1, color=Type, linetype=Type))+
 # geom_line(size=1)+
  #labs(x="Correlation Parameter",y="Conditional Choice Prob. Discount=-2")+
  #theme(axis.text=element_text(size=8),
   #     axis.title=element_text(size=10,face="bold"),
    #    legend.title = element_blank(),
     #   legend.position = "top")+
  #scale_linetype_manual(values=c("solid","dashed"))+
  #scale_color_manual(values=c( "#F8766D","#00BFC4"))+
  #coord_cartesian( ylim = c(0, 1))


#test2<-ggplot(data=frame, aes(x=Correlation,y=Probs2, color=Type, linetype=Type))+
 # geom_line(size=1)+
  #labs(x="Correlation Parameter",y="Conditional Choice Prob. Discount= 0")+
  #theme(axis.text=element_text(size=8),
   #     axis.title=element_text(size=10,face="bold"),
    #    legend.title = element_blank(),
     #   legend.position = "top")+
  #scale_linetype_manual(values=c("solid","dashed"))+
  #scale_color_manual(values=c( "#F8766D","#00BFC4"))+
  #coord_cartesian( ylim = c(0, 1))


Correlation<-rep(theta,2)
Probs1<-condprobinside[,1]
Probs2<-condprobinside[,2]
g<-list("Correlation"=Correlation,"Probs1"=Probs1,"Probs2"=Probs2)
frame<-data.frame(g)


test1<-ggplot(data=frame, aes(x=Correlation,y=Probs1))+
  geom_line(size=1, color ="#F8766D")+
  ggtitle("Discount d = 2")+
  geom_segment( x=0.78, y =0,aes(xend=0.78, yend=1.5), size=0.6, color="black", linetype="solid")+
  geom_segment( x=-0.78, y =0,aes(xend=-0.78, yend=1.5), size=0.6, color="black", linetype="solid")+
  annotate("text", x=0.85, y=0.1, label= expression(paste(theta,"''")), size=6) +
  annotate("text", x=-0.65, y=0.1, label= expression(paste(theta,"'")), size=6) +
  labs( x="Correlation Parameter"~theta,y="Conditional Choice Prob.")+
  #labs( x=bquote(bold("Correlation Parameter"~theta)),y="Conditional Choice Prob.")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, size=14, face = "bold"))+
  coord_cartesian( ylim = c(0, 1))


test2<-ggplot(data=frame, aes(x=Correlation,y=Probs2))+
  geom_line(size=1, color ="#F8766D")+
  geom_segment( x=0.78, y =0,aes(xend=0.78, yend=1.5), size=0.6, color="black", linetype="solid")+
  geom_segment( x=-0.78, y =0,aes(xend=-0.78, yend=1.5), size=0.6, color="black", linetype="solid")+
  annotate("text", x=0.85, y=0.1, label= expression(paste(theta,"''")), size=6) +
  annotate("text", x=-0.65, y=0.1, label= expression(paste(theta,"'")), size=6) +
  ggtitle("Discount d = 0")+
  labs( x="Correlation Parameter"~theta,y="Conditional Choice Prob.")+
  #labs(title = "Discount d=0" , x=bquote(bold("Correlation Parameter"~theta)),
  #    y="Conditional Choice Prob.")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, size=14, face = "bold"))+
  coord_cartesian( ylim = c(0, 1))

windows(height=4, width=10)
par(mfrow=c(1,2))
grid.arrange(test2,test1,ncol=2)