rm(list = ls())
setwd("C:/Users/Matteo/Desktop/RI/R-codes/R_Files")

set.seed(66)
library(ggplot2)
library(gridExtra)
source("RICBC_choice_sets_and_incentives.R")
source("RICBC_RI_SC_Choice.R")


#preferences
beta.vec <- c(1, -1)

#processing cost
lambda <- 0.5

#Discount Scheme
B <- 2#0.75#2
rho <- c(2*(-B),-B,-B,0)
rho.prob <- rep(1, length(rho))/length(rho)

#correlation parameter
theta<-seq(-1,1,by=0.02)

#design matrix
X <- matrix(c(1,2.75,#3.3,
              0,0), byrow = TRUE, ncol=2)

#creating space
condprobinside<-matrix(rep(0,length(theta)*length(rho)), ncol=length(rho))
condprobout<-matrix(rep(0,length(theta)*length(rho)), ncol=length(rho))


for (i in 1:length(theta)){

# Calculation of the states.
states.and.prior <- CreateStatesAndPrior(X, beta.vec, rho, rho.prob)
states.and.prior[[2]]<-c(0.25*(1+theta[i]),0.25*(1-theta[i]),0.25*(1-theta[i]),
                         0.25*(1+theta[i]))

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

Correlation<-rep(theta,2)
Probs1<-as.vector(cbind(condprobinside[,1],condprobout[,1]))
Probs2<-as.vector(cbind(condprobinside[,2],condprobout[,2]))
Probs3<-as.vector(cbind(condprobinside[,3],condprobout[,3]))
Probs4<-as.vector(cbind(condprobinside[,4],condprobout[,4]))
type=c(rep("Inside",length(theta)),rep("Outside",length(theta)))
g<-list("Correlation"=Correlation,"Probs1"=Probs1,"Probs2"=Probs2, "Type"=type)
frame<-data.frame(g)


test1<-ggplot(data=frame, aes(x=Correlation,y=Probs1, color=Type, linetype=Type))+
  geom_line(size=1)+
  labs(title = "Discount d[1] = 2 , d[2] = 2", x="Correlation Parameter"~theta,
       y="Conditional Choice Prob. Discount=2/2")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        legend.title = element_blank(),
        legend.position = "top")+
  scale_linetype_manual(values=c("solid","dashed"))+
  scale_color_manual(values=c( "#F8766D","#00BFC4"))+
  coord_cartesian( ylim = c(0, 1))


test2<-ggplot(data=frame, aes(x=Correlation,y=Probs2, color=Type, linetype=Type))+
  geom_line(size=1)+
  labs(title = "Discount d[1] = 0 , d[2] = 2", x="Correlation Parameter"~theta,
       y="Conditional Choice Prob. Discount= 0/2")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10,face="bold"),
        legend.title = element_blank(),
        legend.position = "top")+
  scale_linetype_manual(values=c("solid","dashed"))+
  scale_color_manual(values=c( "#F8766D","#00BFC4"))+
  coord_cartesian( ylim = c(0, 1))


test3<-ggplot(data=frame, aes(x=Correlation,y=Probs3, color=Type, linetype=Type))+
  geom_line(size=1)+
  labs(title = "Discount d[1] = 2 , d[2] = 0", x="Correlation Parameter"~theta,
       y="Conditional Choice Prob.")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10,face="bold"),
        legend.title = element_blank(),
        legend.position = "top")+
  scale_linetype_manual(values=c("solid","dashed"))+
  scale_color_manual(values=c( "#F8766D","#00BFC4"))+
  coord_cartesian( ylim = c(0, 1))

test4<-ggplot(data=frame, aes(x=Correlation,y=Probs4, color=Type, linetype=Type))+
  geom_line(size=1)+
  labs(title = "Discount d[1] = 0 , d[2] = 0",x="Correlation Parameter"~theta,y="Conditional Choice Prob.")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10,face="bold"),
        legend.title = element_blank(),
        legend.position = "top")+
  scale_linetype_manual(values=c("solid","dashed"))+
  scale_color_manual(values=c( "#F8766D","#00BFC4"))+
  coord_cartesian( ylim = c(0, 1))

Correlation<-theta
Probs1<-condprobinside[,1]
Probs2<-condprobinside[,2]
Probs3<-condprobinside[,3]
Probs4<-condprobinside[,4]
g<-list("Correlation"=Correlation,"Probs1"=Probs1,"Probs2"=Probs2,
        "Probs3"=Probs3,"Probs4"=Probs4, "Type"=type)
frame<-data.frame(g)


test1<-ggplot(data=frame, aes(x=Correlation,y=Probs1))+
  geom_line(size=2,color ="#F8766D")+ 
  ggtitle(bquote(bold(paste("Discounts: ", d[1], " = 2, ", d[2], " = 2"))))+
  labs(x="Correlation Parameter"~theta,y="Cond. Choice Prob.")+
  #labs(x=bquote(bold("Correlation Parameter"~theta)),y="Conditional Choice Prob.")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, size=14))+
  coord_cartesian( ylim = c(0, 1))


test2<-ggplot(data=frame, aes(x=Correlation,y=Probs2))+
  geom_line(size=2,color ="#F8766D")+
  ggtitle(bquote(bold(paste("Discounts: ", d[1], " = 0, ", d[2], " = 2"))))+
  labs(x="Correlation Parameter"~theta,y="Cond. Choice Prob.")+
  #labs(x=bquote(bold("Correlation Parameter"~theta)),y="Conditional Choice Prob")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, size= 14))+
  coord_cartesian( ylim = c(0, 1))


test3<-ggplot(data=frame, aes(x=Correlation,y=Probs3))+
  geom_line(size=2,color ="#F8766D")+
  ggtitle(bquote(bold(paste("Discounts: ", d[1], " = 2, ", d[2], " = 0"))))+
  labs(x="Correlation Parameter"~theta,y="Cond. Choice Prob.")+
  #labs(x=bquote(bold("Correlation Parameter"~theta)),y="Conditional Choice Prob.")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, size=14))+
  coord_cartesian( ylim = c(0, 1))

test4<-ggplot(data=frame, aes(x=Correlation,y=Probs4))+
  geom_line(size=2,color ="#F8766D")+
  ggtitle(bquote(bold(paste("Discounts: ", d[1], " = 0, ", d[2], " = 0"))))+
  labs(x="Correlation Parameter"~theta,y="Cond. Choice Prob.")+
  #labs(x=bquote(bold("Correlation Parameter"~theta)),y="Conditional Choice Prob")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 14))+
  coord_cartesian( ylim = c(0, 1))



windows(height=4.5, width=9)
par(mfrow=c(1,4))
grid.arrange(test4,test2,test3,test1,ncol=2)