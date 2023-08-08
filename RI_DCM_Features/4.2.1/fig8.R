rm(list = ls())
setwd("C:/Users/Matteo/Desktop/RI/Plots and Tables")
#setwd("C:/Users/Matteo/Desktop/RI/R-codes/R_Files")

source("RICBC_choice_sets_and_incentives.R")
source("RICBC_RI_SC_Choice.R")

set.seed(66)
library(ggplot2)
library(gridExtra)


X<- matrix(c(3.5, 4,#   3.5,  #we create a choice set in which the first alternative has the highest brand utility
             3.25,  4,#3.3,
             3, 4), byrow = TRUE, ncol = 2) #2.9
colnames(X) <- c("brand ", "price")

#Initializing-------------------------------------------------------------------------------
beta.vec <- c(1, -1,1)
lambda <-seq(0.1,4, by=0.05)
rho <- c(0,2)
rho.prob <- rep(1, length(rho))/length(rho)

#calculating utilities-----------------------------------------------------------------------
states.and.prior <- CreateStatesAndPriorNO(X, beta.vec, rho, rho.prob)
#states.and.prior <- CreateStatesAndPrior(X, beta.vec, rho, rho.prob)

#calculating choice probs and mutual information
probs<-matrix(rep(0,3*length(lambda)),ncol=3)
mutual.info<-rep(0,length(lambda))
condprob1<-matrix(rep(0,length(lambda)*ncol(states.and.prior[[1]])),ncol=ncol(states.and.prior[[1]]))
condprob2<-matrix(rep(0,length(lambda)*ncol(states.and.prior[[1]])),ncol=ncol(states.and.prior[[1]]))
condprob3<-matrix(rep(0,length(lambda)*ncol(states.and.prior[[1]])),ncol=ncol(states.and.prior[[1]]))



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
  
  probs[i,]=choice.probs.output$`Choice Probabilities`
  condprob1[i,]=choice.probs.output$`State Dependent Choice Probabilities`[1,]
  condprob2[i,]=choice.probs.output$`State Dependent Choice Probabilities`[2,]
  condprob3[i,]=choice.probs.output$`State Dependent Choice Probabilities`[3,]
}

g<-list("Lambda"=c(lambda,lambda,lambda),
        "Type"=c(rep("Uncond. Choice Prob. a1",length(lambda)),
                 rep("Uncond. Choice Prob. a2",length(lambda)),
                 rep("Uncond. Choice Prob. a3",length(lambda))),
        "Type2"=c(rep("Alternative 1",length(lambda)),
                 rep("Alternative 2",length(lambda)),
                 rep("Alternative 3",length(lambda))),
        "Unconditional_Choice_Probabilities"=as.vector(cbind(probs[,1],probs[,2],probs[,3])),
        "Conditional_Choice_Probabilities_A1"=as.vector(cbind(condprob1[,1],condprob2[,1],condprob3[,1])),
        "Conditional_Choice_Probabilities_A3"=as.vector(cbind(condprob1[,5],condprob2[,5],condprob3[,5]))
        )
#g$Type2=factor(g$Type2, levels=c("CondProba1","CondProba2","CondProba3"), 
#               labels = expression("Alternative"[1], "Alternative"[2],"Alternative"[3]))

df<-data.frame(g)
#test1<-ggplot(data=df, aes(x=Lambda,y=Unconditional_Choice_Probabilities,color=Type,
#                           linetype=Type))+
#  geom_line( size=1)+
#geom_segment( x=2, y =2,aes(xend=2, yend=0), size=0.5, color="black", linetype="solid")+
#  #geom_segment( x=3.3, y =2,aes(xend=3.3, yend=0), size=0.5, color="black", linetype="solid")+
# geom_vline(xintercept=c(2,3.3), linetype="solid", size=0.5)+
#  annotate("text", x=2.1, y=-0.01, label= expression(paste(lambda,"'")), size=5, face="bold") +
#  annotate("text", x=3.4, y=-0.01, label= expression(paste(lambda,"''")), size=5) +
# # geom_text(x=0.9,y=-0.02, label=expression(paste(lambda,"'")), size=5)+
# # geom_text(x=6.6,y=-0.02, label=expression(paste(lambda,"''")), size=5)+
#  scale_linetype_manual(values=c("solid","dashed", "dotted"))+
#  scale_color_manual(values=c( "#F8766D","#00BFC4","#00BA38"))+#

#  
#  labs(x="Information Cost",y="Unconditional Choice Probability")+
#  theme(axis.text=element_text(size=8),
#        axis.title=element_text(size=10,face="bold"),
#        legend.title = element_blank(),
#        legend.position = "top")+
#  coord_cartesian( ylim = c(0, 1))

#windows(height=3.6, width=6.1)
#test1


#plotting conditional choice probs as a function of lambda in a state where
#alternative 1 is the best:
#test2<-ggplot(data=df, aes(x=Lambda,y=Conditional_Choice_Probabilities_A1,color=Type2,
#                           linetype=Type2))+
#  geom_line( size=1)+
#  #geom_segment( x=2, y =2,aes(xend=2, yend=0), size=0.5, color="black", linetype="solid")+
#  #geom_segment( x=3.3, y =2,aes(xend=3.3, yend=0), size=0.5, color="black", linetype="solid")+
#  geom_vline(xintercept=c(2,3.3), linetype="solid", size=0.5)+
#  annotate("text", x=2.1, y=-0.01, label= expression(paste(lambda,"'")), size=5, face="bold") +
#  annotate("text", x=3.4, y=-0.01, label= expression(paste(lambda,"''")), size=5) +
#  # geom_text(x=0.9,y=-0.02, label=expression(paste(lambda,"'")), size=5)+
  # geom_text(x=6.6,y=-0.02, label=expression(paste(lambda,"''")), size=5)+
#  scale_linetype_manual(values=c("solid","dashed", "dotted"))+
#  scale_color_manual(values=c( "#F8766D","#00BFC4","#00BA38"), 
#                     labels=c("test","test2","test3"))+
#  
#  
#  labs(x="Information Cost",y="Conditional Choice Probability")+
#  theme(axis.text=element_text(size=8),
#        axis.title=element_text(size=10,face="bold"),
#        legend.title = element_blank(),
#        legend.position = "top")+
#  coord_cartesian( ylim = c(0, 1))

#windows(height=3.6, width=6.1)
#test2

#plotting conditional choice probs as a function of lambda in a state where
#alternative 3 is the best:

test3<-ggplot(data=df, aes(x=Lambda,y=Conditional_Choice_Probabilities_A3,color=Type2,
                           linetype=Type2))+
  geom_line( size=2)+
  #geom_segment( x=2, y =2,aes(xend=2, yend=0), size=0.5, color="black", linetype="solid")+
  #geom_segment( x=3.3, y =2,aes(xend=3.3, yend=0), size=0.5, color="black", linetype="solid")+
  geom_vline(xintercept=c(0,1.35,3.95), linetype="solid", size=1)+
  annotate("text", x=1.43, y=0.1, label= expression(paste(lambda,"''")), size=10) +
  annotate("text", x=4.05, y=0.1, label= expression(paste(lambda,"'''")), size=10) +
  annotate("text", x=0.08, y=0.1, label= expression(paste(lambda,"'")), size=10) +

  geom_segment( x=0, y =1,aes(xend=0.2, yend=1), size=2, color="#00BA38", linetype="dotted")+
  geom_segment( x=0, y =0,aes(xend=0.2, yend=0), size=1, color="#00BFC4", linetype="dashed")+
  geom_segment( x=0, y =0,aes(xend=0.1, yend=0), size=1, color="#F8766D", linetype="solid")+
  
  # geom_text(x=0.9,y=-0.02, label=expression(paste(lambda,"'")), size=5)+
  # geom_text(x=6.6,y=-0.02, label=expression(paste(lambda,"''")), size=5)+
  scale_linetype_manual(values=c("solid","dashed", "dotted"))+
  scale_color_manual(values=c( "#F8766D","#00BFC4","#00BA38"))+

  
  labs(x="Information Processing Cost",y="Conditional Choice Probability")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.title = element_blank(),
        legend.position = "top",
        legend.key.size = unit(3, "lines"),
        legend.text = element_text(size=14, face = "bold"))+
  

  coord_cartesian( ylim = c(0, 1))

windows(height=3.5, width=6)
test3




