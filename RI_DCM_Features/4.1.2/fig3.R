rm(list = ls())
setwd("C:/Users/Matteo/Desktop/RI/Plots and Tables")

source("RICBC_choice_sets_and_incentives.R")
source("RICBC_RI_SC_Choice.R")

set.seed(66)
library(ggplot2)
library(gridExtra)


beta.vec <- c(6, -1,1)
lambda <- 0.5 #lambda information processing cost

p<-seq(5,9, by=0.5)#vector of prices to observe changes in choice probs
choice.set.size<-2#size of choice set

N = length(p)

#discount levels
rho<-c(0,0.5,1,1.5,2,2.5,3,3.5,4)
#distribution of discount levels
rho.prob <- rep(1, length(rho))/length(rho)


#creating space for choice probabilities
statedepprobs<-matrix(rep(0,N*length(rho)),ncol=length(rho)) 


#calculating RI choice probs
for(i in 1:N){
  
  #functions can be found in RICBC_choice_sets_and_incentives.R
  X <- matrix(c(1,0,p[i],0), ncol=length(beta.vec)-1) #fixed design matrix
  
  states.and.prior <- CreateStatesAndPrior(X, beta.vec, rho, rho.prob)
  
  choice.probs.output <- CalcChoiceProbsUnderRIWithShannon(
    Omega = states.and.prior[[1]],
    mu = states.and.prior[[2]],
    lambda = lambda,
    max.iter = 10^7,
    precision = 10^(-10)     
  )
  statedepprobs[i,]=choice.probs.output$`State Dependent Choice Probabilities`[1,]
  
}

RIprob=diag(statedepprobs)


#Logit choice probs
util<-1 #fixed utility level of inside good

#function to evaluate logit choice probs with one inside and one outside
f<-function(x){
  (exp(util))/(exp(util)+1)-x
}
logitprob<-uniroot(f,c(-100,100))$root
logitprob<-rep(logitprob,length(rho))

#plotting results
#RI results
g<-list("Probs"=c(RIprob, logitprob),"Discount"=rho,
        "Type"=c(rep("RI-DCM",length(rho)),
                 rep("RU logit",length(rho))
        )
)
df<-data.frame(g)


#plotting
test1<-ggplot(data=df, aes(x=Discount,y=Probs,color=Type,linetype=Type))+
  geom_line( size=1) + 
  ggtitle("Price")+
  geom_point(aes(x=Discount, y=Probs, shape=Type), colour="black", size=3)+
  #ggtitle("Rational Inattention")+
  scale_linetype_manual(values=c("solid","dashed","dotted"))+
  scale_color_manual(values=c( "#F8766D","#00BFC4","#00BA38"))+
  scale_x_continuous(sec.axis= ~.+5)+#, name = "Price")+
  # geom_line(aes(y =diff2, colour="Larger Range"), size=1) + 
  labs(x="Discount",y="Conditional Choice Prob.")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.title = element_blank(),
        legend.position = c(.05, .05),
        legend.justification = c("left", "bottom"),
        legend.box.just = "left",
        legend.margin = margin(7, 7, 7, 7),
        legend.text = element_text(size=14),
        legend.key.size = unit(4, "lines"),
        plot.title = element_text(hjust = 0.5, size= 14))+
  coord_cartesian( ylim = c(0, 1))


windows()
test1