rm(list = ls())
setwd("C:/Users/Matteo/Desktop/RI/Plots and Tables")
#setwd("C:/Users/Matteo/Desktop/RI/R-codes/R_Files")

source("RICBC_choice_sets_and_incentives.R")
source("RICBC_RI_SC_Choice.R")

set.seed(66)
library(ggplot2)
library(gridExtra)




beta.vec <- c(6, -1,1)
lambda <- 0.5 #lambda is the parameter that governs the distance between choice probs

p<-seq(4,12, by=0.01)#vector of prices to observe changes in choice probs
choice.set.size<-2#size of choice set

N = length(p)

#discount levels
rho<-c(0,0.5,1,1.5,2,2.5,3,3.5,4)
rho.prob <- rep(1, length(rho))/length(rho)

#creating space for choice probabilities
choice.probs1<-matrix(rep(0,N*(choice.set.size+1)),ncol=choice.set.size+1) 
choice.probs1[,1]=p #collecting prices and choice probs 
statedepprobs1<-array(rep(0,N*choice.set.size*length(rho)), dim=c(N,choice.set.size,length(rho)))


#calculating RI choice probs
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
  choice.probs1[i,2:3]=choice.probs.output$`Choice Probabilities`
  statedepprobs1[i,,]=choice.probs.output$`State Dependent Choice Probabilities`
  
}

# Collecting price levels to achieve fixed choice probability in a given state
# (level of rho) from RI.

#cond choice prob 0.1
index01<-c()
for (i in 1:length(rho)){
  index01[i]<-which(statedepprobs1[,1,i]<0.11 & statedepprobs1[,1,i]>0.09)
}
#index01[9]<-411
#index01[9]<-447

plotps01<-p[index01] 

#cond choice prob 0.5
index05<-c()
for (i in 1:length(rho)){
  index05[i]<-which(statedepprobs1[,1,i]<0.51 & statedepprobs1[,1,i]>0.49)
  #which(x<0.51 & x > 0.49)
}
plotps05<-p[index05]

#cond choice prob 0.9
index09<-c()
for (i in 1:length(rho)){
  index09[i]<-which(statedepprobs1[,1,i]<0.91 & statedepprobs1[,1,i]>0.89)
  #which(round(statedepprobs1[,1,i],2)==0.9)
}
plotps09<-p[index09]


# Calculating price levels to achieve fixed choice probability in a given state
# (level of rho)

f<-function(x, utility, prob){ #function to evaluate price levels X 
  (exp(utility-x))/(exp(utility-x)+1)-prob
}

#cond choice prob 0.1
logit01<-c()
for (i in 1:length(plotps01)){
  logit01[i]<-uniroot(f,c(-100,100), utility=beta.vec[1]+rho[i], prob=0.1)$root
}


#cond choice prob 0.5
logit05<-c()
for (i in 1:length(plotps05)){
  logit05[i]<-uniroot(f,c(-100,100), utility=beta.vec[1]+rho[i], prob=0.5)$root
}

#cond choice prob 0.9

logit09<-c()
for (i in 1:length(plotps09)){
  logit09[i]<-uniroot(f,c(-100,100), utility=beta.vec[1]+rho[i], prob=0.9)$root
}

#we get negative prices here, we set them to zero
#for (i in 1:length(logit09)){
#if (logit09[i]<=0){
#  logit09[i] = NA
#}
#}  



#plotting results
#RI results
g<-list("Price"=c(plotps01,plotps05,plotps09),"Discount"=rho,
        "Type"=c(rep("Choice Prob. = 0.1",length(plotps01)),
                 #rep("Choice Prob. = 0.25",length(plotps01)),
                 rep("Choice Prob. = 0.5",length(plotps01)),
                 rep("Choice Prob. = 0.9", length(plotps01)) )
        #                 rep("Choice Prob. = 0.75",length(plotps01)))
)
df<-data.frame(g)

test1<-ggplot(data=df, aes(x=Discount,y=Price,color=Type,linetype=Type))+
  geom_line( size=1) + 
  geom_point(aes(x = Discount, y = Price, shape=Type), size=2, color="black")+
  ggtitle("RI-DCM")+
  geom_abline(intercept = 5, slope=1)+
  scale_linetype_manual(values=c("dotdash","dashed","dotted"))+#,"dotdash"))+
  scale_color_manual(values=c( "#F8766D","#00BFC4","#00BA38"))+#,"black"))+
  # geom_line(aes(y =diff2, colour="Larger Range"), size=1) + 
  labs(x="Discount",y="Price")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.title = element_blank(),
        legend.position = c(.95, .05),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(5, 5, 5, 5),
        legend.text = element_text(size=12),
        legend.key.size = unit(1.5, "lines"),
        plot.title = element_text(hjust = 0.5, size= 14, face = "bold"))+
  coord_cartesian( ylim = c(4, 10))+
  scale_x_continuous("Discount", labels = as.character(rho), breaks = rho)

#Logit Results
g<-list("Price"=c(logit01,logit05,logit09),"Discount"=rho,
        "Type"=c(rep("Choice Prob. = 0.1",length(plotps01)),
                 #rep("Choice Prob. = 0.25",length(plotps01)),
                 rep("Choice Prob. = 0.5",length(plotps01)),
                 #rep("Choice Prob. = 0.75",length(plotps01)))
                 rep("Choice Prob. = 0.9",length(plotps01)))
)
df<-data.frame(g)

test2<-ggplot(data=df, aes(x=Discount,y=Price,color=Type,linetype=Type))+
  geom_line( size=1) + 
  geom_point(aes(x = Discount, y = Price, shape=Type), size=2, color="black")+
  
  ggtitle("RU logit")+
  geom_abline(intercept = 5, slope=1)+
  scale_linetype_manual(values=c("dotdash","dashed","dotted"))+
  scale_color_manual(values=c( "#F8766D","#00BFC4","#00BA38"))+
  # geom_line(aes(y =diff2, colour="Larger Range"), size=1) + 
  labs(x="Discount",y="Price")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.title = element_blank(),
        legend.position = c(.95, .05),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(5, 5, 5, 5),
        legend.text = element_text(size=12),
        legend.key.size = unit(1.5, "lines"),
        plot.title = element_text(hjust = 0.5, size= 14, face = "bold"))+
  coord_cartesian( ylim = c(4, 10))+
  scale_x_continuous("Discount", labels = as.character(rho), breaks = rho)





windows(height=4.5, width=9)
par(mfrow=c(1,2))
grid.arrange(test1,test2,ncol=2)



