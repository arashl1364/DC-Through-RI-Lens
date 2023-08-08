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

p<-seq(6,11, by=0.05)#vector of prices to observe changes in choice probs
choice.set.size<-2#size of choice set

N = length(p)


#small range of complex discount
rho <- c(1,2,3,4)
#rho<-c(1,2)
rho.prob <- rep(1, length(rho))/length(rho)

choice.probs1<-matrix(rep(0,N*(choice.set.size+1)),ncol=choice.set.size+1) 
choice.probs1[,1]=p #collecting prices and choice probs 
statedepprobs1<-array(rep(0,N*choice.set.size*length(rho)), dim=c(N,choice.set.size,length(rho)))

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

#next the increased range of the complex discount
rho <- c(0,2,3,5)
#rho<-c(-3,-2,-1,0)
rho.prob <- rep(1, length(rho))/length(rho)



choice.probs2<-matrix(rep(0,N*(choice.set.size+1)),ncol=choice.set.size+1) 
choice.probs2[,1]=p #collecting prices and choice probs 
statedepprobs2<-array(rep(0,N*choice.set.size*length(rho)), dim=c(N,choice.set.size,length(rho)))


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
  choice.probs2[i,2:3]=choice.probs.output$`Choice Probabilities`
  statedepprobs2[i,,]=choice.probs.output$`State Dependent Choice Probabilities`
  
}


#calculating differences between choice probs in two states.
#Note: we calculate differences moving from discount=2 to discount=3 for both cases with 
#smaller and larger range.
diff1<-rep(0,length(p))
diff2<-rep(0,length(p))
for (i in 1:length(p)){
  diff1[i]=abs(statedepprobs1[i,1,2]-statedepprobs1[i,1,3])
  diff2[i]=abs(statedepprobs2[i,1,2]-statedepprobs2[i,1,3])
}

#plot in R
#windows()
#par(mfrow=c(1,2))
#plot(p,diff1,type="l")
#plot(p,diff2,type="l")

#plotting results
g<-list("Price"=c(p,p),"diff"=as.vector(cbind(diff1,diff2)), 
        "probs1"=as.vector(cbind(statedepprobs1[,1,2],statedepprobs2[,1,2])),
        "probs2"=as.vector(cbind(statedepprobs1[,1,3],statedepprobs2[,1,3])),
        "Type"=c(rep("Small Range",length(p)),rep("Large Range",length(p))))
df<-data.frame(g)


test1<-ggplot(data=df, aes(x=Price,y=diff,color=Type,linetype=Type))+
  geom_line( size=1) + 
  scale_linetype_manual(values=c("solid","dashed"))+
  scale_color_manual(values=c( "#F8766D","#00BFC4"))+
 # geom_line(aes(y =diff2, colour="Larger Range"), size=1) + 
  labs(x="Price",y="Choice Probability Difference")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.title = element_blank(),
        legend.position = "top",
        legend.key.size = unit(2, "lines"),
        
        legend.text = element_text(size=12, face = "bold"))+
  coord_cartesian( ylim = c(0, 1))

test2<-ggplot(data=df, aes(x=Price,y=probs1,color=Type,linetype=Type))+
  geom_line( size=1) + 
  scale_linetype_manual(values=c("solid","dashed"))+
  scale_color_manual(values=c( "#F8766D","#00BFC4"))+
  # geom_line(aes(y =diff2, colour="Larger Range"), size=1) + 
  labs(x="Price",y="Choice Probability d = 2")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size=10, face = "bold"))+
  coord_cartesian( ylim = c(0, 1))



#windows(height=4.5, width=12)
#par(mfrow=c(1,2))
#grid.arrange(test1,test2,ncol=2)

windows()
test1


