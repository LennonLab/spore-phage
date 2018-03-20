#Section 1: Basic pre-requisites, description of parameters, and parameter values ####

#Clear environment
rm(list=ls())
#Necessary libraries
library(deSolve)

#Explanation of parameters and values used

# extinction threshold
extinct = 1

#Efficiency of converting resources to bacteria
e_default = 1

#Maximum consumption rate of resources by bacteria
cmax_default = 1

#Half-saturation constant of resource uptake by bacteria
h_default = 5

#Maximum rate at which active bacteria of strain i transition to dormancy.
#Strain w does not exhibit dormancy under these parameters.
rad1_default = 0.5
rad2_default = 0

#Maximum rate at which dormant bacteria of strain i resuscitate
rda1_default = 0.1
rda2_default = 0

#Steepness of transition rate's dependence on resources
w_default = 5

#Half-saturation constant for transition from active to dormant 
#or vice versa for genotype i
K1_default = 5
K2_default = 5 

#Example plot showing the effect of this m parameter
#Example vector of resources
R1 = seq(from=0,to=10,by=0.1)

#plot(R1,1-rad1_default*R1^2/(R1^2+K1_default^2),type="l",ylim=c(.5,1),xlab="Resources",ylab="Rate of Transition from Active to Dormant")
#points(R1,1-rad1_default*R1^4/(R1^4+K1_default^4),type="l",col="blue")
#points(R1,1-rad1_default*R1^8/(R1^8+K1_default^8),type="l",col="red")
#points(R1,1-rad1_default*R1^16/(R1^16+K1_default^16),type="l",col="green")
#legend(7,.9,lty=1,col=c("black","blue","red","green"),legend=c("w=2","w=4","w=6","w=8","w=16"))

#plot(R1,rad1_default*R1^2/(R1^2+K1_default^2),type="l",ylim=c(0,.5),xlab="Resources",ylab="Rate of Transition from Dormant to Active")
#points(R1,rad1_default*R1^4/(R1^4+K1_default^4),type="l",col="blue")
#points(R1,rad1_default*R1^8/(R1^8+K1_default^8),type="l",col="red")
#points(R1,rad1_default*R1^16/(R1^16+K1_default^16),type="l",col="green")
#legend(7,.25,lty=1,col=c("black","blue","red","green"),legend=c("w=2","w=4","w=6","w=8","w=16"))

#Background mortality rate of active bacteria or dormant bacteria
mA_default = 0.05
mD_default = 5.1565e-4

#Maximum rate of inflow of resources
Rmax_default = 60


#This parameter determines the sharpness of resources influx peaks. This parameter
#should be even so that sin(t/200)^n is always non-negative and thus resource influx
#is never negative.
n_default = 100

#Plot to show the effect of the parameter n
#Example time series
times_ex = seq(from=0,to=2000,by=10)

par(mfrow=c(1,1))
plot(times_ex,Rmax_default*sin(times_ex/200)^150,type="l",col="black", ylab = "Resource input", xlab="time")
points(times_ex,Rmax_default*sin(times_ex/200)^100,type="l",col="blue")
points(times_ex,Rmax_default*sin(times_ex/200)^50,type="l",col="red")
points(times_ex,Rmax_default*sin(times_ex/200)^2,type="l",col="green")
legend("topleft",lty=1,col=c("black","blue","red","green"),legend=c("n=150","n=100","n=50","n=2"))

#Define a vector of time values over which to simulate
# times=seq(from=0,to=60000,by=10)
times=seq(from=0,to=1000,by=0.1)

#Section 2: System of differential equations ####

#Define dormancy and competition model. The differential equation is a function
#of time, a vector of state variables (y), and parameters
Dorm = function(t,y,params){
  #Tell R what the 5 state variables are.
  A1 = y[1]; D1 = y[2]; A2 = y[3]; D2 = y[4]; R = y[5]
  with(
    as.list(params),
    {

      #Define the equation governing the rate of change of each state variable
      # dA1 = e*cmax*A1*R/(R+h)-rad1*(1-R^w/(R^w+K1^w))*A1+rda1*R^w/(R^w+K1^w)*D1-mA*A1
      # dD1 = rad1*(1-R^w/(R^w+K1^w))*A1-rda1*R^w/(R^w+K1^w)*D1-mD*D1
      # dA2 = e*cmax*A2*R/(R+h)-rad2*(1-R^w/(R^w+K2^w))*A2+rda2*R^w/(R^w+K2^w)*D2-mA*A2
      # dD2 = rad2*(1-R^w/(R^w+K2^w))*A2-rda2*R^w/(R^w+K2^w)*D2-mD*D2
      dA1 = e*cmax*A1*R/(R+h)-rad1*(1-R^w/(R^w+K1^w))*A1+rda1*R^w/(R^w+K1^w)*D1-mA*A1
      dD1 = rad1*(1-R^w/(R^w+K1^w))*A1-rda1*R^w/(R^w+K1^w)*D1-mD*D1
      dA2 = e*cmax*A2*R/(R+h)-rad2*(1-R^w/(R^w+K2^w))*A2+rda2*R^w/(R^w+K2^w)*D2-mA*A2
      dD2 = rad2*(1-R^w/(R^w+K2^w))*A2-rda2*R^w/(R^w+K2^w)*D2-mD*D2
      #Resource influx oscillates through time. A peak in resource influx occurs 
      #every few hundred time units
      dR = Rmax*sin(t/200)^n-cmax*A1*R/(R+h)-cmax*A2*R/(R+h)
      res = c(dA1,dD1,dA2,dD2,dR)
      list(res)
    }
  )
}
# 
# rootfunc <- function(t, y, params){
#   return(y[3] - 1)
# }
# 
# eventfunc <- function(t, y, params){
#   y[3] = 0
#   return(y)
# }

eventfunc <- function(t, y, parms){
  # with(as.list(y), {
    if(y[1] < 1) y[1] <- 0
    if(y[2] < 1) y[2] <- 0
    if(y[3] < 1) y[3] <- 0
    if(y[4] < 1) y[4] <- 0
    return(y)
  # })
}


#Section 3: Calculation of parameter effect ####
#Vector of parameter values
parameter_values<-c(e_default,cmax_default,h_default,rad1_default,rad2_default,rda1_default,rda2_default,w_default,K1_default,K2_default,mA_default,mD_default,Rmax_default,n_default)

Dorm_output1=lsoda(c(A1=10,D1=0,A2=10,D2=0,R=Rmax_default),times,Dorm,
                   parm = c(e=e_default,cmax=cmax_default,h=h_default,
                             rad1=rad1_default,rad2=rad2_default,rda1=rda1_default,rda2=rda2_default,
                             w=w_default,K1=K1_default,K2=K2_default,mA=mA_default,mD=mD_default,
                             Rmax=Rmax_default,n=4),
                   # rootfun = rootfunc, 
                   events = list(func = eventfunc, root = F, time = times))
plot(Dorm_output1)

tail(Dorm_output1)
#The only input, N, is the number of which parameter you want to explore. The indexes (NOT PARAMETER VALUES) correspond to each parameter are below.
#Indexes: e = 1, cmax = 2, h = 3, rad1 = 4, rad2 = 5, rda1 = 6, rda2 = 7, w = 8, K1 = 9, K2 = 10, mA = 11, mD = 12, Rmax = 13, n = 14
################################
parameter_exploration<-function(N){
#Generate the range of parameter values used
parameter_range<-seq(0,parameter_values[N]*2,parameter_values[N]/5)
frequency_output<-array(NA,dim=length(parameter_range))
for(i in 1:length(parameter_range)){
  print(paste("Run",i,"of",length(parameter_range)))
Dorm_output1=lsoda(c(A1=1,D1=0,A2=1,D2=0,R=Rmax_default),times,Dorm,
                c(e=ifelse(N==1,parameter_range[i],e_default),cmax=ifelse(N==2,parameter_range[i],cmax_default),h=ifelse(N==3,parameter_range[i],h_default),rad1=ifelse(N==4,parameter_range[i],rad1_default),rad2=ifelse(N==5,parameter_range[i],rad2_default),rda1=ifelse(N==6,parameter_range[i],rda1_default),rda2=ifelse(N==7,parameter_range[i],rda2_default),w=ifelse(N==8,parameter_range[i],w_default),K1=ifelse(N==9,parameter_range[i],K1_default),K2=ifelse(N==10,parameter_range[i],K2_default),mA=ifelse(N==11,parameter_range[i],mA_default),mD=ifelse(N==12,parameter_range[i],mD_default),Rmax=ifelse(N==13,parameter_range[i],Rmax_default),n=ifelse(N==14,parameter_range[i],n_default)))
frequency_output[i]=mean(Dorm_output1[round(.7*length(times),digits=0):length(times),2])/(mean(Dorm_output1[round(.7*length(times),digits=0):length(times),2])+mean(Dorm_output1[round(.7*length(times),digits=0):length(times),4]))
}
return(cbind(parameter_range,frequency_output))
}

#Section 4: Input required for final output ####

#The only input, N, is the number of which parameter you want to explore. The indexes (NOT PARAMETER VALUES) correspond to each parameter are below.
#Indexes: e = 1, cmax = 2, h = 3, rad1 = 4, 
#         rad2 = 5, rda1 = 6, rda2 = 7, w = 8, 
#         K1 = 9, K2 = 10, mA = 11, mD = 12, Rmax = 13, n = 14

#This input is passed to the parameter_exploration function which returns a table with two columns.
#The first column is the range of parameter values used. The second column holds the corresponding
#frequency of the sporulator clone in the active population that results from those parameter values.

#Example of exploring the effect of mD, which is parameter 12.
mD_output<-parameter_exploration(14)
par(mfrow=c(1,1),mar=c(4,4,2,2))
plot(mD_output[,1],mD_output[,2],xlab="Resource volatility (n)",ylab="Frequency of genotype",type='l',ylim = c(0,1))
points(mD_output[,1],mD_output[,2])
points(mD_output[,1],1-mD_output[,2], col="red", type='l')
points(mD_output[,1],1-mD_output[,2], col="red")
abline(h=0.5, lty=2)
legend("topright",lty=1,col=c("black","red"),legend=c("sporulator", "nonsporulator"), text.col=c("black","red"))



mD_output<-parameter_exploration(12)
par(mfrow=c(1,1),mar=c(4,4,2,2))
plot(mD_output[,1],mD_output[,2],ylim = c(0,1),
     xlab="Dormant mortality (mD)",ylab="Frequency of the sporulator",type='l')
points(mD_output[,1],mD_output[,2])
abline(h=0.5, lty=2)

####################################
# explore all paramaeters
param.names <- c("e","cmax","h","rad1","rad2","rda1","rda2","w","K1","K2","mA","mD","Rmax","n")
for(i in 12:length(param.names)){
  mD_output<-parameter_exploration(i)
  par(mfrow=c(1,1),mar=c(4,4,2,2))
  pdf(file = paste0("../figures/param_explorer/",param.names[i],".pdf"))
  plot(mD_output[,1],mD_output[,2],ylim = c(0,1),
       xlab=param.names[i],ylab="Frequency of the sporulator",type='l')
  points(mD_output[,1],mD_output[,2])
  points(mD_output[which(mD_output[,1]==parameter_values[i]),1],mD_output[which(mD_output[,1]==parameter_values[i]),2], pch=19)
  abline(h=0.5, lty=2)
  dev.off()
}

######
####################################
# dynamics with default values
Dorm_output1=euler(c(A1=10,D1=0,A2=10,D2=0,R=Rmax_default),times,Dorm,
                   parms = c(e=e_default,cmax=cmax_default,h=h_default,
                             rad1=rad1_default,rad2=rad2_default,rda1=rda1_default,rda2=rda2_default,
                             w=w_default,K1=K1_default,K2=K2_default,mA=mA_default,mD=mD_default,
                             Rmax=Rmax_default,n=6))
plot(Dorm_output1)
# par(mfrow=c(1,1))
# hist(Dorm_output1[,'A1'])
A1A2.sum <- apply(Dorm_output1[,c('A1','A2')], 1, sum)
freq.A1 <- Dorm_output1[,'A1']/A1A2.sum 
par(mfrow=c(2,3))
for(i in 2:6){
  plot(Dorm_output1[,1], Dorm_output1[,i], type="l", main=colnames(Dorm_output1)[i], xlab = "time", ylab=colnames(Dorm_output1)[i])
}
plot(Dorm_output1[,1],freq.A1 , type="l", main='A1/A1+A2', xlab = "time", ylab='A1/A1+A2', ylim=c(0,1))

#plot the last tx time steps
tx <- c(floor(nrow( Dorm_output1)*0.95):nrow( Dorm_output1))
par(mfrow=c(2,3))
for(i in c(2,3,4,6)){
  plot(Dorm_output1[tx,1], Dorm_output1[tx,i], type="l", main=colnames(Dorm_output1)[i], xlab = "time", ylab=colnames(Dorm_output1)[i])
}
plot(Dorm_output1[tx,1], freq.A1[tx] , type="l", main='A1/(A1+A2)', xlab = "time", ylab='A1/(A1+A2)', ylim=c(0,1))
plot(Dorm_output1[tx,1], Dorm_output1[tx,4] , type="l", main='A1 (black) Vs. A2(red)', xlab = "time", ylab='A1 or A2', col="red")
points(Dorm_output1[tx,1], Dorm_output1[tx,2] , type="l")
     
# zoom in on last steps and compare cells and spores
par(mfrow=c(1,1))
plot(Dorm_output1[tx,1], Dorm_output1[tx,3] , type="l", main='A1 (black); A2 (red); D1 (blue)', xlab = "time", ylab='cell num. (log)', col="blue", log='y', ylim=c(1e-20,1e4))
points(Dorm_output1[tx,1], Dorm_output1[tx,2] , type="l")
points(Dorm_output1[tx,1], Dorm_output1[tx,4] , type="l", col="red")


#####################
# ploting the dunctional forms used
w <- w_default
k <- K1_default
R <- R1
t <- seq(from=0,to=6000,by=10)

# plotting the dormancy and resucitation dependence on resources
par(mfrow=c(1,2))
plot(R, 1-R^w/(R^w+k^w), main = "Dormancy dependence on Resources(R)")
plot(R,R^w/(R^w+k^w), main = "Resucitation dependence on Resources(R)")

# ploting resource input in different n values
par(mfrow=c(2,3))
plot(t,cos(t/200)^3, main = "n=2",type='l')
plot(t,sin(t/200)^4, main = "n=4",type='l')
plot(t,sin(t/200)^6, main = "n=6",type='l')
plot(t,sin(t/200)^8, main = "n=8",type='l')
plot(t,sin(t/200)^10, main = "n=10",type='l')
plot(t,sin(t/200)^20, main = "n=20",type='l')

