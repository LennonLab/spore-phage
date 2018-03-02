#Section 1: Basic pre-requisites, description of parameters, and parameter values ####

#Clear environment
rm(list=ls())
#Necessary libraries
library(deSolve)

#Explanation of parameters and values used

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
times_ex = seq(from=0,to=10,by=.001)

par(mfrow=c(1,1))
#plot(times_ex,Rmax_default*sin(times_ex)^2,type="l",col="green")
#points(times_ex,Rmax_default*sin(times_ex)^8,type="l",col="red")
#points(times_ex,Rmax_default*sin(times_ex)^32,type="l",col="blue")

#Define a vector of time values over which to simulate
times=seq(from=0,to=60000,by=10)

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


#Section 3: Calculation of parameter effect ####
#Vector of parameter values
parameter_values<-c(e_default,cmax_default,h_default,rad1_default,rad2_default,rda1_default,rda2_default,w_default,K1_default,K2_default,mA_default,mD_default,Rmax_default,n_default)

#The only input, N, is the number of which parameter you want to explore. The indexes (NOT PARAMETER VALUES) correspond to each parameter are below.
#Indexes: e = 1, cmax = 2, h = 3, rad1 = 4, rad2 = 5, rda1 = 6, rda2 = 7, w = 8, K1 = 9, K2 = 10, mA = 11, mD = 12, Rmax = 13, n = 14

parameter_exploration<-function(N){
#Generate the range of parameter values used
parameter_range<-seq(0,parameter_values[N]*2,parameter_values[N]/5)
frequency_output<-array(NA,dim=length(parameter_range))
for(i in 1:length(parameter_range)){
  print(paste("Run",i,"of",length(parameter_range)))
Dorm_output1=lsoda(c(A1=1,D1=0,A2=1,D2=0,R=1),times,Dorm,
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
plot(mD_output[,1],mD_output[,2],xlab="Transition rate, active -> dormant, rad1",ylab="Frequency of the sporulator",type='l')
