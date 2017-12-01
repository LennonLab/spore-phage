#Clear environment
rm(list=ls())
#Necessary libraries
library(deSolve)
###
#Define dormancy and competition model. The differential equation is a function
#of time, a vector of state variables (y), and parameters
Dorm = function(t,y,params){
  #Tell R what the 7 state variables are.
  A1 = y[1]; D1 = y[2]; A2 = y[3]; D2 = y[4]; I = y[5]; V = y[6]; R = y[7];
  with(
    as.list(params),
{
  #Define the equation governing the rate of change of each state variable
  dA1 = e*cmax*A1*R/(R+h)-rad1*(1-R^w/(R^w+K1^w))*A1+rda1*R^w/(R^w+K1^w)*D1-mA*A1 - adsorption*A1*V
  dD1 = rad1*(1-R^w/(R^w+K1^w))*A1-rda1*R^w/(R^w+K1^w)*D1-mD*D1
  dA2 = e*cmax*A2*R/(R+h)-rad2*(1-R^w/(R^w+K2^w))*A2+rda2*R^w/(R^w+K2^w)*D2-mA*A2 - adsorption*A2*V
  dD2 = rad2*(1-R^w/(R^w+K2^w))*A2-rda2*R^w/(R^w+K2^w)*D2-mD*D2
  # infected cells, assuming infection is equal on A1 and A2
  dI  = adsorption*A1*V + adsorption*A2*V - rlysis*I - mA*I
  # free viruses
  dV  = burst*rlysis*I - adsorption*A1*V - adsorption*A2*V - mV*V
  
  #Resource influx oscillates through time. A peak in resource influx occurs 
  #every few hundred time units
  dR = Rmax*sin(t/200)^n-cmax*A1*R/(R+h)-cmax*A2*R/(R+h)
  res = c(dA1,dD1,dA2,dD2,dI,dV,dR)
  list(res)
}
  )
}

#Explanation of parameters and values used


#Efficiency of converting resources to bacteria
e = 1

#Maximum consumption rate of resources by bacteria
cmax = 1

#Half-saturation constant of resource uptake by bacteria
h = 5

#Maximum rate at which active bacteria of strain i transition to dormancy.
#Strain w does not exhibit dormancy under these parameters.
rad1 = 0.5
rad2 = 0

#Maximum rate at which dormant bacteria of strain i resuscitate
rda1 = 0.1
rda2 = 0

#Steepness of transition rate's dependence on resources
w = 5

#Half-saturation constant for transition from active to dormant 
#or vice versa for genotype i
K1 = 5
K2 = 5 

## parameters added for phage infection.
## Using values used in NWs host-phage-models

# rate of phage adsorption to hosts (phi)
adsorption = 1e-4

# rate of lysis of infected individual (eta, "the length of an infection is exponentially distributed with a mean infection period 1/eta".Weitz QVE 2016, p76)
rlysis = 0.01

# virus burst size (beta)
burst = 50

# viral mortality (omega, originaly a chemostat washout parameter)
mV = 0.05

#####
# #Example plot showing the effect of this w parameter
# #Example vector of resources
# R1 = seq(from=0,to=10,by=0.1)
# 
# plot(R1,1-rad1*R1^2/(R1^2+K1^2),type="l",ylim=c(.5,1),xlab="Resources",ylab="Rate of Transition from Active to Dormant")
# points(R1,1-rad1*R1^4/(R1^4+K1^4),type="l",col="blue")
# points(R1,1-rad1*R1^8/(R1^8+K1^8),type="l",col="red")
# points(R1,1-rad1*R1^16/(R1^16+K1^16),type="l",col="green")
# legend(7,.9,lty=1,col=c("black","blue","red","green"),legend=c("w=2","w=4","w=6","w=8","w=16"))
# 
# plot(R1,rad1*R1^2/(R1^2+K1^2),type="l",ylim=c(0,.5),xlab="Resources",ylab="Rate of Transition from Dormant to Active")
# points(R1,rad1*R1^4/(R1^4+K1^4),type="l",col="blue")
# points(R1,rad1*R1^8/(R1^8+K1^8),type="l",col="red")
# points(R1,rad1*R1^16/(R1^16+K1^16),type="l",col="green")
# legend(7,.25,lty=1,col=c("black","blue","red","green"),legend=c("w=2","w=4","w=6","w=8","w=16"))
#####

#Background mortality rate of active bacteria or dormant bacteria
mA = 0.05
mD = 0.0005

#Maximum rate of inflow of resources
Rmax = 60


#This parameter determines the sharpness of resources influx peaks. This parameter
#should be even so that sin(t/200)^n is always non-negative and thus resource influx
#is never negative.
n = 4
#####
# #Plot to show the effect of the parameter n
# #Example time series
# times_ex = seq(from=0,to=10,by=.001)
# 
# par(mfrow=c(1,1))
# plot(times_ex,Rmax*sin(times_ex)^2,type="l",col="green")
# points(times_ex,Rmax*sin(times_ex)^8,type="l",col="red")
# points(times_ex,Rmax*sin(times_ex)^32,type="l",col="blue")
# legend(8.5,60,lty=1,col=c("green","red","blue"),legend=c("n=2","n=8","n=32"))
#####

#Define a vector of time values over which to simulate
times=seq(from=0,to=10000,by=1)


#Simulate the model. These simulations took about 6 minutes on my computer
Dorm_output1=ode(c(A1=1e4,D1=1,A2=1e4,D2=0,I=0,V=1,R=2),times,Dorm,
                c(e,cmax,h,rad1,rad2,rda1,rda2,w,K1,K2,mA,mD,Rmax,n=2, adsorption, rlysis, burst, mV),
                method="lsoda")
plot(Dorm_output1)
plot(Dorm_output1[,"A1"], Dorm_output1[,"V"], type = 'l', ylab = "V", xlab = "A1", main = "Phase Plane A1")
plot(Dorm_output1[,"A2"], Dorm_output1[,"V"], type = 'l', ylab = "V", xlab = "A2", main = "Phase Plane A2")


Dorm_output2=ode(c(A1=1,D1=0,A2=1,D2=0,I=0,V=1,R=2),times,Dorm,
                 c(e,cmax,h,rad1,rad2,rda1,rda2,w,K1,K2,mA,mD,Rmax,n=500, adsorption, rlysis, burst, mV),
                 method="lsoda")
plot(Dorm_output2)
plot(Dorm_output2[,"A1"], Dorm_output2[,"V"], type = 'l', ylab = "V", xlab = "A1", main = "Phase Plane A1")
plot(Dorm_output2[,"A2"], Dorm_output2[,"V"], type = 'l', ylab = "V", xlab = "A2", main = "Phase Plane A2")

#These two plots above show that the non-sporulating strain wins competition when
#n=2 while the sporulating strain wins competition when n=500

#Here I create a simulation-bifurcation diagram of the system by repeatedly 
#simulating with different n values

#Create array to hold results
Dorm_output_bifurcation=array(NA,dim=c(length(times),8,30))

#Populate each layer of the array with the results of a single simulation
for(i in 1:30){
  Dorm_output_bifurcation[,,i]=ode(c(A1=1,D1=1,A2=1,D2=0,I=0,V=1,R=2),times,Dorm,
                                  c(e,cmax,h,rad1,rad2,rda1,rda2,w,K1,K2,mA,mD,Rmax,n=2*i, adsorption, rlysis, burst, mV),
                                  method="lsoda")
}

#Calculate the frequency of each genotype averaged over the last 1000 time units (10%) 
#for a given simulation
Frequency1=array(NA,dim=c(30))
Frequency2=array(NA,dim=c(30))
for(i in 1:30){
  Frequency1[i]=mean((Dorm_output_bifurcation[9000:1000,2,i]+Dorm_output_bifurcation[9000:1000,3,i])/(Dorm_output_bifurcation[9000:1000,2,i]+Dorm_output_bifurcation[9000:1000,3,i]+Dorm_output_bifurcation[9000:1000,4,i]+Dorm_output_bifurcation[9000:1000,5,i]))
  Frequency2[i]=mean((Dorm_output_bifurcation[9000:1000,4,i]+Dorm_output_bifurcation[9000:1000,5,i])/(Dorm_output_bifurcation[9000:1000,2,i]+Dorm_output_bifurcation[9000:1000,3,i]+Dorm_output_bifurcation[9000:1000,4,i]+Dorm_output_bifurcation[9000:1000,5,i]))
}

#Plot the n parameter vs strain frequency
par(mfrow=c(1,1))
plot(2*seq(1,30,1),Frequency1,type="l",col="blue",xlab="Volatility parameter: n",ylab="Strain frequency",ylim=c(0,1))
points(2*seq(1,30,1),Frequency2,type="l",col="red")
legend(40,.5,col=c("blue","red"),lty=1,legend=c("Sporulator","Non-sporulator"))


#Save plot as a pdf
pdf('figures/Sporulation_Competition_.pdf', width = 6, height = 6, bg = "white")
plot(2*seq(1,30,1), Frequency1, type="l", col="blue", lwd = 2,
     xlab="Volatility parameter: n",
     ylab="Strain frequency",
     ylim=c(0,1))
points(2*seq(1,30,1), Frequency2, type="l", col="red", lwd = 2)
legend(35,.5, col=c("blue","red"), lty=1, lwd = c(2,2),
       legend=c("Sporulator","Non-sporulator"), bty = "n")
dev.off()

