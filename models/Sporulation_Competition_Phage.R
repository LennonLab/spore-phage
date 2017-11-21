#Clear environment
rm(list=ls())
#Necessary libraries
library(deSolve)

#Define dormancy and competition model. The differential equation is a function
#of time, a vector of state variables (y), and parameters
DormPhage = function(t,y,params){
  #Tell R what the 5 state variables are.
  A1 = y[1]; D1 = y[2]; A2 = y[3]; D2 = y[4]; R = y[5]; V = y[6]; I1 = y[7]; I2 = y[8];
  with(
    as.list(params),
    {
      #Define the equation governing the rate of change of each state variable
      dA1 = e*cmax*A1*R/(R+h)-rad1*(1-R^w/(R^w+K1^w))*A1+rda1*R^w/(R^w+K1^w)*D1-mA*A1 - phi*A1*V
      dD1 = rad1*(1-R^w/(R^w+K1^w))*A1-rda1*R^w/(R^w+K1^w)*D1-mD*D1
      dA2 = e*cmax*A2*R/(R+h)-rad2*(1-R^w/(R^w+K2^w))*A2+rda2*R^w/(R^w+K2^w)*D2-mA*A2 - phi*A2*V
      dD2 = rad2*(1-R^w/(R^w+K2^w))*A2-rda2*R^w/(R^w+K2^w)*D2-mD*D2
      
      #Resource influx oscillates through time. A peak in resource influx occurs 
      #every few hundred time units
      dR = Rmax*sin(t/200)^n-cmax*A1*R/(R+h)-cmax*A2*R/(R+h)
      
      # Phage attacks active cells
      dV = beta*eta*(I1 + I2) - phi*(A1+A2)*V - mV*V
      
      # infections
      dI1 <- phi*A1*V - eta*I1 - mA*I1
      dI2 <- phi*A2*V - eta*I2 - mA*I2
      
      # combine results
      res = c(dA1,dD1,dA2,dD2,dR,dV,dI1,dI2)
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


#Background mortality rate of active bacteria or dormant bacteria
mA = 0.05
mD = 0.0005

#Maximum rate of inflow of resources
Rmax = 60


#This parameter determines the sharpness of resources influx peaks. This parameter
#should be even so that sin(t/200)^n is always non-negative and thus resource influx
#is never negative.
n = 4

# Phage parameters
phi = 1e-7 # adsorption rate
eta = 0.1 # latency
beta = 200 # burst size
mV = 0.005 # phage mortality

#Define a vector of time values over which to simulate
times=seq(from=0,to=10000,by=1)


#Simulate the model. These simulations took about 6 minutes on my computer
DormPhage_output1=ode(c(A1=1e4,D1=1,A2=1e4,D2=0,R=2,V=1,I1=0,I2=0),times,DormPhage,
                 c(e,cmax,h,rad1,rad2,rda1,rda2,w,K1,K2,mA,mD,mV,Rmax,phi,eta,beta,n=2),
                 method="lsoda")
plot(DormPhage_output1)

DormPhage_output2=ode(c(A1=1e4,D1=1,A2=1e4,D2=0,R=2,V=0,I1=0,I2=0),times,DormPhage,
                 c(e,cmax,h,rad1,rad2,rda1,rda2,w,K1,K2,mA,mD,mV,Rmax,phi,eta,beta,n=500),
                 method="lsoda")
plot(DormPhage_output2)

#These two plots above show that the non-sporulating strain wins competition when
#n=2 while the sporulating strain wins competition when n=500

#Here I create a simulation-bifurcation diagram of the system by repeatedly 
#simulating with different n values

#Create array to hold results
DormPhage_output_bifurcation=array(NA,dim=c(length(times),9,30))

#Populate each layer of the array with the results of a single simulation
for(i in 1:30){
  DormPhage_output_bifurcation[,,i]=ode(c(A1=1e4,D1=1,A2=1e4,D2=0,R=2,V=1,I1=0,I2=0),times,DormPhage,
                                   c(e,cmax,h,rad1,rad2,rda1,rda2,w,K1,K2,mA,mD,mV,Rmax,phi,eta,beta,n=2*i),
                                   method="lsoda")
}

#Calculate the frequency of each genotype averaged over the last 1000 time units (10%) 
#for a given simulation
Frequency1=array(NA,dim=c(30))
Frequency2=array(NA,dim=c(30))
for(i in 1:30){
  Frequency1[i]=mean((DormPhage_output_bifurcation[9000:1000,2,i]+DormPhage_output_bifurcation[9000:1000,3,i])/(DormPhage_output_bifurcation[9000:1000,2,i]+DormPhage_output_bifurcation[9000:1000,3,i]+DormPhage_output_bifurcation[9000:1000,4,i]+DormPhage_output_bifurcation[9000:1000,5,i]))
  Frequency2[i]=mean((DormPhage_output_bifurcation[9000:1000,4,i]+DormPhage_output_bifurcation[9000:1000,5,i])/(DormPhage_output_bifurcation[9000:1000,2,i]+DormPhage_output_bifurcation[9000:1000,3,i]+DormPhage_output_bifurcation[9000:1000,4,i]+DormPhage_output_bifurcation[9000:1000,5,i]))
}

#Plot the n parameter vs strain frequency
par(mfrow=c(1,1))
plot(2*seq(1,30,1),Frequency1,type="l",col="blue",xlab="Volatility parameter: n",ylab="Strain frequency",ylim=c(0,1))
points(2*seq(1,30,1),Frequency2,type="l",col="red")
legend(40,.5,col=c("blue","red"),lty=1,legend=c("Sporulator","Non-sporulator"))


#Save plot as a pdf
pdf('figures/Sporulation_Competition_Phage.pdf', width = 6, height = 6, bg = "white")
plot(2*seq(1,30,1), Frequency1, type="l", col="blue", lwd = 2,
     xlab="Volatility parameter: n",
     ylab="Strain frequency",
     ylim=c(0,1))
points(2*seq(1,30,1), Frequency2, type="l", col="red", lwd = 2)
legend(35,.5, col=c("blue","red"), lty=1, lwd = c(2,2),
       legend=c("Sporulator","Non-sporulator"), bty = "n")
dev.off()

