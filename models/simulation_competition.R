rm(list = ls())

library(tidyverse)

#Simulating competition for resource between spore fromer and non spore former under resource fluctutations

#=================
# update function
#=================
g.plus1 <- function(step,r,A1,D1,A2){
  #debug at step i
  if (dbg==step) browser()
  # update resources
  #------------------
  #resource washout
  r <- r * r.dilut
  # periodic resource injection
  if (step %% r.step == 0)
    r <- r+r.input
  
  # Bacterial death
  #------------------
  # death of active bacteria 1
  A1 <- death(A1, A1.death )
  # death of dormant bacteria 1
  D1 <- death(D1, D1.death )
  # death of active bacteria 2
  A2 <- death(A2, A2.death )
  
  
  # Dormancy and resucitation
  # -----------------------------
  # resource threshold entrance to dormancy
  D1.new <- dorm(r = r,A = A1, D = D1)
  
  # Bacterial reproduction
  #------------------------
  # partition resource by relative abundance of active bacteria
  # newly dormant cells do not count for resouce partitioning
  A1.accounted <- ifelse (D1.new > 0, A1 - D1.new, A1)
  # ifelse condition to prevent division by 0
  rel.A1 <-ifelse(A1.accounted+A2, A1.accounted/(A1.accounted+A2),0.5)
  r1 <- r*rel.A1
  r2 <- r-r1
  # bacteria 1 # 
  A1.new <- birth(r = r1, N = A1, d = A1.doubling )
  # bacteria 2 # 
  A2.new <- birth(r = r2, N = A2, d = A2.doubling )
  
  # updating variables
  #-------------------
  A1 <- A1 + A1.new - D1.new
  D1 <- D1 + D1.new
  A2 <- A2 + A2.new
  r <- r - A1.new - A2.new
  # force fractions to zero, no atto-fox!
  # if (r<1) r <- 0
  if (A1<1) A1 <- 0
  if (D1<1) D1 <- 0
  if (A2<1) A2 <- 0
  
  
  return(list(r,A1,D1,A2))
}

#=================
# processes functions
#=================
death <- function(N, m){
  return (N-m*N)
}

birth <- function(r, N, d){
  # bacteria require a single and whole resource unit to double
  r <- floor(min(N,r))
  return (r/d)
}

dorm <- function (r, A, D){
  # resource threshold entrance to dormancy
  d.new <- 0
  if (r<r.low & A>0)
    d.new <- d.fract * A
  # resource threshold exit from dormancy
  if (r>r.high & D>0)
    d.new <-  -1 * d.fract * D
  return(d.new)
}
#=================
# parameters
#=================

# define time steps
steps <- 1000

#############
# resources #
#############

# r.input resource are injected every r.step steps
r.input <- 1000
r.step <- 40

# constant resource dilution by washout (well mixed)
# r.dilut <- 0.98
# currently in this model resources are only removed by use
r.dilut <- 1

##############
# bacteria 1 #
##############

# doubling time (in simulation steps)
A1.doubling <- 2
# bacterial mortality
A1.death <- 0.05
# initial baterial pop size
A1.initial <- 1000

##############
# bacteria 2 #
##############

# doubling time (in simulation steps)
A2.doubling <- 2
# bacterial mortality
A2.death <- 0.05
# initial baterial pop size
A2.initial <- 1000

##############
# dormant 1  #
##############
D1.initial <-  0
D1.death <- 1e-4 * A1.death
#resource threshold to enter dormancy
r.low <- 0.1 * r.input
#resource threshold to exit dormancy
r.high <- 0.5 * r.input
# frction of cells entering or exiting dormancy when thresholds are met
d.fract <- 0.5


#==============
#  THE MODEL  #
#==============

# results matrix
g <- matrix(nrow = steps, ncol = 4)
colnames(g) <- c('r','A1','D1', 'A2')

#initialise
g[1,'r'] <- r.input
g[1,'A1'] <- A1.initial
g[1,'D1'] <- D1.initial
g[1,'A2'] <- A2.initial

dbg <- 0 #debug update function at this step
r.step <- 25
for (i in 2:nrow(g)){
  g[i,] <- unlist(g.plus1(step = i,
                          r = as.numeric(g[i-1,'r']),
                          A1 = as.numeric(g[i-1,'A1']), 
                          D1 =as.numeric(g[i-1,'D1']),  
                          A2=as.numeric(g[i-1,'A2'])))
}

matplot(1:nrow(g), g, type = 'l', log='y', ylab = 'abundance (log-scale)', xlab = 'step', main=paste("resoucre replenished every",r.step,'time steps'), ylim = c(0.1,5e3))
grid()
legend("bottomright", colnames(g), lty=1:ncol(g), col=1:ncol(g))

# # zoom in on some part of the plot
# s <-190; e <- 205
# matplot(s:e, g[s:e,], type = 'l', log='')
# grid()
# abline(h = c(r.low,r.high), lty=2)
# legend("right", colnames(g), lty=1:ncol(g), col=1:ncol(g))

# # phase plane plots
# matplot(g[,'r'], g[,c('A1', 'A2', 'D1')], type = 'l', log='')
# grid()
# legend("topright", colnames(g[,c('A1', 'A2', 'D1')]), lty=1:ncol(g[,c('A1', 'A2', 'D1')]), col=1:ncol(g[,c('A1', 'A2', 'D1')]))
# 
# plot(g[,'A2'], g[,'A1']+g[,'D1'], type = 'l')


#Compete species over a range of resource pulse periods
R.STEPS <- matrix(nrow = 50, ncol = 3 )
colnames(R.STEPS) <- c('r.step','A1+D1','A2')
R.STEPS[,'r.step'] <- floor(seq(1,50,along.with = R.STEPS[,'r.step']))
for (j in 1:length(R.STEPS[,'r.step'])){
  r.step <- R.STEPS[j,'r.step']
  for (i in 2:nrow(g)){
    g[i,] <- unlist(g.plus1(step = i,r = g[i-1,'r'], A1 = g[i-1,'A1'], D1 =g[i-1,'D1'],  A2=g[i-1,'A2']))
  }
  # average species abundance over last 10% of time steps
  R.STEPS[j,'A1+D1'] <- mean(rowSums(g[(0.9*steps):steps, c('A1','D1')]))
  R.STEPS[j,'A2'] <-  mean(g[(0.9*steps):steps, 'A2'])
}


matplot(R.STEPS[,'r.step'], R.STEPS[,2:3], type = 'l', log='y',
        xlab = 'resource replenish period', ylab = 'abundance')
grid()
legend("topright", colnames(R.STEPS[,2:3]), lty=1:2, col=1:2)


# plot with gglot
d <- gather(as.data.frame(R.STEPS), key = 'sp', value = 'abundance', 'A1+D1', 'A2')
ggplot(d, aes(x=r.step, y=abundance, fill=sp))+
  xlab('resource replenish period')+ ylab('relative abundance')+xlim(1,50)+
  geom_area(position='fill', alpha=0.5)+theme_bw()

ggplot(d, aes(x=r.step, y=abundance, col=sp))+
  geom_line()+theme_bw()

