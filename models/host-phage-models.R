require(deSolve)
require(rootSolve)

HPmod <- function(t, x, parms){
  with(as.list(c(parms, x)), {
    
    dN <- r*N*(1-(N+I)/K) - phi*N*V - omega*N
    dI <- phi*N*N - eta*I - omega*I
    dV <- beta*eta*I - phi*N*V - omega*V
    res <- c(dN, dI, dV)
    return(list(res))
  })
}

times <- seq(0, 1000, by = .1)
x0 <- c(N = 10e4, I = 0, V = 1)
parms <- c(eta = 1, r = 0.16, K = 2.2e7, phi = 10e-9, beta = 50, omega = 0.05)

out <- lsoda(x0, times, HPmod, parms)
plot(out)
plot(out[,"N"], out[,"V"], type = 'l', ylab = "V", xlab = "N", main = "Phase Plane")

stode(y = (tail(out, n = 1)[,-1]), func = HPmod, parms = parms, pos = T)
