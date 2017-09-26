require(Matrix)
require(numDeriv)

# Functions to simulate SI/SIS models, numerically
# calculate solution to master equation and calculate
# likelihood based on final size after a fixed exposure time 

SIforward<-function(io, ntot, beta, tend)
{
  t <- 0
  s <- (ntot - io)
  i <- io
  while(t < tend & i < ntot)	
  {
   rate <- (beta * s * i) / ntot
   tnext <- rexp(1,rate)
   t <- t + tnext
   if(t < tend)
   {
   	 i = i + 1
   	 s = s - 1
   }
  	
  }
  return(i)

}

SISforward<-function(io, ntot, beta, u, tend)
{
  t <- 0
  s <- (ntot - io)
  i <- io
 rate <- (beta * s * i) + u*i
  while(t < tend & rate > 0)	
  {
   
   tnext <- rexp(1,rate)
   t <- t + tnext
   if(t < tend)
   {
     nexte=runif(1,0,rate)
     
     if(nexte < (beta * s * i))
     {
	 # Infection
   	 i = i + 1
   	 s = s - 1
   	 }else
   	 {
   	 # Replacement
   	 i = i - 1
   	 s = s + 1
   	 }
   	 
   }
  	rate <- (beta * s * i) + u*i
  }
  return(i)

}


generate_SIp <- function(l,n)
{

 SIp <-matrix(0,n,1)
 SIp[l] <- 1
 return((SIp))
	
}

generate_SISp <- function(l,n)
{

 SIp <-matrix(0,n+1,1)
 SIp[l+1] <- 1
 return((SIp))
	
}


generate_SIQ <- function(beta_t,N)
{
 Q <- matrix(0,N,N)
 for(i in 1:(N-1))
 {
 	Q[i,i] = -beta_t*(i)*(N-i)/N
 	Q[i+1,i] = beta_t*(i)*(N-i)/N
}
return(Q)
}

generate_SISQ <- function(beta,u,N)
{
 Q <- matrix(0,N+1,N+1)
 for(j in 1:(N+1))
 {
 	i = j-1
 	s = N-i
 	Q[j,j] = -beta*(i)*(s) - u*i
 	if(j<(N+1)){
 	Q[j,j+1] = beta*(i)*(s)}
 	if(j>1){
 	Q[j,j-1] = u*i}
 	
 	
}
return(t(Q))
}

posteriorSI <- function(aito)
{
  T = aito[1]
 
  if( T < 0.0)
 {
 return(-Inf)	
 }
 
  Q <- generate_SIQ(T,N)
  pSI <- generate_SIp(iota,N)
  Evolve <- expm(Q)
  pSI <- Evolve %*% pSI
 
  accum = 0
 
  for(i in 1:length(expI))
  {
  	 
    accum = accum + log(pSI[expI[i]])
    
  	
  }
  
  return(accum)
  
}

solveSIS <- function(beta,N,mu,tau,iota)
{
  Q <- generate_SISQ(beta,mu,N)
  pSI <- generate_SISp(iota,N)
  Evolve <- expm(Q*tau)
  pSI <- Evolve %*% pSI
 
  return(pSI)
  
}

solveSISk <- function(beta,N,mu,tau,iota,k)
{

return(solveSIS(beta,N,mu,tau,iota)[k])

}


fisher_info <- function(beta,N,mu,tau,iota)
{

 fisher = 0
 for(k in 1:(N+1))
 {
 pk <- solveSISk(beta,N,mu,tau,iota,k)
 partial_thetas <- grad(solveSISk,beta,N=N,mu=mu,tau=Tau,iota=iota,k=k)
 fisher = fisher + (1/pk)*partial_thetas*partial_thetas
 
 }
 	
 return(fisher)
 
}

likelihoodSIS <- function(aito,datap)
{
  beta = aito[1]
  N = datap[1]
  mu = datap[2]
  tau = datap[3]
  iota = datap[4]
  expI = datap[5]
  
  if( beta < 0.0 || expI > N)
  {
  return(-Inf)	
  }
 
  Q <- generate_SISQ(beta,mu,N)
  pSI <- generate_SISp(iota,N)
  Evolve <- expm(Q*tau)
  pSI <- Evolve %*% pSI
 
  return(pSI[expI+1])
  
}

