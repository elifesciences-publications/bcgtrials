source('SISLikelihood.r')
require(parallel)
require(coda)
require(MASS)

mu = 0.0

# In contact time of 10 months
Tau = 364.0*(10/12)

# Number of infected animals at end of in contact period in each group of ten animals
vals=c(3,1,1,0,0,0,0,0,2,1)

lik_it <- function(beta)
{
beta=exp(beta)
if(beta < 0.0 | beta > 0.1)
{return(-Inf)}
else{
p=solveSIS(beta,10,0.0,Tau,4)
l=p[5:11]
sum(log(l[vals+1]))
}
}

prior_prob <- function(x)
{

x=exp(x)

pgamma(x,1,1)

}

# Maximum Likelihood estimate
mleTrans=optimise(lik_it,log(c(0.0000001,0.1)),maximum=T)$maximum

# Estimate Bayesian posterior distribution by MCMC
# Code adapted from skeleton provided by TJ McKinley

n = 100000
accept = 0
npar = 1
# set up output vector of length 'niter' to record chain
# (append extra column for unnormalised posterior)

# Mixing proportion for adaptive MCMC
scale = 0.05

psamples=matrix(NA,n,npar+1)

x=mleTrans

        lik=lik_it(x)
		acc.curr = (lik)-log(prior_prob(x))

        psamples[1,]=c(x,acc.curr)
		pars = x
		
prop.var.ini <- diag(npar) * (0.1 ^ 2) / npar
        
         # run chain
    	nacc <- 0
        for (i in 2:n) {
         				# propose new value
        				if (i <= (2 * npar)) {pars.prop <- mvrnorm(1, pars, prop.var.ini)
        				}else
        				{
            			u <- runif(1, 0, 1)
            			if (u < scale) {pars.prop <- mvrnorm(1, pars, prop.var.ini) 
            			}else{pars.prop <- mvrnorm(1, pars, prop.var)}
        				
        				}
        		if(prior_prob(pars.prop)!=0)
        		{
                oldlik=lik
                lik=lik_it(pars.prop)
                
                acc.prop=(log(prior_prob(pars.prop))+(lik))
               
               #accept-reject proposal
        		if(is.finite(acc.prop))
       			 {
          	 	 	 acc <- acc.prop -oldlik-log(prior_prob(pars))
           			 acc <- exp(acc)
           		 	 u <- runif(1, 0, 1)
           		 	if (u < acc)
            		{
                	pars <- pars.prop
                	acc.curr <- acc.prop
                	nacc <- nacc + 1
            		}
        		}
               	}
               	
                psamples[i,]=c(pars, acc.curr)

	    # calculations for adaptive proposal 
        if (i >= (2 * npar)) prop.var <- var(psamples[1:i, 1:npar]) * (2.38 ^ 2) / npar
        
        # print some output to screen for book-keeping
        if (i %% 100 == 0)
        {
            print(paste("i =", i, "acceptance rate =", nacc / 100))
            nacc <- 0
            #save.image('PostTransMCMC.RData')

        }

        }
        
        
psamples[,1]=exp(psamples[,1])
pm <- mcmc(data=psamples,thin=1)   
     
#Posterior comparison to mle     
truehist(log10(pm[,1]),col='black')    
abline(v=log10(exp(mleTrans)),col='red',lwd=2)     
   
ppd <- function(beta)
{
p=solveSIS(beta,10,0.0,Tau,4)
sum(p[5:11]*seq(0,6))
}   

truehist(sapply(pm[,1],ppd),col='black')
    
save.image('PostTransMCMC.RData')
        