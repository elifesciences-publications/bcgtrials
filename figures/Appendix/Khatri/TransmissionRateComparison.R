# Load ABC posterior distributions
samples <- read.csv('../../../parameters/samplesSOR.csv')[,2:12]
weight  <- scan('../../../parameters/weightSOR.dat')

samples = resample(samples,weight)

beta_t = (samples[,6])/364.0
q = samples[,7]

mu=scan('../../../data/DarthTurnover.csv')

ExperimentR0 <- function(incontact, beta, q, H, Hm)
{
return(incontact*beta*H/((H/Hm)^q))
}

Hm=165

herd_size = seq(10,400,5)

R0 = sapply(herd_size,function(N){summary(ExperimentR0(10*30,beta_t,q,N,Hm))})
matplot(herd_size,t(R0),type='l',ylim=c(0,4))

N=10
betaPred10 = sapply(1:dim(samples)[1],function(i){beta_t[i]/(N/Hm)^(q[i])})

N=50
betaPred50 = sapply(1:dim(samples)[1],function(i){beta_t[i]/(N/Hm)^(q[i])})


N=100
betaPred100 = sapply(1:dim(samples)[1],function(i){beta_t[i]/(N/Hm)^(q[i])})

N=Hm
betaPredHm = sapply(1:dim(samples)[1],function(i){beta_t[i]*(N/Hm)^(q[i])})

Nupper=390
betaPredUpper = sapply(1:dim(samples)[1],function(i){beta_t[i]*(N/Hm)^(q[i])})


 load('PostTransMCMC.RData')
 
 cfactor <- 0.0393700787402

 
pdf('KhatriEstimate.pdf',width=(135)*cfactor,height=(135/1.75)*cfactor,family='Helvetica',pointsize=8) 
par(mfrow=c(1,2)) 
plot(density(log10(pm[,1])),ylim=c(0,0.15),xlab=expression(log10(beta)),main='Experimental Estimate (SI model)')
abline(v= log10(exp(mleTrans)),lty=2,lwd=2)
legend('topleft',c('Posterior','MLE'),col=c('black'),lwd=2,lty=c(1,2))

plot(density((betaPred10*10)),col='red',lwd=2,xlab=expression(beta),main='Field Estimates (ABC)',xlim=c(0,0.0025))
lines(density((betaPred50*50)),col='green',lwd=2)
lines(density((betaPredHm*Hm)),col='blue',lwd=2)
abline(v= (exp(mleTrans)),lty=2,lwd=2)
legend('topright',c('SOR (H=10)','SOR (H=50)','SOR (H=165)'),col=c('red','green','blue'),lwd=2)
dev.off()
