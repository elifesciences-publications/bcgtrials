resample <- function(samples,weight)
{

pick <- replicate(length(weight),which(rmultinom(1,prob=weight,size=1)==1))

return(samples[pick,])

}


require(ggplot2)

df <- data.frame(model='Barlow 1997',beta=(3.0e-5),q=0.0,TO=22.0,TR=180)
df <- rbind(df,data.frame(model='Fischer 2005',beta=(6e-3),q=1.0,TO=41.17647,TR=233.3333))
df <- rbind(df,data.frame(model='O\'Hare 2014',beta=(2.747253e-05),q=0.0,TO=100.0,TR=190.0))

samples <- read.csv('../../../parameters/samplesSOR2012.csv')[,2:13]
weight  <- scan('../../../parameters/weightSOR2012.dat')

samples<-resample(samples,weight)

df <- rbind(df,data.frame(model='Conlan SOR',beta=median(samples[,6])/364.0,q=median(samples[,7]),TO=364.0*median(samples[,5]),TR=NA))

samples <- read.csv('../../../parameters/samplesSORI2012.csv')[,2:14]
weight  <- scan('../../../parameters/weightSORI2012.dat')


samples<-resample(samples,weight)

df <- rbind(df,data.frame(model='Conlan SORI',beta=median(samples[,6])/364.0,q=median(samples[,7]),TO=364.0*median(samples[,5]),TR=364.0*median(samples[,13])))


samples <- read.csv('../../../parameters/samplesSOR.csv')[,2:12]
weight  <- scan('../../../parameters/weightSOR.dat')

samples<-resample(samples,weight)

df <- rbind(df,data.frame(model='Conlan SORD',beta=median(samples[,6])/364.0,q=median(samples[,7]),TO=364.0*median(samples[,5]),TR=NA))

samples <- read.csv('../../../parameters/samplesSORI.csv')[,2:13]
weight  <- scan('../../../parameters/weightSORI.dat')


samples<-resample(samples,weight)

df <- rbind(df,data.frame(model='Conlan SORID',beta=median(samples[,6])/364.0,q=median(samples[,7]),TO=364.0*median(samples[,5]),TR=364.0*median(samples[,12])))


ExperimentR0 <- function(incontact, beta, q, H, Hm)
{
return(incontact*beta*H/((H/Hm)^q))
}

SORIExperimentR0 <- function(incontact, beta, q, H, Hm,so,sr)
{

return((beta*H/((H/Hm)^q))
* ((incontact+(so/(sr*(so-sr)))*(exp(-sr*incontact)-1) - (sr/(so*(so-sr)))*(exp(-so*incontact)-1))))

}

df2 = data.frame(Initial='Infective',Model=df$model[1],H=seq(1,200,1),R0 = ExperimentR0(364.0,df$beta[1],df$q[1],seq(1,200,1),1))
df2 = rbind(df2,data.frame(Initial='Infective',Model=df$model[2],H=seq(1,200,1),R0 = ExperimentR0(364,df$beta[2],df$q[2],seq(1,200,1),1)))
df2 = rbind(df2,data.frame(Initial='Infective',Model=df$model[3],H=seq(1,200,1),R0 = ExperimentR0(364.0,df$beta[3],df$q[3],seq(1,200,1),1)))

for(i in 4:7)
{

df2 = rbind(df2,data.frame(Initial='Infective',Model=df$model[i],H=seq(1,200,1),R0 = ExperimentR0(364.0,df$beta[i],df$q[i],seq(1,200,1),165)))

}

df2 = rbind(df2,data.frame(Initial='Occult',Model=df$model[1],H=seq(1,200,1),R0 = SORIExperimentR0(364.0,df$beta[1],df$q[1],seq(1,200,1),1,1/df$TO[1],1/df$TR[1])))
df2 = rbind(df2,data.frame(Initial='Occult',Model=df$model[2],H=seq(1,200,1),R0 = SORIExperimentR0(364,df$beta[2],df$q[2],seq(1,200,1),1,1/df$TO[2],1/df$TR[2])))
df2 = rbind(df2,data.frame(Initial='Occult',Model=df$model[3],H=seq(1,200,1),R0 = SORIExperimentR0(364.0,df$beta[3],df$q[3],seq(1,200,1),1,1/df$TO[3],1/df$TR[3])))

for(i in c(4,6))
{

df2 = rbind(df2,data.frame(Initial='Occult',Model=df$model[i],H=seq(1,200,1),R0 = ExperimentR0(364.0,df$beta[i],df$q[i],seq(1,200,1),165)))

}

for(i in c(5,7))
{

df2 = rbind(df2,data.frame(Initial='Occult',Model=df$model[i],H=seq(1,200,1),R0 = SORIExperimentR0(364.0,df$beta[i],df$q[i],seq(1,200,1),165,1/df$TO[i],1/df$TR[i])))

}


pdf('CompareR.pdf',family='Helvetica',points=12,height=5/2,width=10/2)
print(ggplot(df2,aes(x=H,y=R0,col=Model,linetype=Model))+geom_line()+xlab('Group Size')+ylab('R')+facet_wrap(~Initial)+scale_linetype_discrete(name="Model",labels=c("Barlow 1997", "Fischer 2005", "O\'Hare 2014", "SOR 2012", "SORI 2012", "SOR 2015", "SORI 2015"))+scale_colour_discrete(name="Model",labels=c("Barlow 1997", "Fischer 2005", "O\'Hare 2014", "SOR 2012", "SORI 2012", "SOR 2015", "SORI 2015")))
dev.off()