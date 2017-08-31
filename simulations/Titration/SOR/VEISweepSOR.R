# Load ABC posterior distributions
samples <- read.csv('../../../parameters/samplesSOR.csv')[,2:12]
weight  <- scan('../../../parameters/weightSOR.dat')

reps <- 5000

require('parallel')

PTI<-scan('../../../data/DarthPTI.csv')
HerdSize<-scan('../../../data/DarthHerds.csv')

make_unique_filename <- function(x)
{
	return(sprintf('output%d.csv',x))
}

make_unique_filenameFu <- function(x)
{
	return(sprintf('output%d.csvFu',x))
}

single_runIndivSORI <- function(vacc_target, batch_interval,retain_reactors,vacc_effS,vacc_effI,DIVA, Immunity=1.0/1.650)
{

trial_duration = 5
batch_interval = 180

herd_no = sample(which(PTI==0),1)-1
params <- samples[which(rmultinom(1,prob=weight,size=1)==1),]
ranx <- round(runif(1,0,.Machine$integer.max))
filename <- make_unique_filename(ranx)
filenameFu <- make_unique_filenameFu(ranx)
print(filename)


#DIVA = c(73.3/100.0,aito[2])

#DIVA = c(aito[1],1-99.85/100.0)
#DIVA = c(aito[1],aito[2])

doVacc=1
suspend_break=0
Eligible=0.0

jazzle<-system(sprintf('../../../bin/bTBField %g %g %g %g %g %g %g %g %g %g %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %d %s %g %d %d %d %g %g %g %d %d',
DIVA[1], DIVA[2],DIVA[1],DIVA[2],
0.0,0.0,DIVA[1],DIVA[2],
1.0-vacc_effS,1.0,1,
params[5],-1.0,Immunity,1.0,params[5],params[5],-1.0,
params[6],params[6],params[6],(1.0-vacc_effI)*(params[6]),(1.0-vacc_effI)*(params[6]),(1.0-vacc_effI)*(params[6]),(1.0-vacc_effI)*(params[6]),params[7],
params[8],params[9],params[10],
1.80,2.3,1.5,1.3,
0.0,0.0,0.0,
params[11], 
24, 0.0,1.0,1,filename,batch_interval,doVacc,1,suspend_break,Eligible,vacc_target,trial_duration,retain_reactors,herd_no),intern=F)

#disco<-read.csv('TestFuDiscloseH0.dat',sep=' ')
disco <- read.csv(sprintf('output%d.csvFuDiscloseH0.dat',ranx),sep=' ')
slaughter <- try(read.csv(sprintf('output%d.csvFuReactorsSH0.dat',ranx),sep=' ',header=FALSE),silent=TRUE)
Fu <- read.csv(filenameFu,sep=',')
indiv <- try(read.csv(sprintf('output%d.csvFuIndividual0.dat',ranx)),silent=TRUE)

system(sprintf('rm %s',filename))
system(sprintf('rm %s*',filenameFu))
system(sprintf('rm *.dat',filenameFu))

agg=data.frame(ranx=ranx,Herd=herd_no,Vaccinated=0,Reactor=sum(indiv$Reactor[indiv$Vaccinated==0]), VL=sum(indiv$VL[indiv$Vaccinated==0]),AtRisk=sum(indiv$Vaccinated==0))
agg=rbind(agg,data.frame(ranx=ranx,Herd=herd_no,Vaccinated=1,Reactor=sum(indiv$Reactor[indiv$Vaccinated==1]), VL=sum(indiv$VL[indiv$Vaccinated==1]),AtRisk=sum(indiv$Vaccinated==1)))

return(agg)
}

DIVA1 = c(73.3/100.0,1-99.85/100.0)
DIVA2 = c(0.25,0.0)
DIVA3 = c(1.0,1-99.85/100.0)
DIVA4 = c(1.0,0.0)

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(0.0,0,0.0,0.0,DIVA1,(1.0/1.650))},mc.set.seed=TRUE,mc.preschedule=FALSE,mc.cores=4) 
SimsDf = cbind(as.data.frame(do.call(rbind,outlist)),immune=(1.0/1.650),remove=T,treatment=F)

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(0.50,0,0.5,0.5,DIVA1,(1.0/1.650))},mc.set.seed=TRUE,mc.preschedule=FALSE,mc.cores=4) 
SimsDf = rbind(SimsDf, cbind(as.data.frame(do.call(rbind,outlist)),immune=(1.0/1.650),remove=T,treatment=T))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(0.50,0,0.5,0.5,DIVA1,(1.0/1.650)*3)},mc.set.seed=TRUE,mc.preschedule=FALSE,mc.cores=4) 
SimsDf = rbind(SimsDf, cbind(as.data.frame(do.call(rbind,outlist)),immune=(1.0/1.650)*3,remove=T,treatment=T))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(0.50,0,0.5,0.5,DIVA1,(1.0/1.650)*6)},mc.set.seed=TRUE,mc.preschedule=FALSE,mc.cores=4) 
SimsDf = rbind(SimsDf, cbind(as.data.frame(do.call(rbind,outlist)),immune=(1.0/1.650)*6,remove=T,treatment=T))


outlist <- mclapply(1:reps,function(x){single_runIndivSORI(0.0,1,0.0,0.0,DIVA1,(1.0/1.650))},mc.set.seed=TRUE,mc.preschedule=FALSE,mc.cores=4) 
SimsDfR = cbind(as.data.frame(do.call(rbind,outlist)),immune=(1.0/1.650),remove=F,treatment=F)

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(0.50,1,0.5,0.5,DIVA1,(1.0/1.650))},mc.set.seed=TRUE,mc.preschedule=FALSE,mc.cores=4) 
SimsDfR = rbind(SimsDfR, cbind(as.data.frame(do.call(rbind,outlist)),immune=(1.0/1.650),remove=F,treatment=T))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(0.50,1,0.5,0.5,DIVA1,(1.0/1.650)*3)},mc.set.seed=TRUE,mc.preschedule=FALSE,mc.cores=4) 
SimsDfR = rbind(SimsDfR, cbind(as.data.frame(do.call(rbind,outlist)),immune=(1.0/1.650)*3,remove=F,treatment=T))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(0.50,1,0.5,0.5,DIVA1,(1.0/1.650)*6)},mc.set.seed=TRUE,mc.preschedule=FALSE,mc.cores=4) 
SimsDfR = rbind(SimsDfR, cbind(as.data.frame(do.call(rbind,outlist)),immune=(1.0/1.650)*6,remove=F,treatment=T))

save.image('ImmunityRetain.RData')