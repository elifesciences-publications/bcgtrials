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


single_runIndivSORI <- function(trial_duration,vacc_target, batch_interval,retain_reactors,vacc_effS,vacc_effI,DIVA)
{

herd_no = sample(which(PTI==0),1)-1

aito <- samples[which(rmultinom(1,prob=weight,size=1)==1),]

params <- log10(aito)

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
10^params[5],-1.0,1.0/1.650,1.0,10^params[5],10^params[5],-1.0,
10^params[6],10^params[6],10^params[6],(1.0-vacc_effI)*10^params[6],(1.0-vacc_effI)*10^params[6],(1.0-vacc_effI)*10^params[6],(1.0-vacc_effI)*10^params[6],10^params[7],
10^params[8],10^params[9],10^params[10],
1.80,2.3,1.5,1.3,
0.0,0.0,0.0,
10^params[11], 
24, 0.0,1.0,1,filename,batch_interval,doVacc,1,suspend_break,Eligible,vacc_target,trial_duration,retain_reactors,herd_no),intern=F)

#disco<-read.csv('TestFuDiscloseH0.dat',sep=' ')
disco <- read.csv(sprintf('output%d.csvFuDiscloseH0.dat',ranx),sep=' ')
slaughter <- try(read.csv(sprintf('output%d.csvFuReactorsSH0.dat',ranx),sep=' ',header=FALSE),silent=TRUE)
Fu <- read.csv(filenameFu,sep=',')
#indiv <- try(read.csv(sprintf('output%d.csvFuIndividual0.dat',ranx)),silent=TRUE)

system(sprintf('rm %s',filename))
system(sprintf('rm %s*',filenameFu))
system(sprintf('rm *.dat',filenameFu))
return(Fu)
}

DIVA1 = c(73.3/100.0,1-99.85/100.0)
DIVA2 = c(0.25,0.0)
DIVA3 = c(1.0,1-99.85/100.0)
DIVA4 = c(1.0,0.0)

ICDIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.0,180,0,0.0,0.0,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 

VES30VEI0IVP0DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.0,180,0,0.3,0.0,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI0IVP10DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.1,180,0,0.3,0.0,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI0IVP25DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.25,180,0,0.3,0.0,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI0IVP50DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.50,180,0,0.3,0.0,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI0IVP75DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.75,180,0,0.3,0.0,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI0IVP90DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.90,180,0,0.3,0.0,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI0IVP100DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,1.0,180,0,0.3,0.0,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 

VES30VEI30IVP0DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.0,180,0,0.3,0.3,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI30IVP10DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.1,180,0,0.3,0.3,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI30IVP25DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.25,180,0,0.3,0.3,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI30IVP50DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.50,180,0,0.3,0.3,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI30IVP75DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.75,180,0,0.3,0.3,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI30IVP90DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.90,180,0,0.3,0.3,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI30IVP100DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,1.0,180,0,0.3,0.3,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 

VES30VEI60IVP0DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.0,180,0,0.3,0.6,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI60IVP10DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.1,180,0,0.3,0.6,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI60IVP25DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.25,180,0,0.3,0.6,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI60IVP50DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.50,180,0,0.3,0.6,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI60IVP75DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.75,180,0,0.3,0.6,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI60IVP90DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.90,180,0,0.3,0.6,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI60IVP100DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,1.0,180,0,0.3,0.6,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 

VES30VEI90IVP0DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.0,180,0,0.3,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI90IVP10DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.1,180,0,0.3,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI90IVP25DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.25,180,0,0.3,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI90IVP50DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.50,180,0,0.3,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI90IVP75DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.75,180,0,0.3,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI90IVP90DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.90,180,0,0.3,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES30VEI90IVP100DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,1.0,180,0,0.3,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 

ICDIVA1df = as.data.frame(do.call(rbind,ICDIVA1))

VES30VEI0IVP0DIVA1df = as.data.frame(do.call(rbind,VES30VEI0IVP0DIVA1))
VES30VEI0IVP10DIVA1df = as.data.frame(do.call(rbind,VES30VEI0IVP10DIVA1))
VES30VEI0IVP25DIVA1df = as.data.frame(do.call(rbind,VES30VEI0IVP25DIVA1))
VES30VEI0IVP50DIVA1df = as.data.frame(do.call(rbind,VES30VEI0IVP50DIVA1))
VES30VEI0IVP75DIVA1df = as.data.frame(do.call(rbind,VES30VEI0IVP75DIVA1))
VES30VEI0IVP90DIVA1df = as.data.frame(do.call(rbind,VES30VEI0IVP90DIVA1))
VES30VEI0IVP100DIVA1df = as.data.frame(do.call(rbind,VES30VEI0IVP100DIVA1))

VES30VEI30IVP0DIVA1df = as.data.frame(do.call(rbind,VES30VEI30IVP0DIVA1))
VES30VEI30IVP10DIVA1df = as.data.frame(do.call(rbind,VES30VEI30IVP10DIVA1))
VES30VEI30IVP25DIVA1df = as.data.frame(do.call(rbind,VES30VEI30IVP25DIVA1))
VES30VEI30IVP50DIVA1df = as.data.frame(do.call(rbind,VES30VEI30IVP50DIVA1))
VES30VEI30IVP75DIVA1df = as.data.frame(do.call(rbind,VES30VEI30IVP75DIVA1))
VES30VEI30IVP90DIVA1df = as.data.frame(do.call(rbind,VES30VEI30IVP90DIVA1))
VES30VEI30IVP100DIVA1df = as.data.frame(do.call(rbind,VES30VEI30IVP100DIVA1))

VES30VEI60IVP0DIVA1df = as.data.frame(do.call(rbind,VES30VEI60IVP0DIVA1))
VES30VEI60IVP10DIVA1df = as.data.frame(do.call(rbind,VES30VEI60IVP10DIVA1))
VES30VEI60IVP25DIVA1df = as.data.frame(do.call(rbind,VES30VEI60IVP25DIVA1))
VES30VEI60IVP50DIVA1df = as.data.frame(do.call(rbind,VES30VEI60IVP50DIVA1))
VES30VEI60IVP75DIVA1df = as.data.frame(do.call(rbind,VES30VEI60IVP75DIVA1))
VES30VEI60IVP90DIVA1df = as.data.frame(do.call(rbind,VES30VEI60IVP90DIVA1))
VES30VEI60IVP100DIVA1df = as.data.frame(do.call(rbind,VES30VEI60IVP100DIVA1))

VES30VEI90IVP0DIVA1df = as.data.frame(do.call(rbind,VES30VEI90IVP0DIVA1))
VES30VEI90IVP10DIVA1df = as.data.frame(do.call(rbind,VES30VEI90IVP10DIVA1))
VES30VEI90IVP25DIVA1df = as.data.frame(do.call(rbind,VES30VEI90IVP25DIVA1))
VES30VEI90IVP50DIVA1df = as.data.frame(do.call(rbind,VES30VEI90IVP50DIVA1))
VES30VEI90IVP75DIVA1df = as.data.frame(do.call(rbind,VES30VEI90IVP75DIVA1))
VES30VEI90IVP90DIVA1df = as.data.frame(do.call(rbind,VES30VEI90IVP90DIVA1))
VES30VEI90IVP100DIVA1df = as.data.frame(do.call(rbind,VES30VEI90IVP100DIVA1))

save.image('VEISweepSORP1.RData')

VES60VEI0IVP0DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.0,180,0,0.6,0.0,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI0IVP10DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.1,180,0,0.6,0.0,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI0IVP25DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.25,180,0,0.6,0.0,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI0IVP50DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.50,180,0,0.6,0.0,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI0IVP75DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.75,180,0,0.6,0.0,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI0IVP90DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.90,180,0,0.6,0.0,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI0IVP100DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,1.0,180,0,0.6,0.0,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 

VES60VEI30IVP0DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.0,180,0,0.6,0.3,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI30IVP10DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.1,180,0,0.6,0.3,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI30IVP25DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.25,180,0,0.6,0.3,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI30IVP50DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.50,180,0,0.6,0.3,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI30IVP75DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.75,180,0,0.6,0.3,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI30IVP90DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.90,180,0,0.6,0.3,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI30IVP100DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,1.0,180,0,0.6,0.3,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 

VES60VEI60IVP0DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.0,180,0,0.6,0.6,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI60IVP10DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.1,180,0,0.6,0.6,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI60IVP25DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.25,180,0,0.6,0.6,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI60IVP50DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.50,180,0,0.6,0.6,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI60IVP75DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.75,180,0,0.6,0.6,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI60IVP90DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.90,180,0,0.6,0.6,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI60IVP100DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,1.0,180,0,0.6,0.6,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 

VES60VEI90IVP0DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.0,180,0,0.6,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI90IVP10DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.1,180,0,0.6,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI90IVP25DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.25,180,0,0.6,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI90IVP50DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.50,180,0,0.6,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI90IVP75DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.75,180,0,0.6,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI90IVP90DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.90,180,0,0.6,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES60VEI90IVP100DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,1.0,180,0,0.6,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 


VES60VEI0IVP0DIVA1df = as.data.frame(do.call(rbind,VES60VEI0IVP0DIVA1))
VES60VEI0IVP10DIVA1df = as.data.frame(do.call(rbind,VES60VEI0IVP10DIVA1))
VES60VEI0IVP25DIVA1df = as.data.frame(do.call(rbind,VES60VEI0IVP25DIVA1))
VES60VEI0IVP50DIVA1df = as.data.frame(do.call(rbind,VES60VEI0IVP50DIVA1))
VES60VEI0IVP75DIVA1df = as.data.frame(do.call(rbind,VES60VEI0IVP75DIVA1))
VES60VEI0IVP90DIVA1df = as.data.frame(do.call(rbind,VES60VEI0IVP90DIVA1))
VES60VEI0IVP100DIVA1df = as.data.frame(do.call(rbind,VES60VEI0IVP100DIVA1))

VES60VEI30IVP0DIVA1df = as.data.frame(do.call(rbind,VES60VEI30IVP0DIVA1))
VES60VEI30IVP10DIVA1df = as.data.frame(do.call(rbind,VES60VEI30IVP10DIVA1))
VES60VEI30IVP25DIVA1df = as.data.frame(do.call(rbind,VES60VEI30IVP25DIVA1))
VES60VEI30IVP50DIVA1df = as.data.frame(do.call(rbind,VES60VEI30IVP50DIVA1))
VES60VEI30IVP75DIVA1df = as.data.frame(do.call(rbind,VES60VEI30IVP75DIVA1))
VES60VEI30IVP90DIVA1df = as.data.frame(do.call(rbind,VES60VEI30IVP90DIVA1))
VES60VEI30IVP100DIVA1df = as.data.frame(do.call(rbind,VES60VEI30IVP100DIVA1))

VES60VEI60IVP0DIVA1df = as.data.frame(do.call(rbind,VES60VEI60IVP0DIVA1))
VES60VEI60IVP10DIVA1df = as.data.frame(do.call(rbind,VES60VEI60IVP10DIVA1))
VES60VEI60IVP25DIVA1df = as.data.frame(do.call(rbind,VES60VEI60IVP25DIVA1))
VES60VEI60IVP50DIVA1df = as.data.frame(do.call(rbind,VES60VEI60IVP50DIVA1))
VES60VEI60IVP75DIVA1df = as.data.frame(do.call(rbind,VES60VEI60IVP75DIVA1))
VES60VEI60IVP90DIVA1df = as.data.frame(do.call(rbind,VES60VEI60IVP90DIVA1))
VES60VEI60IVP100DIVA1df = as.data.frame(do.call(rbind,VES60VEI60IVP100DIVA1))

VES60VEI90IVP0DIVA1df = as.data.frame(do.call(rbind,VES60VEI90IVP0DIVA1))
VES60VEI90IVP10DIVA1df = as.data.frame(do.call(rbind,VES60VEI90IVP10DIVA1))
VES60VEI90IVP25DIVA1df = as.data.frame(do.call(rbind,VES60VEI90IVP25DIVA1))
VES60VEI90IVP50DIVA1df = as.data.frame(do.call(rbind,VES60VEI90IVP50DIVA1))
VES60VEI90IVP75DIVA1df = as.data.frame(do.call(rbind,VES60VEI90IVP75DIVA1))
VES60VEI90IVP90DIVA1df = as.data.frame(do.call(rbind,VES60VEI90IVP90DIVA1))
VES60VEI90IVP100DIVA1df = as.data.frame(do.call(rbind,VES60VEI90IVP100DIVA1))

save.image('VEISweepSORP2.RData')

VES90VEI0IVP0DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.0,180,0,0.9,0.0,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI0IVP10DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.1,180,0,0.9,0.0,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI0IVP25DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.25,180,0,0.9,0.0,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI0IVP50DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.50,180,0,0.9,0.0,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI0IVP75DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.75,180,0,0.9,0.0,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI0IVP90DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.90,180,0,0.9,0.0,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI0IVP100DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,1.0,180,0,0.9,0.0,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 

VES90VEI30IVP0DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.0,180,0,0.9,0.3,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI30IVP10DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.1,180,0,0.9,0.3,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI30IVP25DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.25,180,0,0.9,0.3,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI30IVP50DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.50,180,0,0.9,0.3,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI30IVP75DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.75,180,0,0.9,0.3,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI30IVP90DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.90,180,0,0.9,0.3,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI30IVP100DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,1.0,180,0,0.9,0.3,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 

VES90VEI60IVP0DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.0,180,0,0.9,0.6,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI60IVP10DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.1,180,0,0.9,0.6,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI60IVP25DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.25,180,0,0.9,0.6,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI60IVP50DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.50,180,0,0.9,0.6,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI60IVP75DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.75,180,0,0.9,0.6,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI60IVP90DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.90,180,0,0.9,0.6,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI60IVP100DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,1.0,180,0,0.9,0.6,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 

VES90VEI90IVP0DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.0,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI90IVP10DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.1,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI90IVP25DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.25,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI90IVP50DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.50,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI90IVP75DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.75,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI90IVP90DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.90,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
VES90VEI90IVP100DIVA1 <- mclapply(1:reps,function(x){single_runIndivSORI(3,1.0,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 


VES90VEI0IVP0DIVA1df = as.data.frame(do.call(rbind,VES90VEI0IVP0DIVA1))
VES90VEI0IVP10DIVA1df = as.data.frame(do.call(rbind,VES90VEI0IVP10DIVA1))
VES90VEI0IVP25DIVA1df = as.data.frame(do.call(rbind,VES90VEI0IVP25DIVA1))
VES90VEI0IVP50DIVA1df = as.data.frame(do.call(rbind,VES90VEI0IVP50DIVA1))
VES90VEI0IVP75DIVA1df = as.data.frame(do.call(rbind,VES90VEI0IVP75DIVA1))
VES90VEI0IVP90DIVA1df = as.data.frame(do.call(rbind,VES90VEI0IVP90DIVA1))
VES90VEI0IVP100DIVA1df = as.data.frame(do.call(rbind,VES90VEI0IVP100DIVA1))

VES90VEI30IVP0DIVA1df = as.data.frame(do.call(rbind,VES90VEI30IVP0DIVA1))
VES90VEI30IVP10DIVA1df = as.data.frame(do.call(rbind,VES90VEI30IVP10DIVA1))
VES90VEI30IVP25DIVA1df = as.data.frame(do.call(rbind,VES90VEI30IVP25DIVA1))
VES90VEI30IVP50DIVA1df = as.data.frame(do.call(rbind,VES90VEI30IVP50DIVA1))
VES90VEI30IVP75DIVA1df = as.data.frame(do.call(rbind,VES90VEI30IVP75DIVA1))
VES90VEI30IVP90DIVA1df = as.data.frame(do.call(rbind,VES90VEI30IVP90DIVA1))
VES90VEI30IVP100DIVA1df = as.data.frame(do.call(rbind,VES90VEI30IVP100DIVA1))

VES90VEI60IVP0DIVA1df = as.data.frame(do.call(rbind,VES90VEI60IVP0DIVA1))
VES90VEI60IVP10DIVA1df = as.data.frame(do.call(rbind,VES90VEI60IVP10DIVA1))
VES90VEI60IVP25DIVA1df = as.data.frame(do.call(rbind,VES90VEI60IVP25DIVA1))
VES90VEI60IVP50DIVA1df = as.data.frame(do.call(rbind,VES90VEI60IVP50DIVA1))
VES90VEI60IVP75DIVA1df = as.data.frame(do.call(rbind,VES90VEI60IVP75DIVA1))
VES90VEI60IVP90DIVA1df = as.data.frame(do.call(rbind,VES90VEI60IVP90DIVA1))
VES90VEI60IVP100DIVA1df = as.data.frame(do.call(rbind,VES90VEI60IVP100DIVA1))

VES90VEI90IVP0DIVA1df = as.data.frame(do.call(rbind,VES90VEI90IVP0DIVA1))
VES90VEI90IVP10DIVA1df = as.data.frame(do.call(rbind,VES90VEI90IVP10DIVA1))
VES90VEI90IVP25DIVA1df = as.data.frame(do.call(rbind,VES90VEI90IVP25DIVA1))
VES90VEI90IVP50DIVA1df = as.data.frame(do.call(rbind,VES90VEI90IVP50DIVA1))
VES90VEI90IVP75DIVA1df = as.data.frame(do.call(rbind,VES90VEI90IVP75DIVA1))
VES90VEI90IVP90DIVA1df = as.data.frame(do.call(rbind,VES90VEI90IVP90DIVA1))
VES90VEI90IVP100DIVA1df = as.data.frame(do.call(rbind,VES90VEI90IVP100DIVA1))

save.image('VEISweepSOR.RData')
