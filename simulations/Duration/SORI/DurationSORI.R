# Load ABC posterior distributions
samples <- read.csv('../../../parameters/samplesSORI.csv')[,2:13]
weight  <- scan('../../../parameters/weightSORI.dat')

reps <- 1000

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
params[5],params[12],1.0,1.0,params[5],params[5],params[12],
0.0,0.0,params[6],0.0,0.0,0.0,(1.0-vacc_effI)*(params[6]),params[7],
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

agg=data.frame(ranx=ranx,Herd=herd_no,Vaccinated=0,Reactor=sum(indiv$Reactor[indiv$Vaccinated==0]), VL=sum(indiv$VL[indiv$Vaccinated==0]),AtRisk=sum(indiv$Vaccinated==0))
agg=rbind(agg,data.frame(ranx=ranx,Herd=herd_no,Vaccinated=1,Reactor=sum(indiv$Reactor[indiv$Vaccinated==1]), VL=sum(indiv$VL[indiv$Vaccinated==1]),AtRisk=sum(indiv$Vaccinated==1)))

return(agg)
}

DIVA1 = c(73.3/100.0,1-99.85/100.0)
DIVA2 = c(0.25,0.0)
DIVA3 = c(1.0,1-99.85/100.0)
DIVA4 = c(1.0,0.0)

# 3 Years

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.0,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = cbind(as.data.frame(do.call(rbind,outlist)),prop=0,duration=3,retain=F)

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.50,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=50,duration=3,retain=F))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(3,1.0,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=100,duration=3,retain=F))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.0,180,1,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=0,duration=3,retain=T))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(3,0.50,180,1,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=50,duration=3,retain=T))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(3,1.0,180,1,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=100,duration=3,retain=T))

save.image('DurationSORIP1.RData')

# 6 Years

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(6,0.0,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=0,duration=6,retain=F))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(6,0.50,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=50,duration=6,retain=F))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(6,1.0,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=100,duration=6,retain=F))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(6,0.0,180,1,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=0,duration=6,retain=T))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(6,0.50,180,1,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=50,duration=6,retain=T))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(6,1.0,180,1,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=100,duration=6,retain=T))

save.image('DurationSORIP2.RData')

# 9 Years

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(9,0.0,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=0,duration=9,retain=F))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(9,0.50,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=50,duration=9,retain=F))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(9,1.0,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=100,duration=9,retain=F))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(9,0.0,180,1,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=0,duration=9,retain=T))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(9,0.50,180,1,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=50,duration=9,retain=T))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(9,1.0,180,1,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=100,duration=9,retain=T))

save.image('DurationSORIP3.RData')

# 12 Years

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(12,0.0,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=0,duration=12,retain=F))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(12,0.50,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=50,duration=12,retain=F))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(12,1.0,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=100,duration=12,retain=F))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(12,0.0,180,1,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=0,duration=12,retain=T))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(12,0.50,180,1,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=50,duration=12,retain=T))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(12,1.0,180,1,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=100,duration=12,retain=T))

save.image('DurationSORIP4.RData')

# 15 Years

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(15,0.0,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=0,duration=15,retain=F))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(15,0.50,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=50,duration=15,retain=F))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(15,1.0,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=100,duration=15,retain=F))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(15,0.0,180,1,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=0,duration=15,retain=T))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(15,0.50,180,1,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=50,duration=15,retain=T))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(15,1.0,180,1,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=100,duration=15,retain=T))

save.image('DurationSORIP5.RData')

# 19 Years

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(19,0.0,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=0,duration=19,retain=F))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(19,0.50,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=50,duration=19,retain=F))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(19,1.0,180,0,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=100,duration=19,retain=F))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(19,0.0,180,1,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=0,duration=19,retain=T))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(19,0.50,180,1,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=50,duration=19,retain=T))

outlist <- mclapply(1:reps,function(x){single_runIndivSORI(19,1.0,180,1,0.9,0.9,DIVA1)},mc.set.seed=TRUE,mc.preschedule=TRUE,mc.cores=24) 
out_df = rbind(out_df,cbind(as.data.frame(do.call(rbind,outlist)),prop=100,duration=19,retain=T))

save.image('DurationSORIP6.RData')




save.image('DurationSORI.RData')
