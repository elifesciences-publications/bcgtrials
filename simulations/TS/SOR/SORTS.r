# Load ABC posterior distributions
samples <- read.csv('../../../parameters/samplesSOR.csv')[,2:12]
weight  <- scan('../../../parameters/weightSOR.dat')

reps <- 5

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

single_runSORTS <- function(vacc_target,retain_reactors,vacc_effS,vacc_effI,DIVA,herd_no,Immunity=1.0/1.650)
{
trial_duration = 50
batch_interval = 180

#params <- samples[which(rmultinom(1,prob=weight,size=1)==1),]
params <- apply(resample(samples,weight),2,median)


ranx <- round(runif(1,0,.Machine$integer.max))
filename <- make_unique_filename(ranx)
filenameFu <- make_unique_filenameFu(ranx)
print(filename)

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
#slaughter <- try(read.csv(sprintf('output%d.csvFuReactorsSH0.dat',ranx),sep=' ',header=FALSE),silent=TRUE)
#Fu <- read.csv(filenameFu,sep=',')
#indiv <- try(read.csv(sprintf('output%d.csvFuIndividual0.dat',ranx)),silent=TRUE)

system(sprintf('rm %s',filename))
system(sprintf('rm %s*',filenameFu))

return(disco)
}

DIVA1 = c(73.3/100.0,1-99.85/100.0)

baseline<-single_runSORTS(0.0,0,0.0,0.0,DIVA1,4657-1)
treatment<-single_runSORTS(1.0,0,0.75,0.75,DIVA1,4657-1, 2.0)

plot(baseline$Time/360.0,baseline$R,pch=19)
points(treatment$Time/360.0,treatment$RU,pch=19,col='red')
points(treatment$Time/360.0,treatment$RV,pch=19,col='blue')

plot(cumsum(baseline$R),pch=19)
points(cumsum(treatment$RU),pch=19,col='red')
points(cumsum(treatment$RV),pch=19,col='blue')

par(mfrow=c(1,3))
plot(1-(treatment$RV/(treatment$OnHerd*0.5))/(treatment$RU/(treatment$OnHerd*0.5)),pch=19,ylim=c(0,1))
plot(1-(treatment$R/(treatment$OnHerd))/(baseline$R/baseline$OnHerd),pch=19,ylim=c(0,1))
plot(1-(treatment$RU/(treatment$OnHerd*0.5))/(baseline$R/baseline$OnHerd),pch=19,ylim=c(0,1))




