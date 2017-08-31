require(ggplot2)
require(grid)

SweepByHerd <-function(cases,pop,controlc,controlp,nherds)
{

ARV = numeric(10000)

ARV = sapply(rep(nherds,10000),function(x,cases,pop,controlc,controlp)
{
k<-sample(1:length(cases),x,replace=TRUE)
j<-sample(1:length(controlc),x,replace=TRUE)
return(cbind(sum(cases[k])/sum(pop[k]) / (sum(controlc[j])/sum(controlp[j])),sqrt((1/sum(cases[k]))-(1/sum(pop[k]))+(1/sum(controlc[j]))-(1/sum(controlp[j]))),sum(pop[k]),sum(controlp[j])))
},cases=cases,pop=pop,controlc=controlc,controlp=controlp)
ARV=t(ARV)

zcr = qnorm(p = 1-0.025, mean = 0, sd = 1)
#power = mean((1 - pnorm(q = zcr, mean = abs(log(ARV[,1])/ARV[,2]), sd = 1)))
power = mean((log(ARV[,1])/ARV[,2] < (-zcr)))
#power=1-mean(pnorm(-abs(log(ARV[,1])/ARV[,2]))>0.05)
return(data.frame(nherds,power))
}
SweepByHerd_MixedD1 <-function(cases,pop,baselinec,baselinep,controlc,controlp,nherds,pmix)
{

ARV = numeric(10000)

ARV = sapply(rep(nherds,10000),function(x,cases,pop,baselinec,baselinep,controlc,controlp,pmix)
{

y=as.integer(pmix*x)
z=x-y


m<-sample(1:length(cases),y,replace=TRUE)
n<-sample(1:length(controlc),y,replace=TRUE)

oot = sum(cases[m])/sum(pop[m]) / (sum(controlc[n])/sum(controlp[n]))
poot = sqrt((1/sum(cases[m]))-(1/sum(pop[m]))+(1/sum(controlc[n]))-(1/sum(controlp[n])))


return(cbind(oot,poot))
},cases=cases,pop=pop,baselinec=baselinec,baselinep=baselinep,controlc=controlc,controlp=controlp,pmix=pmix)
ARV=t(ARV)

zcr = qnorm(p = 1-0.025, mean = 0, sd = 1)
#power = mean((1 - pnorm(q = zcr, mean = abs(log(ARV[,1])/ARV[,2]), sd = 1)))
power = mean((log(ARV[,1])/ARV[,2] < (-zcr)))
#power=1-mean(pnorm(-abs(log(ARV[,1])/ARV[,2]))>0.05)
return(data.frame(nherds,power))
}



SweepByHerd_MixedD2 <-function(cases,pop,baselinec,baselinep,controlc,controlp,nherds,pmix)
{

ARV = numeric(10000)

ARV = sapply(rep(nherds,10000),function(x,cases,pop,baselinec,baselinep,controlc,controlp,pmix)
{

y=as.integer(pmix*x)
z=x-y


k<-sample(1:length(baselinec),z,replace=TRUE)
m<-sample(1:length(cases),y,replace=TRUE)
n<-sample(1:length(controlc),y,replace=TRUE)

oot = sum(cases[m])/sum(pop[m]) / (sum(baselinec[k])/sum(baselinep[k]))
poot = sqrt((1/sum(cases[m]))-(1/sum(pop[m]))+(1/sum(baselinec[k]))-(1/sum(baselinep[k])))

return(cbind(oot,poot))
},cases=cases,pop=pop,baselinec=baselinec,baselinep=baselinep,controlc=controlc,controlp=controlp,pmix=pmix)
ARV=t(ARV)

zcr = qnorm(p = 1-0.025, mean = 0, sd = 1)
#power = mean((1 - pnorm(q = zcr, mean = abs(log(ARV[,1])/ARV[,2]), sd = 1)))
power = mean((log(ARV[,1])/ARV[,2] < (-zcr)))
#power=1-mean(pnorm(-abs(log(ARV[,1])/ARV[,2]))>0.05)
return(data.frame(nherds,power))

}

SweepByHerd_MixedD1E <-function(cases,pop,baselinec,baselinep,controlc,controlp,nherds,pmix)
{

ARV = numeric(10000)

ARV = sapply(rep(nherds,10000),function(x,cases,pop,baselinec,baselinep,controlc,controlp,pmix)
{

y=as.integer(pmix*x)
z=x-y


m<-sample(1:length(cases),y,replace=TRUE)
n<-sample(1:length(controlc),y,replace=TRUE)

oot = 1-sum(cases[m])/sum(pop[m]) / (sum(controlc[n])/sum(controlp[n]))

return(oot)
},cases=cases,pop=pop,baselinec=baselinec,baselinep=baselinep,controlc=controlc,controlp=controlp,pmix=pmix)
#return(ARV)
power=mean(ARV)
return(data.frame(nherds,power))
}



SweepByHerd_MixedD2E <-function(cases,pop,baselinec,baselinep,controlc,controlp,nherds,pmix)
{

ARV = numeric(10000)

ARV = sapply(rep(nherds,10000),function(x,cases,pop,baselinec,baselinep,controlc,controlp,pmix)
{

y=as.integer(pmix*x)
z=x-y


k<-sample(1:length(baselinec),z,replace=TRUE)
m<-sample(1:length(cases),y,replace=TRUE)
n<-sample(1:length(controlc),y,replace=TRUE)

oot = 1-sum(cases[m])/sum(pop[m]) / (sum(baselinec[k])/sum(baselinep[k]))

return(oot)
},cases=cases,pop=pop,baselinec=baselinec,baselinep=baselinep,controlc=controlc,controlp=controlp,pmix=pmix)
#return(ARV)
power=mean(ARV)
return(data.frame(nherds,power))
}

# Herd-by-herd controls

split_p = 0.75

# Vaccinated herds DIVA 1 / VL
# Direct Efficacy

baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES30VEI0IVP50DIVA1df$VL[VES30VEI0IVP50DIVA1df$Vaccinated==0]
controlp<-VES30VEI0IVP50DIVA1df$AtRisk[VES30VEI0IVP50DIVA1df$Vaccinated==0]

cases<-VES30VEI0IVP50DIVA1df$VL[VES30VEI0IVP50DIVA1df$Vaccinated==1]
pop<-VES30VEI0IVP50DIVA1df$AtRisk[VES30VEI0IVP50DIVA1df$Vaccinated==1]

df=rbind(data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,VE=0.0,measure='Direct(WH)',Herds=50,VP=50))
df=rbind(df,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,VE=0.0,measure='Direct(BH)',Herds=50,VP=50))

dfE=rbind(data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,VE=0.0,measure='Direct(WH)',Herds=50,VP=50))
dfE=rbind(dfE,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,VE=0.0,measure='Direct(BH)',Herds=50,VP=50))

for(h in seq(60,300,10))
{

df=rbind(df,data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Direct(WH)',Herds=h,VP=50))
df=rbind(df,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Direct(BH)',Herds=h,VP=50))

# Effect

dfE=rbind(dfE,data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Direct(WH)',Herds=h,VP=50))
dfE=rbind(dfE,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Direct(BH)',Herds=h,VP=50))

}

# Vaccinated herds DIVA 1 / VL
# Indirect Efficacy

baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES30VEI0IVP50DIVA1df$VL[VES30VEI0IVP50DIVA1df$Vaccinated==0]
controlp<-VES30VEI0IVP50DIVA1df$AtRisk[VES30VEI0IVP50DIVA1df$Vaccinated==0]

cases<-VES30VEI0IVP50DIVA1df$VL[VES30VEI0IVP50DIVA1df$Vaccinated==0]
pop<-VES30VEI0IVP50DIVA1df$AtRisk[VES30VEI0IVP50DIVA1df$Vaccinated==0]

for(h in seq(50,300,10))
{

df=rbind(df,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Indirect',Herds=h,VP=50))

# Effect

dfE=rbind(dfE,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Indirect',Herds=h,VP=50))
}


# Vaccinated herds DIVA 1 / VL
# Total Efficacy


baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES30VEI0IVP50DIVA1df$VL[VES30VEI0IVP50DIVA1df$Vaccinated==0]
controlp<-VES30VEI0IVP50DIVA1df$AtRisk[VES30VEI0IVP50DIVA1df$Vaccinated==0]

cases<-VES30VEI0IVP50DIVA1df$VL
pop<-VES30VEI0IVP50DIVA1df$AtRisk

for(h in seq(50,300,10))
{
df=rbind(df,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Total',Herds=h,VP=50))

# Effect

dfE=rbind(dfE,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Total',Herds=h,VP=50))

}

#plot(dfT$nherds,dfT$power,pch=19,type='b',ylim=c(0,1))

# Direct Efficacy

baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES30VEI30IVP50DIVA1df$VL[VES30VEI30IVP50DIVA1df$Vaccinated==0]
controlp<-VES30VEI30IVP50DIVA1df$AtRisk[VES30VEI30IVP50DIVA1df$Vaccinated==0]

cases<-VES30VEI30IVP50DIVA1df$VL[VES30VEI30IVP50DIVA1df$Vaccinated==1]
pop<-VES30VEI30IVP50DIVA1df$AtRisk[VES30VEI30IVP50DIVA1df$Vaccinated==1]

for(h in seq(50,300,10))
{

df=rbind(df,data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Direct(WH)',Herds=h,VP=50))
df=rbind(df,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Direct(BH)',Herds=h,VP=50))

# Effect

dfE=rbind(dfE,data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Direct(WH)',Herds=h,VP=50))
dfE=rbind(dfE,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Direct(BH)',Herds=h,VP=50))

}

# Vaccinated herds DIVA 1 / VL
# Indirect Efficacy

baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES30VEI30IVP50DIVA1df$VL[VES30VEI30IVP50DIVA1df$Vaccinated==0]
controlp<-VES30VEI30IVP50DIVA1df$AtRisk[VES30VEI30IVP50DIVA1df$Vaccinated==0]

cases<-VES30VEI30IVP50DIVA1df$VL[VES30VEI30IVP50DIVA1df$Vaccinated==0]
pop<-VES30VEI30IVP50DIVA1df$AtRisk[VES30VEI30IVP50DIVA1df$Vaccinated==0]

for(h in seq(50,300,10))
{

df=rbind(df,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Indirect',Herds=h,VP=50))

# Effect

dfE=rbind(dfE,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Indirect',Herds=h,VP=50))

}

#plot(df$nherds,df$power,pch=19,type='b',ylim=c(0,1))

# Vaccinated herds DIVA 1 / VL
# Total Efficacy

baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES30VEI30IVP50DIVA1df$VL[VES30VEI30IVP50DIVA1df$Vaccinated==0]
controlp<-VES30VEI30IVP50DIVA1df$AtRisk[VES30VEI30IVP50DIVA1df$Vaccinated==0]

cases<-VES30VEI30IVP50DIVA1df$VL
pop<-VES30VEI30IVP50DIVA1df$AtRisk

for(h in seq(50,300,10))
{

df=rbind(df,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Total',Herds=h,VP=50))

# Effect

dfE=rbind(dfE,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Total',Herds=h,VP=50))

}

# Direct Efficacy

baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES30VEI60IVP50DIVA1df$VL[VES30VEI60IVP50DIVA1df$Vaccinated==0]
controlp<-VES30VEI60IVP50DIVA1df$AtRisk[VES30VEI60IVP50DIVA1df$Vaccinated==0]

cases<-VES30VEI60IVP50DIVA1df$VL[VES30VEI60IVP50DIVA1df$Vaccinated==1]
pop<-VES30VEI60IVP50DIVA1df$AtRisk[VES30VEI60IVP50DIVA1df$Vaccinated==1]

for(h in seq(50,300,10))
{

df=rbind(df,data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Direct(WH)',Herds=50,VP=50))
df=rbind(df,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Direct(BH)',Herds=50,VP=50))

# Effect

dfE=rbind(dfE,data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Direct(WH)',Herds=50,VP=50))
dfE=rbind(dfE,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Direct(BH)',Herds=50,VP=50))

}

# Vaccinated herds DIVA 1 / VL
# Indirect Efficacy

baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES30VEI60IVP50DIVA1df$VL[VES30VEI60IVP50DIVA1df$Vaccinated==0]
controlp<-VES30VEI60IVP50DIVA1df$AtRisk[VES30VEI60IVP50DIVA1df$Vaccinated==0]

cases<-VES30VEI60IVP50DIVA1df$VL[VES30VEI60IVP50DIVA1df$Vaccinated==0]
pop<-VES30VEI60IVP50DIVA1df$AtRisk[VES30VEI60IVP50DIVA1df$Vaccinated==0]

for(h in seq(50,300,10))
{

df=rbind(df,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Indirect',Herds=h,VP=50))

# EFFECT

dfE=rbind(dfE,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Indirect',Herds=h,VP=50))

}

#plot(df$nherds,df$power,pch=19,type='b',ylim=c(0,1))

# Vaccinated herds DIVA 1 / VL
# Total Efficacy

baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES30VEI60IVP50DIVA1df$VL[VES30VEI60IVP50DIVA1df$Vaccinated==0]
controlp<-VES30VEI60IVP50DIVA1df$AtRisk[VES30VEI60IVP50DIVA1df$Vaccinated==0]

cases<-VES30VEI60IVP50DIVA1df$VL
pop<-VES30VEI60IVP50DIVA1df$AtRisk

for(h in seq(50,300,10))
{

df=rbind(df,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Total',Herds=h,VP=50))

# EFFECT

dfE=rbind(dfE,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Total',Herds=h,VP=50))

}

# Direct Efficacy

baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES30VEI90IVP50DIVA1df$VL[VES30VEI90IVP50DIVA1df$Vaccinated==0]
controlp<-VES30VEI90IVP50DIVA1df$AtRisk[VES30VEI90IVP50DIVA1df$Vaccinated==0]

cases<-VES30VEI90IVP50DIVA1df$VL[VES30VEI90IVP50DIVA1df$Vaccinated==1]
pop<-VES30VEI90IVP50DIVA1df$AtRisk[VES30VEI90IVP50DIVA1df$Vaccinated==1]


for(h in seq(50,300,10))
{

df=rbind(df,data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Direct(WH)',Herds=h,VP=50))
df=rbind(df,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Direct(BH)',Herds=h,VP=50))

# EFFECT

dfE=rbind(dfE,data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Direct(WH)',Herds=h,VP=50))
dfE=rbind(dfE,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Direct(BH)',Herds=h,VP=50))
}

# Vaccinated herds DIVA 1 / VL
# Indirect Efficacy

baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES30VEI90IVP50DIVA1df$VL[VES30VEI90IVP50DIVA1df$Vaccinated==0]
controlp<-VES30VEI90IVP50DIVA1df$AtRisk[VES30VEI90IVP50DIVA1df$Vaccinated==0]

cases<-VES30VEI90IVP50DIVA1df$VL[VES30VEI90IVP50DIVA1df$Vaccinated==0]
pop<-VES30VEI90IVP50DIVA1df$AtRisk[VES30VEI90IVP50DIVA1df$Vaccinated==0]


for(h in seq(50,300,10))
{

df=rbind(df,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Indirect',Herds=h,VP=50))

# EFFECT

dfE=rbind(dfE,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Indirect',Herds=h,VP=50))

}


#plot(df$nherds,df$power,pch=19,type='b',ylim=c(0,1))

# Vaccinated herds DIVA 1 / VL
# Total Efficacy

baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES30VEI90IVP50DIVA1df$VL[VES30VEI90IVP50DIVA1df$Vaccinated==0]
controlp<-VES30VEI90IVP50DIVA1df$AtRisk[VES30VEI90IVP50DIVA1df$Vaccinated==0]

cases<-VES30VEI90IVP50DIVA1df$VL
pop<-VES30VEI90IVP50DIVA1df$AtRisk

for(h in seq(50,300,10))
{

df=rbind(df,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Total',Herds=h,VP=50))

#EFFECT

dfE=rbind(dfE,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Total',Herds=h,VP=50))

}

df$VE = as.factor(df$VE)
df$VP = as.factor(df$VP)
df$measure = as.factor(df$measure)

dfE$VE = as.factor(dfE$VE)
dfE$VP = as.factor(dfE$VP)
dfE$measure = as.factor(dfE$measure)

p4=(ggplot(dfE,aes(x=Herds,y=power,col=measure,linetype=VE,shape=measure)) + geom_line()  +  scale_x_continuous('Herds',limits=c(51,300)) + scale_y_continuous('Efficacy',limits=c(-0.2,1))) + scale_colour_manual(values=cbPaletteB) + geom_point() + guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Measures'),shape=guide_legend(order=2,title='Measures'))

p3=(ggplot(df,aes(x=Herds,y=power,col=measure,linetype=VE,shape=measure)) + geom_line()  +  scale_x_continuous('Herds',limits=c(51,300)) + scale_y_continuous('Power')) + scale_colour_manual(values=cbPaletteB) + geom_point() + guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Measures'),shape=guide_legend(order=2,title='Measures'))

pdf('PowerCompareEfficacyVES30byVI.pdf',height=2.5,width=8.5,family='Helvetica',pointsize=12)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,2)))
vplayout <- function(x,y)
{viewport(layout.pos.row=x,layout.pos.col=y)}
print(p3,vp=vplayout(1,2))
print(p4,vp=vplayout(1,1))
dev.off()



