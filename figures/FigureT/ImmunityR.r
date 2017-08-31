require(ggplot2)
require(grid)

load('../../simulations/Titration/SORI/ImmunityRetain.RData')

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

uni_imm = unique(SimsDf$immune)

split_p = 0.75

# Vaccinated herds DIVA 1 / VL
# Direct Efficacy

baselinec<-SimsDf$VL[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[1] & SimsDf$treatment==FALSE]
baselinep<-SimsDf$AtRisk[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[1] & SimsDf$treatment==FALSE ]

controlc<-SimsDf$VL[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[1] & SimsDf$treatment==TRUE]
controlp<-SimsDf$AtRisk[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[1] & SimsDf$treatment==TRUE ]

cases<-SimsDf$VL[SimsDf$Vaccinated==1 & SimsDf$immune == uni_imm[1] & SimsDf$treatment==TRUE]
pop<-SimsDf$AtRisk[SimsDf$Vaccinated==1 & SimsDf$immune == uni_imm[1] & SimsDf$treatment==TRUE]

df=rbind(data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,Im=uni_imm[1],measure='Direct(WH)',Herds=50,VP=50))
df=rbind(df,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,Im=uni_imm[1],measure='Direct(BH)',Herds=50,VP=50))

dfE=rbind(data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,Im=uni_imm[1],measure='Direct(WH)',Herds=50,VP=50))
dfE=rbind(dfE,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,Im=uni_imm[1],measure='Direct(BH)',Herds=50,VP=50))

for(h in seq(60,300,10))
{

df=rbind(df,data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,Im=uni_imm[1],measure='Direct(WH)',Herds=h,VP=50))
df=rbind(df,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,Im=uni_imm[1],measure='Direct(BH)',Herds=h,VP=50))

# Effect

dfE=rbind(dfE,data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,Im=uni_imm[1],measure='Direct(WH)',Herds=h,VP=50))
dfE=rbind(dfE,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,Im=uni_imm[1],measure='Direct(BH)',Herds=h,VP=50))

}


# Vaccinated herds DIVA 1 / VL
# Indirect Efficacy



baselinec<-SimsDf$VL[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[1] & SimsDf$treatment==FALSE]
baselinep<-SimsDf$AtRisk[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[1] & SimsDf$treatment==FALSE ]

controlc<-SimsDf$VL[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[1] & SimsDf$treatment==TRUE]
controlp<-SimsDf$AtRisk[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[1] & SimsDf$treatment==TRUE ]

cases<-SimsDf$VL[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[1] & SimsDf$treatment==TRUE]
pop<-SimsDf$AtRisk[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[1] & SimsDf$treatment==TRUE ]

for(h in seq(50,300,10))
{

df=rbind(df,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,Im=uni_imm[1],measure='Indirect',Herds=h,VP=50))

# Effect

dfE=rbind(dfE,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,Im=uni_imm[1],measure='Indirect',Herds=h,VP=50))
}




# Vaccinated herds DIVA 1 / VL
# Total Efficacy


baselinec<-SimsDf$VL[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[1] & SimsDf$treatment==FALSE]
baselinep<-SimsDf$AtRisk[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[1] & SimsDf$treatment==FALSE ]

controlc<-SimsDf$VL[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[1] & SimsDf$treatment==TRUE]
controlp<-SimsDf$AtRisk[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[1] & SimsDf$treatment==TRUE ]

cases<-SimsDf$VL[SimsDf$immune == uni_imm[1] & SimsDf$treatment==TRUE]
pop<-SimsDf$AtRisk[SimsDf$immune == uni_imm[1] & SimsDf$treatment==TRUE ]

for(h in seq(50,300,10))
{
df=rbind(df,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,Im=uni_imm[1],measure='Total',Herds=h,VP=50))

# Effect

dfE=rbind(dfE,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,Im=uni_imm[1],measure='Total',Herds=h,VP=50))

}

for(i in 2:3)
{

# Vaccinated herds DIVA 1 / VL
# Direct Efficacy

baselinec<-SimsDf$VL[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[i] & SimsDf$treatment==FALSE]
baselinep<-SimsDf$AtRisk[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[i] & SimsDf$treatment==FALSE ]

controlc<-SimsDf$VL[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[i] & SimsDf$treatment==TRUE]
controlp<-SimsDf$AtRisk[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[i] & SimsDf$treatment==TRUE ]

cases<-SimsDf$VL[SimsDf$Vaccinated==1 & SimsDf$immune == uni_imm[i] & SimsDf$treatment==TRUE]
pop<-SimsDf$AtRisk[SimsDf$Vaccinated==1 & SimsDf$immune == uni_imm[i] & SimsDf$treatment==TRUE]

df=rbind(df,data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,Im=uni_imm[i],measure='Direct(WH)',Herds=50,VP=50))
df=rbind(df,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,Im=uni_imm[i],measure='Direct(BH)',Herds=50,VP=50))

dfE=rbind(dfE,data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,Im=uni_imm[i],measure='Direct(WH)',Herds=50,VP=50))
dfE=rbind(dfE,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,Im=uni_imm[i],measure='Direct(BH)',Herds=50,VP=50))

for(h in seq(60,300,10))
{

df=rbind(df,data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,Im=uni_imm[i],measure='Direct(WH)',Herds=h,VP=50))
df=rbind(df,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,Im=uni_imm[i],measure='Direct(BH)',Herds=h,VP=50))

# Effect

dfE=rbind(dfE,data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,Im=uni_imm[i],measure='Direct(WH)',Herds=h,VP=50))
dfE=rbind(dfE,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,Im=uni_imm[i],measure='Direct(BH)',Herds=h,VP=50))

}





# Vaccinated herds DIVA 1 / VL
# Indirect Efficacy



baselinec<-SimsDf$VL[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[i] & SimsDf$treatment==FALSE]
baselinep<-SimsDf$AtRisk[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[i] & SimsDf$treatment==FALSE ]

controlc<-SimsDf$VL[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[i] & SimsDf$treatment==TRUE]
controlp<-SimsDf$AtRisk[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[i] & SimsDf$treatment==TRUE ]

cases<-SimsDf$VL[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[i] & SimsDf$treatment==TRUE]
pop<-SimsDf$AtRisk[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[i] & SimsDf$treatment==TRUE ]

for(h in seq(50,300,10))
{

df=rbind(df,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,Im=uni_imm[i],measure='Indirect',Herds=h,VP=50))

# Effect

dfE=rbind(dfE,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,Im=uni_imm[i],measure='Indirect',Herds=h,VP=50))
}




# Vaccinated herds DIVA 1 / VL
# Total Efficacy


baselinec<-SimsDf$VL[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[i] & SimsDf$treatment==FALSE]
baselinep<-SimsDf$AtRisk[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[i] & SimsDf$treatment==FALSE ]

controlc<-SimsDf$VL[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[i] & SimsDf$treatment==TRUE]
controlp<-SimsDf$AtRisk[SimsDf$Vaccinated==0 & SimsDf$immune == uni_imm[i] & SimsDf$treatment==TRUE ]

cases<-SimsDf$VL[SimsDf$immune == uni_imm[i] & SimsDf$treatment==TRUE]
pop<-SimsDf$AtRisk[SimsDf$immune == uni_imm[i] & SimsDf$treatment==TRUE ]

for(h in seq(50,300,10))
{
df=rbind(df,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,Im=uni_imm[i],measure='Total',Herds=h,VP=50))

# Effect

dfE=rbind(dfE,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,Im=uni_imm[i],measure='Total',Herds=h,VP=50))

}



}






p4=(ggplot(dfE,aes(x=Herds,y=power,col=measure)) + geom_line()  +  scale_x_continuous('Herds',limits=c(50,300)) + scale_y_continuous('Efficacy',limits=c(-0.2,1)))
p3=(ggplot(df,aes(x=Herds,y=power,col=measure)) + geom_line()  +  scale_x_continuous('Herds',limits=c(50,300)) + scale_y_continuous('Power'))

















