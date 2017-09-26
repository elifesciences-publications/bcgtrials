

# Herd-by-herd controls

split_p = 0.75

# Vaccinated herds DIVA 1 / VL
# Direct Efficacy

baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES60VEI0IVP50DIVA1df$VL[VES60VEI0IVP50DIVA1df$Vaccinated==0]
controlp<-VES60VEI0IVP50DIVA1df$AtRisk[VES60VEI0IVP50DIVA1df$Vaccinated==0]

cases<-VES60VEI0IVP50DIVA1df$VL[VES60VEI0IVP50DIVA1df$Vaccinated==1]
pop<-VES60VEI0IVP50DIVA1df$AtRisk[VES60VEI0IVP50DIVA1df$Vaccinated==1]

df60=rbind(data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,VE=0.0,measure='Direct(WH)',Herds=50,VP=50))
df60=rbind(df60,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,VE=0.0,measure='Direct(BH)',Herds=50,VP=50))

dfE60=rbind(data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,VE=0.0,measure='Direct(WH)',Herds=50,VP=50))
dfE60=rbind(dfE60,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,VE=0.0,measure='Direct(BH)',Herds=50,VP=50))

for(h in seq(60,300,10))
{

df60=rbind(df60,data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Direct(WH)',Herds=h,VP=50))
df60=rbind(df60,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Direct(BH)',Herds=h,VP=50))

# Effect

dfE60=rbind(dfE60,data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Direct(WH)',Herds=h,VP=50))
dfE60=rbind(dfE60,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Direct(BH)',Herds=h,VP=50))

}

# Vaccinated herds DIVA 1 / VL
# Indirect Efficacy

baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES60VEI0IVP50DIVA1df$VL[VES60VEI0IVP50DIVA1df$Vaccinated==0]
controlp<-VES60VEI0IVP50DIVA1df$AtRisk[VES60VEI0IVP50DIVA1df$Vaccinated==0]

cases<-VES60VEI0IVP50DIVA1df$VL[VES60VEI0IVP50DIVA1df$Vaccinated==0]
pop<-VES60VEI0IVP50DIVA1df$AtRisk[VES60VEI0IVP50DIVA1df$Vaccinated==0]

for(h in seq(50,300,10))
{

df60=rbind(df60,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Indirect',Herds=h,VP=50))

# Effect

dfE60=rbind(dfE60,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Indirect',Herds=h,VP=50))
}


# Vaccinated herds DIVA 1 / VL
# Total Efficacy


baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES60VEI0IVP50DIVA1df$VL[VES60VEI0IVP50DIVA1df$Vaccinated==0]
controlp<-VES60VEI0IVP50DIVA1df$AtRisk[VES60VEI0IVP50DIVA1df$Vaccinated==0]

cases<-VES60VEI0IVP50DIVA1df$VL
pop<-VES60VEI0IVP50DIVA1df$AtRisk

for(h in seq(50,300,10))
{
df60=rbind(df60,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Total',Herds=h,VP=50))

# Effect

dfE60=rbind(dfE60,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Total',Herds=h,VP=50))

}

#plot(dfT$nherds,dfT$power,pch=19,type='b',ylim=c(0,1))

# Direct Efficacy

baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES60VEI30IVP50DIVA1df$VL[VES60VEI30IVP50DIVA1df$Vaccinated==0]
controlp<-VES60VEI30IVP50DIVA1df$AtRisk[VES60VEI30IVP50DIVA1df$Vaccinated==0]

cases<-VES60VEI30IVP50DIVA1df$VL[VES60VEI30IVP50DIVA1df$Vaccinated==1]
pop<-VES60VEI30IVP50DIVA1df$AtRisk[VES60VEI30IVP50DIVA1df$Vaccinated==1]

for(h in seq(50,300,10))
{

df60=rbind(df60,data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Direct(WH)',Herds=h,VP=50))
df60=rbind(df60,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Direct(BH)',Herds=h,VP=50))

# Effect

dfE60=rbind(dfE60,data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Direct(WH)',Herds=h,VP=50))
dfE60=rbind(dfE60,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Direct(BH)',Herds=h,VP=50))

}

# Vaccinated herds DIVA 1 / VL
# Indirect Efficacy

baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES60VEI30IVP50DIVA1df$VL[VES60VEI30IVP50DIVA1df$Vaccinated==0]
controlp<-VES60VEI30IVP50DIVA1df$AtRisk[VES60VEI30IVP50DIVA1df$Vaccinated==0]

cases<-VES60VEI30IVP50DIVA1df$VL[VES60VEI30IVP50DIVA1df$Vaccinated==0]
pop<-VES60VEI30IVP50DIVA1df$AtRisk[VES60VEI30IVP50DIVA1df$Vaccinated==0]

for(h in seq(50,300,10))
{

df60=rbind(df60,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Indirect',Herds=h,VP=50))

# Effect

dfE60=rbind(dfE60,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Indirect',Herds=h,VP=50))

}

#plot(df60$nherds,df60$power,pch=19,type='b',ylim=c(0,1))

# Vaccinated herds DIVA 1 / VL
# Total Efficacy

baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES60VEI30IVP50DIVA1df$VL[VES60VEI30IVP50DIVA1df$Vaccinated==0]
controlp<-VES60VEI30IVP50DIVA1df$AtRisk[VES60VEI30IVP50DIVA1df$Vaccinated==0]

cases<-VES60VEI30IVP50DIVA1df$VL
pop<-VES60VEI30IVP50DIVA1df$AtRisk

for(h in seq(50,300,10))
{

df60=rbind(df60,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Total',Herds=h,VP=50))

# Effect

dfE60=rbind(dfE60,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Total',Herds=h,VP=50))

}

# Direct Efficacy

baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES60VEI60IVP50DIVA1df$VL[VES60VEI60IVP50DIVA1df$Vaccinated==0]
controlp<-VES60VEI60IVP50DIVA1df$AtRisk[VES60VEI60IVP50DIVA1df$Vaccinated==0]

cases<-VES60VEI60IVP50DIVA1df$VL[VES60VEI60IVP50DIVA1df$Vaccinated==1]
pop<-VES60VEI60IVP50DIVA1df$AtRisk[VES60VEI60IVP50DIVA1df$Vaccinated==1]

for(h in seq(50,300,10))
{

df60=rbind(df60,data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Direct(WH)',Herds=50,VP=50))
df60=rbind(df60,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Direct(BH)',Herds=50,VP=50))

# Effect

dfE60=rbind(dfE60,data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Direct(WH)',Herds=50,VP=50))
dfE60=rbind(dfE60,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Direct(BH)',Herds=50,VP=50))

}

# Vaccinated herds DIVA 1 / VL
# Indirect Efficacy

baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES60VEI60IVP50DIVA1df$VL[VES60VEI60IVP50DIVA1df$Vaccinated==0]
controlp<-VES60VEI60IVP50DIVA1df$AtRisk[VES60VEI60IVP50DIVA1df$Vaccinated==0]

cases<-VES60VEI60IVP50DIVA1df$VL[VES60VEI60IVP50DIVA1df$Vaccinated==0]
pop<-VES60VEI60IVP50DIVA1df$AtRisk[VES60VEI60IVP50DIVA1df$Vaccinated==0]

for(h in seq(50,300,10))
{

df60=rbind(df60,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Indirect',Herds=h,VP=50))

# EFFECT

dfE60=rbind(dfE60,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Indirect',Herds=h,VP=50))

}

#plot(df60$nherds,df60$power,pch=19,type='b',ylim=c(0,1))

# Vaccinated herds DIVA 1 / VL
# Total Efficacy

baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES60VEI60IVP50DIVA1df$VL[VES60VEI60IVP50DIVA1df$Vaccinated==0]
controlp<-VES60VEI60IVP50DIVA1df$AtRisk[VES60VEI60IVP50DIVA1df$Vaccinated==0]

cases<-VES60VEI60IVP50DIVA1df$VL
pop<-VES60VEI60IVP50DIVA1df$AtRisk

for(h in seq(50,300,10))
{

df60=rbind(df60,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Total',Herds=h,VP=50))

# EFFECT

dfE60=rbind(dfE60,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Total',Herds=h,VP=50))

}

# Direct Efficacy

baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES60VEI90IVP50DIVA1df$VL[VES60VEI90IVP50DIVA1df$Vaccinated==0]
controlp<-VES60VEI90IVP50DIVA1df$AtRisk[VES60VEI90IVP50DIVA1df$Vaccinated==0]

cases<-VES60VEI90IVP50DIVA1df$VL[VES60VEI90IVP50DIVA1df$Vaccinated==1]
pop<-VES60VEI90IVP50DIVA1df$AtRisk[VES60VEI90IVP50DIVA1df$Vaccinated==1]


for(h in seq(50,300,10))
{

df60=rbind(df60,data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Direct(WH)',Herds=h,VP=50))
df60=rbind(df60,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Direct(BH)',Herds=h,VP=50))

# EFFECT

dfE60=rbind(dfE60,data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Direct(WH)',Herds=h,VP=50))
dfE60=rbind(dfE60,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Direct(BH)',Herds=h,VP=50))
}

# Vaccinated herds DIVA 1 / VL
# Indirect Efficacy

baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES60VEI90IVP50DIVA1df$VL[VES60VEI90IVP50DIVA1df$Vaccinated==0]
controlp<-VES60VEI90IVP50DIVA1df$AtRisk[VES60VEI90IVP50DIVA1df$Vaccinated==0]

cases<-VES60VEI90IVP50DIVA1df$VL[VES60VEI90IVP50DIVA1df$Vaccinated==0]
pop<-VES60VEI90IVP50DIVA1df$AtRisk[VES60VEI90IVP50DIVA1df$Vaccinated==0]


for(h in seq(50,300,10))
{

df60=rbind(df60,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Indirect',Herds=h,VP=50))

# EFFECT

dfE60=rbind(dfE60,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Indirect',Herds=h,VP=50))

}


#plot(df60$nherds,df60$power,pch=19,type='b',ylim=c(0,1))

# Vaccinated herds DIVA 1 / VL
# Total Efficacy

baselinec<-ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]
baselinep<-ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]

controlc<-VES60VEI90IVP50DIVA1df$VL[VES60VEI90IVP50DIVA1df$Vaccinated==0]
controlp<-VES60VEI90IVP50DIVA1df$AtRisk[VES60VEI90IVP50DIVA1df$Vaccinated==0]

cases<-VES60VEI90IVP50DIVA1df$VL
pop<-VES60VEI90IVP50DIVA1df$AtRisk

for(h in seq(50,300,10))
{

df60=rbind(df60,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Total',Herds=h,VP=50))

#EFFECT

dfE60=rbind(dfE60,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Total',Herds=h,VP=50))

}

df60$VE = as.factor(df60$VE)
df60$VP = as.factor(df60$VP)
df60$measure = as.factor(df60$measure)

dfE60$VE = as.factor(dfE60$VE)
dfE60$VP = as.factor(dfE60$VP)
dfE60$measure = as.factor(dfE60$measure)





