# Herd-by-herd controls

split_p = 0.75

# Vaccinated herds DIVA 1 / VL
# Direct Efficacy

baselinec<-na.omit(as.numeric(ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]))
baselinep<-na.omit(as.numeric(ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]))

controlc<-na.omit(as.numeric(VES90VEI0IVP50DIVA1df$VL[VES90VEI0IVP50DIVA1df$Vaccinated==0]))
controlp<-na.omit(as.numeric(VES90VEI0IVP50DIVA1df$AtRisk[VES90VEI0IVP50DIVA1df$Vaccinated==0]))

cases<-na.omit(as.numeric(VES90VEI0IVP50DIVA1df$VL[VES90VEI0IVP50DIVA1df$Vaccinated==1]))
pop<-na.omit(as.numeric(VES90VEI0IVP50DIVA1df$AtRisk[VES90VEI0IVP50DIVA1df$Vaccinated==1]))

df90=rbind(data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,VE=0.0,measure='Direct(WH)',Herds=50,VP=50))
df90=rbind(df90,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,VE=0.0,measure='Direct(BH)',Herds=50,VP=50))

dfE90=rbind(data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,VE=0.0,measure='Direct(WH)',Herds=50,VP=50))
dfE90=rbind(dfE90,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,VE=0.0,measure='Direct(BH)',Herds=50,VP=50))

for(h in seq(60,300,10))
{

df90=rbind(df90,data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Direct(WH)',Herds=h,VP=50))
df90=rbind(df90,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Direct(BH)',Herds=h,VP=50))

# Effect

dfE90=rbind(dfE90,data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Direct(WH)',Herds=h,VP=50))
dfE90=rbind(dfE90,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Direct(BH)',Herds=h,VP=50))

}

# Vaccinated herds DIVA 1 / VL
# Indirect Efficacy

baselinec<-na.omit(as.numeric(ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]))
baselinep<-na.omit(as.numeric(ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]))

controlc<-na.omit(as.numeric(VES90VEI0IVP50DIVA1df$VL[VES90VEI0IVP50DIVA1df$Vaccinated==0]))
controlp<-na.omit(as.numeric(VES90VEI0IVP50DIVA1df$AtRisk[VES90VEI0IVP50DIVA1df$Vaccinated==0]))

cases<-na.omit(as.numeric(VES90VEI0IVP50DIVA1df$VL[VES90VEI0IVP50DIVA1df$Vaccinated==0]))
pop<-na.omit(as.numeric(VES90VEI0IVP50DIVA1df$AtRisk[VES90VEI0IVP50DIVA1df$Vaccinated==0]))

for(h in seq(50,300,10))
{

df90=rbind(df90,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Indirect',Herds=h,VP=50))

# Effect

dfE90=rbind(dfE90,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Indirect',Herds=h,VP=50))
}


# Vaccinated herds DIVA 1 / VL
# Total Efficacy

baselinec<-na.omit(as.numeric(ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]))
baselinep<-na.omit(as.numeric(ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0])
)
controlc<-na.omit(as.numeric(VES90VEI0IVP50DIVA1df$VL[VES90VEI0IVP50DIVA1df$Vaccinated==0]))
controlp<-na.omit(as.numeric(VES90VEI0IVP50DIVA1df$AtRisk[VES90VEI0IVP50DIVA1df$Vaccinated==0]))

cases<-na.omit(as.numeric(VES90VEI0IVP50DIVA1df$VL))
pop<-na.omit(as.numeric(VES90VEI0IVP50DIVA1df$AtRisk))

for(h in seq(50,300,10))
{
df90=rbind(df90,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Total',Herds=h,VP=50))

# Effect

dfE90=rbind(dfE90,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Total',Herds=h,VP=50))

}

#plot(dfT$nherds,dfT$power,pch=19,type='b',ylim=c(0,1))

# Direct Efficacy

baselinec<-na.omit(as.numeric(ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]))
baselinep<-na.omit(as.numeric(ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]))

controlc<-na.omit(as.numeric(VES90VEI30IVP50DIVA1df$VL[VES90VEI30IVP50DIVA1df$Vaccinated==0]))
controlp<-na.omit(as.numeric(VES90VEI30IVP50DIVA1df$AtRisk[VES90VEI30IVP50DIVA1df$Vaccinated==0]))

cases<-na.omit(as.numeric(VES90VEI30IVP50DIVA1df$VL[VES90VEI30IVP50DIVA1df$Vaccinated==1]))
pop<-na.omit(as.numeric(VES90VEI30IVP50DIVA1df$AtRisk[VES90VEI30IVP50DIVA1df$Vaccinated==1]))

for(h in seq(50,300,10))
{

df90=rbind(df90,data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Direct(WH)',Herds=h,VP=50))
df90=rbind(df90,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Direct(BH)',Herds=h,VP=50))

# Effect

dfE90=rbind(dfE90,data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Direct(WH)',Herds=h,VP=50))
dfE90=rbind(dfE90,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Direct(BH)',Herds=h,VP=50))

}

# Vaccinated herds DIVA 1 / VL
# Indirect Efficacy

baselinec<-na.omit(as.numeric(ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]))
baselinep<-na.omit(as.numeric(ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]))

controlc<-na.omit(as.numeric(VES90VEI30IVP50DIVA1df$VL[VES90VEI30IVP50DIVA1df$Vaccinated==0]))
controlp<-na.omit(as.numeric(VES90VEI30IVP50DIVA1df$AtRisk[VES90VEI30IVP50DIVA1df$Vaccinated==0]))

cases<-na.omit(as.numeric(VES90VEI30IVP50DIVA1df$VL[VES90VEI30IVP50DIVA1df$Vaccinated==0]))
pop<-na.omit(as.numeric(VES90VEI30IVP50DIVA1df$AtRisk[VES90VEI30IVP50DIVA1df$Vaccinated==0]))

for(h in seq(50,300,10))
{

df90=rbind(df90,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Indirect',Herds=h,VP=50))

# Effect

dfE90=rbind(dfE90,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Indirect',Herds=h,VP=50))

}

#plot(df90$nherds,df90$power,pch=19,type='b',ylim=c(0,1))

# Vaccinated herds DIVA 1 / VL
# Total Efficacy

baselinec<-na.omit(as.numeric(ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]))
baselinep<-na.omit(as.numeric(ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]))

controlc<-na.omit(as.numeric(VES90VEI30IVP50DIVA1df$VL[VES90VEI30IVP50DIVA1df$Vaccinated==0]))
controlp<-na.omit(as.numeric(VES90VEI30IVP50DIVA1df$AtRisk[VES90VEI30IVP50DIVA1df$Vaccinated==0]))

cases<-na.omit(as.numeric(VES90VEI30IVP50DIVA1df$VL))
pop<-na.omit(as.numeric(VES90VEI30IVP50DIVA1df$AtRisk))

for(h in seq(50,300,10))
{

df90=rbind(df90,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Total',Herds=h,VP=50))

# Effect

dfE90=rbind(dfE90,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Total',Herds=h,VP=50))

}

# Direct Efficacy

baselinec<-na.omit(as.numeric(ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]))
baselinep<-na.omit(as.numeric(ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]))

controlc<-na.omit(as.numeric(VES90VEI60IVP50DIVA1df$VL[VES90VEI60IVP50DIVA1df$Vaccinated==0]))
controlp<-na.omit(as.numeric(VES90VEI60IVP50DIVA1df$AtRisk[VES90VEI60IVP50DIVA1df$Vaccinated==0]))

cases<-na.omit(as.numeric(VES90VEI60IVP50DIVA1df$VL[VES90VEI60IVP50DIVA1df$Vaccinated==1]))
pop<-na.omit(as.numeric(VES90VEI60IVP50DIVA1df$AtRisk[VES90VEI60IVP50DIVA1df$Vaccinated==1]))

for(h in seq(50,300,10))
{

df90=rbind(df90,data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Direct(WH)',Herds=50,VP=50))
df90=rbind(df90,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Direct(BH)',Herds=50,VP=50))

# Effect

dfE90=rbind(dfE90,data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Direct(WH)',Herds=50,VP=50))
dfE90=rbind(dfE90,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Direct(BH)',Herds=50,VP=50))

}

# Vaccinated herds DIVA 1 / VL
# Indirect Efficacy

baselinec<-na.omit(as.numeric(ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]))
baselinep<-na.omit(as.numeric(ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]))

controlc<-na.omit(as.numeric(VES90VEI60IVP50DIVA1df$VL[VES90VEI60IVP50DIVA1df$Vaccinated==0]))
controlp<-na.omit(as.numeric(VES90VEI60IVP50DIVA1df$AtRisk[VES90VEI60IVP50DIVA1df$Vaccinated==0]))

cases<-na.omit(as.numeric(VES90VEI60IVP50DIVA1df$VL[VES90VEI60IVP50DIVA1df$Vaccinated==0]))
pop<-na.omit(as.numeric(VES90VEI60IVP50DIVA1df$AtRisk[VES90VEI60IVP50DIVA1df$Vaccinated==0]))

for(h in seq(50,300,10))
{

df90=rbind(df90,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Indirect',Herds=h,VP=50))

# EFFECT

dfE90=rbind(dfE90,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Indirect',Herds=h,VP=50))

}

#plot(df90$nherds,df90$power,pch=19,type='b',ylim=c(0,1))

# Vaccinated herds DIVA 1 / VL
# Total Efficacy

baselinec<-na.omit(as.numeric(ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]))
baselinep<-na.omit(as.numeric(ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]))

controlc<-na.omit(as.numeric(VES90VEI60IVP50DIVA1df$VL[VES90VEI60IVP50DIVA1df$Vaccinated==0]))
controlp<-na.omit(as.numeric(VES90VEI60IVP50DIVA1df$AtRisk[VES90VEI60IVP50DIVA1df$Vaccinated==0]))

cases<-na.omit(as.numeric(VES90VEI60IVP50DIVA1df$VL))
pop<-na.omit(as.numeric(VES90VEI60IVP50DIVA1df$AtRisk))

for(h in seq(50,300,10))
{

df90=rbind(df90,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Total',Herds=h,VP=50))

# EFFECT

dfE90=rbind(dfE90,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Total',Herds=h,VP=50))

}

# Direct Efficacy

baselinec<-na.omit(as.numeric(ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]))
baselinep<-na.omit(as.numeric(ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]))

controlc<-na.omit(as.numeric(VES90VEI90IVP50DIVA1df$VL[VES90VEI90IVP50DIVA1df$Vaccinated==0]))
controlp<-na.omit(as.numeric(VES90VEI90IVP50DIVA1df$AtRisk[VES90VEI90IVP50DIVA1df$Vaccinated==0]))

cases<-na.omit(as.numeric(VES90VEI90IVP50DIVA1df$VL[VES90VEI90IVP50DIVA1df$Vaccinated==1]))
pop<-na.omit(as.numeric(VES90VEI90IVP50DIVA1df$AtRisk[VES90VEI90IVP50DIVA1df$Vaccinated==1]))


for(h in seq(50,300,10))
{

df90=rbind(df90,data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Direct(WH)',Herds=h,VP=50))
df90=rbind(df90,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Direct(BH)',Herds=h,VP=50))

# EFFECT

dfE90=rbind(dfE90,data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Direct(WH)',Herds=h,VP=50))
dfE90=rbind(dfE90,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Direct(BH)',Herds=h,VP=50))
}

# Vaccinated herds DIVA 1 / VL
# Indirect Efficacy

baselinec<-na.omit(as.numeric(ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]))
baselinep<-na.omit(as.numeric(ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]))

controlc<-na.omit(as.numeric(VES90VEI90IVP50DIVA1df$VL[VES90VEI90IVP50DIVA1df$Vaccinated==0]))
controlp<-na.omit(as.numeric(VES90VEI90IVP50DIVA1df$AtRisk[VES90VEI90IVP50DIVA1df$Vaccinated==0]))

cases<-na.omit(as.numeric(VES90VEI90IVP50DIVA1df$VL[VES90VEI90IVP50DIVA1df$Vaccinated==0]))
pop<-na.omit(as.numeric(VES90VEI90IVP50DIVA1df$AtRisk[VES90VEI90IVP50DIVA1df$Vaccinated==0]))


for(h in seq(50,300,10))
{

df90=rbind(df90,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Indirect',Herds=h,VP=50))

# EFFECT

dfE90=rbind(dfE90,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Indirect',Herds=h,VP=50))

}


#plot(df90$nherds,df90$power,pch=19,type='b',ylim=c(0,1))

# Vaccinated herds DIVA 1 / VL
# Total Efficacy

baselinec<-na.omit(as.numeric(ICDIVA1df$VL[ICDIVA1df$Vaccinated==0]))
baselinep<-na.omit(as.numeric(ICDIVA1df$AtRisk[ICDIVA1df$Vaccinated==0]))

controlc<-na.omit(as.numeric(VES90VEI90IVP50DIVA1df$VL[VES90VEI90IVP50DIVA1df$Vaccinated==0]))
controlp<-na.omit(as.numeric(VES90VEI90IVP50DIVA1df$AtRisk[VES90VEI90IVP50DIVA1df$Vaccinated==0]))

cases<-na.omit(as.numeric(VES90VEI90IVP50DIVA1df$VL))
pop<-na.omit(as.numeric(VES90VEI90IVP50DIVA1df$AtRisk))

for(h in seq(50,300,10))
{

df90=rbind(df90,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Total',Herds=h,VP=50))

#EFFECT

dfE90=rbind(dfE90,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Total',Herds=h,VP=50))

}

df90$VE = as.factor(df90$VE)
df90$VP = as.factor(df90$VP)
df90$measure = as.factor(df90$measure)

dfE90$VE = as.factor(dfE90$VE)
dfE90$VP = as.factor(dfE90$VP)
dfE90$measure = as.factor(dfE90$measure)


