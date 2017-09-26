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

df30=rbind(data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,VE=0.0,measure='Direct(WH)',Herds=50,VP=50))
df30=rbind(df30,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,VE=0.0,measure='Direct(BH)',Herds=50,VP=50))

dfE30=rbind(data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,VE=0.0,measure='Direct(WH)',Herds=50,VP=50))
dfE30=rbind(dfE30,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,50,split_p),DIVA=1,VE=0.0,measure='Direct(BH)',Herds=50,VP=50))

for(h in seq(60,300,10))
{

df30=rbind(df30,data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Direct(WH)',Herds=h,VP=50))
df30=rbind(df30,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Direct(BH)',Herds=h,VP=50))

# Effect

dfE30=rbind(dfE30,data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Direct(WH)',Herds=h,VP=50))
dfE30=rbind(dfE30,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Direct(BH)',Herds=h,VP=50))

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

df30=rbind(df30,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Indirect',Herds=h,VP=50))

# Effect

dfE30=rbind(dfE30,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Indirect',Herds=h,VP=50))
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
df30=rbind(df30,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Total',Herds=h,VP=50))

# Effect

dfE30=rbind(dfE30,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=0.0,measure='Total',Herds=h,VP=50))

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

df30=rbind(df30,data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Direct(WH)',Herds=h,VP=50))
df30=rbind(df30,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Direct(BH)',Herds=h,VP=50))

# Effect

dfE30=rbind(dfE30,data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Direct(WH)',Herds=h,VP=50))
dfE30=rbind(dfE30,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Direct(BH)',Herds=h,VP=50))

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

df30=rbind(df30,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Indirect',Herds=h,VP=50))

# Effect

dfE30=rbind(dfE30,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Indirect',Herds=h,VP=50))

}

#plot(df30$nherds,df30$power,pch=19,type='b',ylim=c(0,1))

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

df30=rbind(df30,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Total',Herds=h,VP=50))

# Effect

dfE30=rbind(dfE30,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=30.0,measure='Total',Herds=h,VP=50))

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

df30=rbind(df30,data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Direct(WH)',Herds=50,VP=50))
df30=rbind(df30,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Direct(BH)',Herds=50,VP=50))

# Effect

dfE30=rbind(dfE30,data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Direct(WH)',Herds=50,VP=50))
dfE30=rbind(dfE30,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Direct(BH)',Herds=50,VP=50))

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

df30=rbind(df30,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Indirect',Herds=h,VP=50))

# EFFECT

dfE30=rbind(dfE30,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Indirect',Herds=h,VP=50))

}

#plot(df30$nherds,df30$power,pch=19,type='b',ylim=c(0,1))

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

df30=rbind(df30,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Total',Herds=h,VP=50))

# EFFECT

dfE30=rbind(dfE30,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=60.0,measure='Total',Herds=h,VP=50))

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

df30=rbind(df30,data.frame(SweepByHerd_MixedD1(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Direct(WH)',Herds=h,VP=50))
df30=rbind(df30,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Direct(BH)',Herds=h,VP=50))

# EFFECT

dfE30=rbind(dfE30,data.frame(SweepByHerd_MixedD1E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Direct(WH)',Herds=h,VP=50))
dfE30=rbind(dfE30,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Direct(BH)',Herds=h,VP=50))
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

df30=rbind(df30,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Indirect',Herds=h,VP=50))

# EFFECT

dfE30=rbind(dfE30,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Indirect',Herds=h,VP=50))

}


#plot(df30$nherds,df30$power,pch=19,type='b',ylim=c(0,1))

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

df30=rbind(df30,data.frame(SweepByHerd_MixedD2(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Total',Herds=h,VP=50))

#EFFECT

dfE30=rbind(dfE30,data.frame(SweepByHerd_MixedD2E(cases,pop,baselinec,baselinep,controlc,controlp,h,split_p),DIVA=1,VE=90.0,measure='Total',Herds=h,VP=50))

}

df30$VE = as.factor(df30$VE)
df30$VP = as.factor(df30$VP)
df30$measure = as.factor(df30$measure)

dfE30$VE = as.factor(dfE30$VE)
dfE30$VP = as.factor(dfE30$VP)
dfE30$measure = as.factor(dfE30$measure)







