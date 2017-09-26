
control <- ICDIVA1df$Breaklength[ICDIVA1df$Breaklength>0]>60
casesW <- VES30VEI0IVP100DIVA1df$Breaklength[VES30VEI0IVP100DIVA1df$Breaklength>0]>60
casesM <- VES30VEI0IVP50DIVA1df$Breaklength[VES30VEI0IVP50DIVA1df$Breaklength>0]>60

targus = SweepByHerdP(casesW,casesM,control,50/2)

df30=rbind(data.frame(median=targus[1,1],lower=targus[1,2],upper=targus[1,3],power=targus[1,4],VEI=0.0,coverage='100%',Herds=50,VP=100))
df30=rbind(df30,data.frame(median=targus[2,1],lower=targus[2,2],upper=targus[2,3],power=targus[2,4],VEI=0.0,coverage='50%',Herds=50,VP=100))

for(h in seq(50,2000,50))
{

targus = SweepByHerdP(casesW,casesM,control,h/2)

df30=rbind(df30,data.frame(median=targus[1,1],lower=targus[1,2],upper=targus[1,3],power=targus[1,4],VEI=0.0,coverage='100%',Herds=h,VP=100))
df30=rbind(df30,data.frame(median=targus[2,1],lower=targus[2,2],upper=targus[2,3],power=targus[2,4],VEI=0.0,coverage='50%',Herds=h,VP=100))

}


control <- ICDIVA1df$Breaklength[ICDIVA1df$Breaklength>0]>60
casesW <- VES30VEI30IVP100DIVA1df$Breaklength[VES30VEI30IVP100DIVA1df$Breaklength>0]>60
casesM <- VES30VEI30IVP50DIVA1df$Breaklength[VES30VEI30IVP50DIVA1df$Breaklength>0]>60

targus = SweepByHerdP(casesW,casesM,control,50/2)

df30=rbind(df30,data.frame(median=targus[1,1],lower=targus[1,2],upper=targus[1,3],power=targus[1,4],VEI=30.0,coverage='100%',Herds=50,VP=100))
df30=rbind(df30,data.frame(median=targus[2,1],lower=targus[2,2],upper=targus[2,3],power=targus[2,4],VEI=30.0,coverage='50%',Herds=50,VP=100))

for(h in seq(50,2000,50))
{

targus = SweepByHerdP(casesW,casesM,control,h/2)

df30=rbind(df30,data.frame(median=targus[1,1],lower=targus[1,2],upper=targus[1,3],power=targus[1,4],VEI=30.0,coverage='100%',Herds=h,VP=100))
df30=rbind(df30,data.frame(median=targus[2,1],lower=targus[2,2],upper=targus[2,3],power=targus[2,4],VEI=30.0,coverage='50%',Herds=h,VP=100))

}

control <- ICDIVA1df$Breaklength[ICDIVA1df$Breaklength>0]>60
casesW <- VES30VEI60IVP100DIVA1df$Breaklength[VES30VEI60IVP100DIVA1df$Breaklength>0]>60
casesM <- VES30VEI60IVP50DIVA1df$Breaklength[VES30VEI60IVP50DIVA1df$Breaklength>0]>60

targus = SweepByHerdP(casesW,casesM,control,50/2)

df30=rbind(df30,data.frame(median=targus[1,1],lower=targus[1,2],upper=targus[1,3],power=targus[1,4],VEI=60.0,coverage='100%',Herds=50,VP=100))
df30=rbind(df30,data.frame(median=targus[2,1],lower=targus[2,2],upper=targus[2,3],power=targus[2,4],VEI=60.0,coverage='50%',Herds=50,VP=100))

for(h in seq(50,2000,50))
{

targus = SweepByHerdP(casesW,casesM,control,h/2)

df30=rbind(df30,data.frame(median=targus[1,1],lower=targus[1,2],upper=targus[1,3],power=targus[1,4],VEI=60.0,coverage='100%',Herds=h,VP=100))
df30=rbind(df30,data.frame(median=targus[2,1],lower=targus[2,2],upper=targus[2,3],power=targus[2,4],VEI=60.0,coverage='50%',Herds=h,VP=100))

}

control <- ICDIVA1df$Breaklength[ICDIVA1df$Breaklength>0]>60
casesW <- VES30VEI90IVP100DIVA1df$Breaklength[VES30VEI90IVP100DIVA1df$Breaklength>0]>60
casesM <- VES30VEI90IVP50DIVA1df$Breaklength[VES30VEI90IVP50DIVA1df$Breaklength>0]>60

targus = SweepByHerdP(casesW,casesM,control,50/2)

df30=rbind(df30,data.frame(median=targus[1,1],lower=targus[1,2],upper=targus[1,3],power=targus[1,4],VEI=90.0,coverage='100%',Herds=50,VP=100))
df30=rbind(df30,data.frame(median=targus[2,1],lower=targus[2,2],upper=targus[2,3],power=targus[2,4],VEI=90.0,coverage='50%',Herds=50,VP=100))

for(h in seq(50,2000,50))
{

targus = SweepByHerdP(casesW,casesM,control,h/2)

df30=rbind(df30,data.frame(median=targus[1,1],lower=targus[1,2],upper=targus[1,3],power=targus[1,4],VEI=90.0,coverage='100%',Herds=h,VP=100))
df30=rbind(df30,data.frame(median=targus[2,1],lower=targus[2,2],upper=targus[2,3],power=targus[2,4],VEI=90.0,coverage='50%',Herds=h,VP=100))

}

df30$VEI=as.factor(df30$VEI)


