
SweepByHerdP1 <-function(casesW,casesM,control,nherds)
{

#ARV = numeric(10000)

ARV = sapply(rep(nherds,10000),function(x,casesW,casesM,control)
{
i<-sample(1:length(casesW),x,replace=TRUE)
k<-sample(1:length(casesM),x,replace=TRUE)
j<-sample(1:length(control),x,replace=TRUE)
return(cbind(mean(control[j]-casesW[i]),mean(control[j]-casesM[k])))
},casesW=casesW,casesM=casesM,control=control)

power=mean(ARV>0)
effectW=quantile(ARV[1,],c(0.025,0.5,0.975))
effectM=quantile(ARV[2,],c(0.025,0.5,0.975))
effect=rbind(effectW,effectM)
return(data.frame(median=rowMeans(ARV),lower=effect[,1],upper=effect[,3],power=(rowMeans(ARV>0))))
}


SweepByHerdP <-function(casesW,casesM,control,nherds)
{

#ARV = numeric(10000)

ARV = sapply(rep(nherds,10000),function(x,casesW,casesM,control)
{
i<-sample(1:length(casesW),x,replace=TRUE)
k<-sample(1:length(casesM),x,replace=TRUE)
j<-sample(1:length(control),x,replace=TRUE)

return(cbind(mean(control[j])-mean(casesW[i]),mean(c(control[j],casesW[i])),
mean(control[j])-mean(casesM[k]),mean(c(control[j],casesM[k]))))
},casesW=casesW,casesM=casesM,control=control)

ARV = t(ARV)
zcr = qnorm(p = 1-0.025, mean = 0, sd = 1)

ZW=ARV[,1]/sqrt(ARV[,2]*(1-ARV[,2])*((1/nherds)+(1/nherds)))
ZM=ARV[,3]/sqrt(ARV[,4]*(1-ARV[,4])*((1/nherds)+(1/nherds)))

power=cbind(mean((ZW > (zcr))),mean((ZM > (zcr))))

#power=cbind(mean(1-pnorm(-abs(ZW))),mean(1-pnorm(-abs(ZM))))

effectW=quantile(ARV[,1],c(0.025,0.5,0.975))
effectM=quantile(ARV[,3],c(0.025,0.5,0.975))
effect=rbind(effectW,effectM)
return(data.frame(median=colMeans(ARV[,c(1,3)]),lower=effect[,1],upper=effect[,3],power=t(power)))
}


control <- ICDIVA1df$Break_Recurr[ICDIVA1df$Breaklength>0]
casesW <- VES30VEI0IVP100DIVA1df$Break_Recurr[VES30VEI0IVP100DIVA1df$Breaklength>0]
casesM <- VES30VEI0IVP50DIVA1df$Break_Recurr[VES30VEI0IVP50DIVA1df$Breaklength>0]

control = control!=0
casesW = casesW!=0
casesM = casesM!=0

targus = SweepByHerdP(casesW,casesM,control,50/2)

df30=rbind(data.frame(median=targus[1,1],lower=targus[1,2],upper=targus[1,3],power=targus[1,4],VEI=0.0,coverage='100%',Herds=50,VP=100))
df30=rbind(df30,data.frame(median=targus[2,1],lower=targus[2,2],upper=targus[2,3],power=targus[2,4],VEI=0.0,coverage='50%',Herds=50,VP=100))

for(h in seq(50,2000,50))
{

targus = SweepByHerdP(casesW,casesM,control,h/2)

df30=rbind(df30,data.frame(median=targus[1,1],lower=targus[1,2],upper=targus[1,3],power=targus[1,4],VEI=0.0,coverage='100%',Herds=h,VP=100))
df30=rbind(df30,data.frame(median=targus[2,1],lower=targus[2,2],upper=targus[2,3],power=targus[2,4],VEI=0.0,coverage='50%',Herds=h,VP=100))

}

control <- ICDIVA1df$Break_Recurr[ICDIVA1df$Breaklength>0]
casesW <- VES30VEI30IVP100DIVA1df$Break_Recurr[VES30VEI30IVP100DIVA1df$Breaklength>0]
casesM <- VES30VEI30IVP50DIVA1df$Break_Recurr[VES30VEI30IVP50DIVA1df$Breaklength>0]

control = control!=0
casesW = casesW!=0
casesM = casesM!=0

targus = SweepByHerdP(casesW,casesM,control,50/2)

df30=rbind(df30,data.frame(median=targus[1,1],lower=targus[1,2],upper=targus[1,3],power=targus[1,4],VEI=30.0,coverage='100%',Herds=50,VP=100))
df30=rbind(df30,data.frame(median=targus[2,1],lower=targus[2,2],upper=targus[2,3],power=targus[2,4],VEI=30.0,coverage='50%',Herds=50,VP=100))

for(h in seq(50,2000,50))
{

targus = SweepByHerdP(casesW,casesM,control,h/2)

df30=rbind(df30,data.frame(median=targus[1,1],lower=targus[1,2],upper=targus[1,3],power=targus[1,4],VEI=30.0,coverage='100%',Herds=h,VP=100))
df30=rbind(df30,data.frame(median=targus[2,1],lower=targus[2,2],upper=targus[2,3],power=targus[2,4],VEI=30.0,coverage='50%',Herds=h,VP=100))

}

control <- ICDIVA1df$Break_Recurr[ICDIVA1df$Breaklength>0]
casesW <- VES30VEI60IVP100DIVA1df$Break_Recurr[VES30VEI60IVP100DIVA1df$Breaklength>0]
casesM <- VES30VEI60IVP50DIVA1df$Break_Recurr[VES30VEI60IVP50DIVA1df$Breaklength>0]

control = control!=0
casesW = casesW!=0
casesM = casesM!=0

targus = SweepByHerdP(casesW,casesM,control,50/2)

df30=rbind(df30,data.frame(median=targus[1,1],lower=targus[1,2],upper=targus[1,3],power=targus[1,4],VEI=60.0,coverage='100%',Herds=50,VP=100))
df30=rbind(df30,data.frame(median=targus[2,1],lower=targus[2,2],upper=targus[2,3],power=targus[2,4],VEI=60.0,coverage='50%',Herds=50,VP=100))

for(h in seq(50,2000,50))
{

targus = SweepByHerdP(casesW,casesM,control,h/2)

df30=rbind(df30,data.frame(median=targus[1,1],lower=targus[1,2],upper=targus[1,3],power=targus[1,4],VEI=60.0,coverage='100%',Herds=h,VP=100))
df30=rbind(df30,data.frame(median=targus[2,1],lower=targus[2,2],upper=targus[2,3],power=targus[2,4],VEI=60.0,coverage='50%',Herds=h,VP=100))

}

control <- ICDIVA1df$Break_Recurr[ICDIVA1df$Breaklength>0]
casesW <- VES30VEI90IVP100DIVA1df$Break_Recurr[VES30VEI90IVP100DIVA1df$Breaklength>0]
casesM <- VES30VEI90IVP50DIVA1df$Break_Recurr[VES30VEI90IVP50DIVA1df$Breaklength>0]

control = control!=0
casesW = casesW!=0
casesM = casesM!=0

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


