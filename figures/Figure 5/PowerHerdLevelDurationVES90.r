require(ggplot2)
require(grid)

load('VEISweepSORIP.RData')

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


control <- ICDIVA1df$Breaklength[ICDIVA1df$Breaklength>0]>60
casesW <- VES90VEI0IVP100DIVA1df$Breaklength[VES90VEI0IVP100DIVA1df$Breaklength>0]>60
casesM <- VES90VEI0IVP50DIVA1df$Breaklength[VES90VEI0IVP50DIVA1df$Breaklength>0]>60

targus = SweepByHerdP(casesW,casesM,control,50/2)

df=rbind(data.frame(median=targus[1,1],lower=targus[1,2],upper=targus[1,3],power=targus[1,4],VEI=0.0,coverage='100%',Herds=50,VP=100))
df=rbind(df,data.frame(median=targus[2,1],lower=targus[2,2],upper=targus[2,3],power=targus[2,4],VEI=0.0,coverage='50%',Herds=50,VP=100))

for(h in seq(50,2000,50))
{

targus = SweepByHerdP(casesW,casesM,control,h/2)

df=rbind(df,data.frame(median=targus[1,1],lower=targus[1,2],upper=targus[1,3],power=targus[1,4],VEI=0.0,coverage='100%',Herds=h,VP=100))
df=rbind(df,data.frame(median=targus[2,1],lower=targus[2,2],upper=targus[2,3],power=targus[2,4],VEI=0.0,coverage='50%',Herds=h,VP=100))

}


control <- ICDIVA1df$Breaklength[ICDIVA1df$Breaklength>0]>60
casesW <- VES90VEI30IVP100DIVA1df$Breaklength[VES90VEI30IVP100DIVA1df$Breaklength>0]>60
casesM <- VES90VEI30IVP50DIVA1df$Breaklength[VES90VEI30IVP50DIVA1df$Breaklength>0]>60

targus = SweepByHerdP(casesW,casesM,control,50/2)

df=rbind(df,data.frame(median=targus[1,1],lower=targus[1,2],upper=targus[1,3],power=targus[1,4],VEI=30.0,coverage='100%',Herds=50,VP=100))
df=rbind(df,data.frame(median=targus[2,1],lower=targus[2,2],upper=targus[2,3],power=targus[2,4],VEI=30.0,coverage='50%',Herds=50,VP=100))

for(h in seq(50,2000,50))
{

targus = SweepByHerdP(casesW,casesM,control,h/2)

df=rbind(df,data.frame(median=targus[1,1],lower=targus[1,2],upper=targus[1,3],power=targus[1,4],VEI=30.0,coverage='100%',Herds=h,VP=100))
df=rbind(df,data.frame(median=targus[2,1],lower=targus[2,2],upper=targus[2,3],power=targus[2,4],VEI=30.0,coverage='50%',Herds=h,VP=100))

}

control <- ICDIVA1df$Breaklength[ICDIVA1df$Breaklength>0]>60
casesW <- VES90VEI60IVP100DIVA1df$Breaklength[VES90VEI60IVP100DIVA1df$Breaklength>0]>60
casesM <- VES90VEI60IVP50DIVA1df$Breaklength[VES90VEI60IVP50DIVA1df$Breaklength>0]>60

targus = SweepByHerdP(casesW,casesM,control,50/2)

df=rbind(df,data.frame(median=targus[1,1],lower=targus[1,2],upper=targus[1,3],power=targus[1,4],VEI=60.0,coverage='100%',Herds=50,VP=100))
df=rbind(df,data.frame(median=targus[2,1],lower=targus[2,2],upper=targus[2,3],power=targus[2,4],VEI=60.0,coverage='50%',Herds=50,VP=100))

for(h in seq(50,2000,50))
{

targus = SweepByHerdP(casesW,casesM,control,h/2)

df=rbind(df,data.frame(median=targus[1,1],lower=targus[1,2],upper=targus[1,3],power=targus[1,4],VEI=60.0,coverage='100%',Herds=h,VP=100))
df=rbind(df,data.frame(median=targus[2,1],lower=targus[2,2],upper=targus[2,3],power=targus[2,4],VEI=60.0,coverage='50%',Herds=h,VP=100))

}

control <- ICDIVA1df$Breaklength[ICDIVA1df$Breaklength>0]>60
casesW <- VES90VEI90IVP100DIVA1df$Breaklength[VES90VEI90IVP100DIVA1df$Breaklength>0]>60
casesM <- VES90VEI90IVP50DIVA1df$Breaklength[VES90VEI90IVP50DIVA1df$Breaklength>0]>60

targus = SweepByHerdP(casesW,casesM,control,50/2)

df=rbind(df,data.frame(median=targus[1,1],lower=targus[1,2],upper=targus[1,3],power=targus[1,4],VEI=90.0,coverage='100%',Herds=50,VP=100))
df=rbind(df,data.frame(median=targus[2,1],lower=targus[2,2],upper=targus[2,3],power=targus[2,4],VEI=90.0,coverage='50%',Herds=50,VP=100))

for(h in seq(50,2000,50))
{

targus = SweepByHerdP(casesW,casesM,control,h/2)

df=rbind(df,data.frame(median=targus[1,1],lower=targus[1,2],upper=targus[1,3],power=targus[1,4],VEI=90.0,coverage='100%',Herds=h,VP=100))
df=rbind(df,data.frame(median=targus[2,1],lower=targus[2,2],upper=targus[2,3],power=targus[2,4],VEI=90.0,coverage='50%',Herds=h,VP=100))

}

df$VEI=as.factor(df$VEI)

p5=(ggplot(df,aes(x=Herds,y=100*power,linetype=VEI,col=coverage)) + geom_line()  +  scale_x_continuous('Herds') + scale_y_continuous('Power',limits=c(0.0,100)))+ guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Coverage'))
p6=(ggplot(df,aes(x=Herds,y=100*median,ymin=lower,ymax=upper,linetype=VEI,col=coverage)) + geom_line()  +  scale_x_continuous('Herds') + scale_y_continuous('% Effect Size',limits=c(0,15)))+ guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Coverage'))

pdf('VES90PowerDurationbyVI.pdf',height=2.5,width=6.5,family='Helvetica',pointsize=12)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,2)))
vplayout <- function(x,y)
{viewport(layout.pos.row=x,layout.pos.col=y)}
print(p5,vp=vplayout(1,2))
print(p6,vp=vplayout(1,1))
dev.off()
