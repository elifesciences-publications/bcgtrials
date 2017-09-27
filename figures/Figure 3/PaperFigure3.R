require(gridExtra)
require(ggplot2)
require(grid)

SweepByHerd_MixedD1 <-function(cases,pop,baselinec,baselinep,controlc,controlp,nherds,pmix)
{

ARV = numeric(10000)

ARV = sapply(rep(nherds,10000),function(x,cases,pop,baselinec,baselinep,controlc,controlp,pmix)
{

y=as.integer(pmix*x)
z=x-y


m<-sample(1:length(cases),y,replace=TRUE)
n<-sample(1:length(controlc),y,replace=TRUE)

# Relative Risk Ratio
oot = sum(cases[m])/sum(pop[m]) / (sum(controlc[n])/sum(controlp[n]))
# Standard error of (log) relative risk ratio
poot = sqrt((1/sum(cases[m]))-(1/sum(pop[m]))+(1/sum(controlc[n]))-(1/sum(controlp[n])))


return(cbind(oot,poot))
},cases=cases,pop=pop,baselinec=baselinec,baselinep=baselinep,controlc=controlc,controlp=controlp,pmix=pmix)
ARV=t(ARV)

# Aggregate over all simulations to calculated expected effect size
expected_effect = sum(cases)/sum(pop) / (sum(controlc)/sum(controlp))
powerE = mean(ARV[,1] > (expected_effect - 0.05*expected_effect) & ARV[,1] < (expected_effect + 0.05*expected_effect))

zcr = qnorm(p = 1-0.025, mean = 0, sd = 1)
power = mean((log(ARV[,1])/ARV[,2] < (-zcr)))

return(data.frame(nherds,power,powerE))
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

# Aggregate over all simulations to calculated expected effect size
expected_effect = sum(cases)/sum(pop) / (sum(baselinec)/sum(baselinep))
powerE = mean(ARV[,1] > (expected_effect - 0.05*expected_effect) & ARV[,1] < (expected_effect + 0.05*expected_effect))

zcr = qnorm(p = 1-0.025, mean = 0, sd = 1)
power = mean((log(ARV[,1])/ARV[,2] < (-zcr)))

return(data.frame(nherds,power,powerE))

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
effect=mean(ARV)
rangeARV=quantile(ARV,c(0.025,0.975))
return(data.frame(nherds,effect=effect,effectL=rangeARV[1],effectU=rangeARV[2]))
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
effect=mean(ARV)
rangeARV=quantile(ARV,c(0.025,0.975))
return(data.frame(nherds,effect=effect,effectL=rangeARV[1],effectU=rangeARV[2]))
}


load('../../simulations/Efficacy/SORI/VEISweepSORI.RData')

# Colour blind palette with grey:
cbPaletteG <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
# Colour blind palette with black:
cbPaletteB <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


source('PowerCompareEfficacyMeasuresWithVI.R')
source('PowerCompareEfficacyMeasuresWithVIVES30.R')
source('PowerCompareEfficacyMeasuresWithVIVES90.R')

save.image('PaperFigure3.RData')

require(cowplot)

p2=(ggplot(dfE60,aes(x=Herds,y=100*effect,col=measure,shape=measure,linetype=VE)) + geom_line() + geom_point(size=0.75) +  scale_x_continuous('Herds',limits=c(51,300)) + scale_y_continuous('Efficacy',limits=c(-20,100))) + scale_colour_manual(values=cbPaletteB)  + guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Measures'))+ background_grid(major = "xy", minor = "xy")


p1=(ggplot(df60,aes(x=Herds,y=100*power,col=measure,shape=measure,linetype=VE)) + geom_line() + geom_point(size=0.75) +  scale_x_continuous('Herds',limits=c(51,300)) + scale_y_continuous('Power',limits=c(0,100))) + scale_colour_manual(values=cbPaletteB) + guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Measures')) + background_grid(major = "xy", minor = "xy")


p1e=(ggplot(df60,aes(x=Herds,y=100*powerE,col=measure,shape=measure,linetype=VE)) + geom_line() + geom_point(size=0.75) +  scale_x_continuous('Herds',limits=c(51,300)) + scale_y_continuous('Power',limits=c(0,100.0))) + scale_colour_manual(values=cbPaletteB) + guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Measures')) + background_grid(major = "xy", minor = "xy")



p4=(ggplot(dfE30,aes(x=Herds,y=100*effect,col=measure,shape=measure,linetype=VE)) + geom_line() + geom_point(size=0.75) +  scale_x_continuous('Herds',limits=c(51,300)) + scale_y_continuous('Efficacy',limits=c(-20,100))) + scale_colour_manual(values=cbPaletteB)  + guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Measures')) + background_grid(major = "xy", minor = "xy")


p3=(ggplot(df30,aes(x=Herds,y=100*power,col=measure,shape=measure,linetype=VE)) + geom_line() + geom_point(size=0.75) +  scale_x_continuous('Herds',limits=c(51,300)) + scale_y_continuous('Power',limits=c(0,100))) + scale_colour_manual(values=cbPaletteB)  + guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Measures')) + background_grid(major = "xy", minor = "xy")


p3e=(ggplot(df30,aes(x=Herds,y=100*powerE,col=measure,shape=measure,linetype=VE)) + geom_line() + geom_point(size=0.75) +  scale_x_continuous('Herds',limits=c(51,300)) + scale_y_continuous('Power',limits=c(0,100))) + scale_colour_manual(values=cbPaletteB)  + guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Measures')) + background_grid(major = "xy", minor = "xy")



p6=(ggplot(dfE90,aes(x=Herds,y=100*effect,col=measure,shape=measure,linetype=VE,shape=measure)) + geom_line() + geom_point(size=0.75)  +  scale_x_continuous('Herds',limits=c(51,300)) + scale_y_continuous('Efficacy',limits=c(-20,100))) + scale_colour_manual(values=cbPaletteB) + guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Measures')) + background_grid(major = "xy", minor = "xy")


p5=(ggplot(df90,aes(x=Herds,y=power,col=measure,shape=measure,linetype=VE,shape=measure)) + geom_line() + geom_point(size=0.75) +  scale_x_continuous('Herds',limits=c(51,300)) + scale_y_continuous('Power')) + scale_colour_manual(values=cbPaletteB)+ guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Measures'))+ background_grid(major = "xy", minor = "xy")


p5e=(ggplot(df90,aes(x=Herds,y=100*powerE,col=measure,shape=measure,linetype=VE,shape=measure)) + geom_point(size=0.75) + geom_line()  +  scale_x_continuous('Herds',limits=c(51,300)) + scale_y_continuous('Power',limits=c(0,100))) + scale_colour_manual(values=cbPaletteB)+ guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Measures'))+ background_grid(major = "xy", minor = "xy")

 
legend <- get_legend(p1 + theme(legend.position="right"))
 
pbox <- plot_grid(p4 + theme(legend.position = "none"),p3 + theme(legend.position = "none"),p2 + theme(legend.position = "none"),p1 + theme(legend.position = "none"),p6 + theme(legend.position = "none"),p5 + theme(legend.position = "none"),labels=c('A','B','C','D','E','F'),hjust=-1,nrow=3)


pdf('Figure3.pdf',height=2.5*4,width=8.5,family='Helvetica',pointsize=12)
print(plot_grid(pbox,legend,ncol=2,rel_widths=c(1,0.2)))
dev.off()

pdf('Figure3-supplement1.pdf',height=2.5*2,width=2.5*2,family='Helvetica',pointsize=12)
(ggplot(dfE90[dfE90$VE==90,],aes(x=Herds,y=100*effect,ymin=100*effectL,ymax=100*effectU,fill=measure,col=measure))  +  scale_x_continuous('Herds',limits=c(50,300)) + scale_y_continuous('Efficacy')) + scale_colour_manual(values=cbPaletteB) + scale_fill_manual(values=cbPaletteB)+ geom_ribbon(alpha=0.5,col=NA) + facet_wrap(~measure) + geom_line() + guides(linetype = FALSE,colour=FALSE,fill=FALSE)+ background_grid(major = "xy", minor = "xy")
dev.off()



