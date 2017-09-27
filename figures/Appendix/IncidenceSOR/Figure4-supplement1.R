require(ggplot2)
require(cowplot)

load('../../../simulations/Persistence/SOR/VEISweepSOR.RData')

# Colour blind palette with grey:
cbPaletteG <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
# Colour blind palette with black:
cbPaletteB <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

source('PowerHerdLevelIncidenceVES30.R')
source('PowerHerdLevelIncidenceVES60.R')
source('PowerHerdLevelIncidenceVES90.R')

save.image('PaperFigure4-supplement1.RData')

p1=(ggplot(df30,aes(x=Herds,y=100*power,linetype=VEI,col=coverage,shape=coverage)) + geom_line()  +  scale_x_continuous('Herds') + scale_y_continuous('Power',limits=c(0.0,100)))+ guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Coverage'),shape=guide_legend(order=2,title='Coverage'))+ scale_colour_manual(values=cbPaletteB) + geom_point(size=0.75)
p2=(ggplot(df30,aes(x=Herds,y=100*median,ymin=lower,ymax=upper,linetype=VEI,col=coverage,shape=coverage)) + geom_line()  +  scale_x_continuous('Herds') + scale_y_continuous('% Effect Size',limits=c(0,15)))+ guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Coverage'),shape=guide_legend(order=2,title='Coverage'))+ scale_colour_manual(values=cbPaletteB) + geom_point(size=0.75)

p3=(ggplot(df60,aes(x=Herds,y=100*power,linetype=VEI,col=coverage,shape=coverage)) + geom_line()  +  scale_x_continuous('Herds') + scale_y_continuous('Power',limits=c(0.0,100)))+ guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Coverage'),shape=guide_legend(order=2,title='Coverage'))+ scale_colour_manual(values=cbPaletteB) + geom_point(size=0.75)
p4=(ggplot(df60,aes(x=Herds,y=100*median,ymin=lower,ymax=upper,linetype=VEI,col=coverage,shape=coverage)) + geom_line()  +  scale_x_continuous('Herds') + scale_y_continuous('% Effect Size',limits=c(0,15)))+ guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Coverage'),shape=guide_legend(order=2,title='Coverage'))+ scale_colour_manual(values=cbPaletteB) + geom_point(size=0.75)

p5=(ggplot(df90,aes(x=Herds,y=100*power,linetype=VEI,col=coverage,shape=coverage)) + geom_line()  +  scale_x_continuous('Herds') + scale_y_continuous('Power',limits=c(0.0,100)))+ guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Coverage'),shape=guide_legend(order=2,title='Coverage'))+ scale_colour_manual(values=cbPaletteB) + geom_point(size=0.75)
p6=(ggplot(df90,aes(x=Herds,y=100*median,ymin=lower,ymax=upper,linetype=VEI,col=coverage,shape=coverage)) + geom_line()  +  scale_x_continuous('Herds') + scale_y_continuous('% Effect Size',limits=c(0,15)))+ guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Coverage'),shape=guide_legend(order=2,title='Coverage'))+ scale_colour_manual(values=cbPaletteB) + geom_point(size=0.75)

legend <- get_legend(p1 + theme(legend.position="right"))
 
pbox <- plot_grid(p4 + theme(legend.position = "none"),p3 + theme(legend.position = "none"),p2 + theme(legend.position = "none"),p1 + theme(legend.position = "none"),p6 + theme(legend.position = "none"),p5 + theme(legend.position = "none"),labels=c('A','B','C','D','E','F'),hjust=-1,nrow=3)

power_plot <- plot_grid(pbox,legend,ncol=2,rel_widths=c(1,0.2))

pdf('Figure-Incidence-Power.pdf',height=2.5*4,width=8.5,family='Helvetica',pointsize=12)
print(power_plot)
dev.off()

pdf('Figure-Incidence-Effect.pdf',height=2.5*2,width=2.5*2,family='Helvetica',pointsize=12)
(ggplot(df90,aes(x=Herds,y=100*median,ymin=100*lower,ymax=100*upper,fill=as.factor(coverage),col=as.factor(coverage)))  +  scale_x_continuous('Herds',limits=c(50,300)) + scale_y_continuous('Efficacy')) + scale_colour_manual(values=cbPaletteB) + scale_fill_manual(values=cbPaletteB)+ geom_ribbon(alpha=0.5,col=NA) + geom_line() + guides(linetype = guide_legend(order=1,title='Coverage'),colour=guide_legend(order=1,title='Coverage'),fill=guide_legend(order=1,title='Coverage'))+ background_grid(major = "xy", minor = "xy") + facet_wrap(~VEI)
dev.off()
       

