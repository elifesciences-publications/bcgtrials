load('../../simulations/Persistence/SORI/VEISweepSORIP.RData')

# Colour blind palette with grey:
cbPaletteG <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
# Colour blind palette with black:
cbPaletteB <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


source('PowerHerdLevelIncidenceVES30.R')
source('PowerHerdLevelIncidenceVES60.R')
source('PowerHerdLevelIncidenceVES90.R')

save.image('PaperFigure4.RData')

p1=(ggplot(df30,aes(x=Herds,y=100*power,linetype=VEI,col=coverage,shape=coverage)) + geom_line()  +  scale_x_continuous('Herds') + scale_y_continuous('Power',limits=c(0.0,100)))+ guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Coverage'),shape=guide_legend(order=2,title='Coverage')) + scale_colour_manual(values=cbPaletteB) + geom_point()

p2=(ggplot(df30,aes(x=Herds,y=100*median,ymin=lower,ymax=upper,linetype=VEI,col=coverage,shape=coverage)) + geom_line()  +  scale_x_continuous('Herds') + scale_y_continuous('% Effect Size',limits=c(0,15)))+ guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Coverage'),shape=guide_legend(order=2,title='Coverage'))+ scale_colour_manual(values=cbPaletteB) + geom_point()

p3=(ggplot(df60,aes(x=Herds,y=100*power,linetype=VEI,col=coverage)) + geom_line()  +  scale_x_continuous('Herds') + scale_y_continuous('Power',limits=c(0.0,100)))+ guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Coverage'),shape=guide_legend(order=2,title='Coverage'))+ scale_colour_manual(values=cbPaletteB) + geom_point()

p4=(ggplot(df60,aes(x=Herds,y=100*median,ymin=lower,ymax=upper,linetype=VEI,col=coverage)) + geom_line()  +  scale_x_continuous('Herds') + scale_y_continuous('% Effect Size',limits=c(0,20)))+ guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Coverage'),shape=guide_legend(order=2,title='Coverage'))+ scale_colour_manual(values=cbPaletteB) + geom_point()

p5=(ggplot(df90,aes(x=Herds,y=100*power,linetype=VEI,col=coverage,shape=coverage)) + geom_line()  +  scale_x_continuous('Herds') + scale_y_continuous('Power',limits=c(0.0,100)))+ guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Coverage')) + scale_colour_manual(values=cbPaletteB) + geom_point()

p6=(ggplot(df90,aes(x=Herds,y=100*median,ymin=lower,ymax=upper,linetype=VEI,col=coverage,shape=coverage)) + geom_line()  +  scale_x_continuous('Herds') + scale_y_continuous('% Effect Size',limits=c(0,15)))+ guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Coverage'),shape=guide_legend(order=2,title='Coverage')) + scale_colour_manual(values=cbPaletteB) + geom_point()



require(gridExtra)

A.grob <- textGrob(
    label = "A",
    x = unit(1, "lines"), 
    y = unit(0, "lines"),
    hjust = 0, vjust = 0,
    gp = gpar(fontsize = 12))

B.grob <- textGrob(
    label = "B",
    x = unit(1, "lines"), 
    y = unit(0, "lines"),
    hjust = 0, vjust = 0,
    gp = gpar(fontsize = 12))
    
C.grob <- textGrob(
    label = "C",
    x = unit(1, "lines"), 
    y = unit(0, "lines"),
    hjust = 0, vjust = 0,
    gp = gpar(fontsize = 12))

D.grob <- textGrob(
    label = "D",
    x = unit(1, "lines"), 
    y = unit(0, "lines"),
    hjust = 0, vjust = 0,
    gp = gpar(fontsize = 12))

E.grob <- textGrob(
    label = "E",
    x = unit(1, "lines"), 
    y = unit(0, "lines"),
    hjust = 0, vjust = 0,
    gp = gpar(fontsize = 12))

F.grob <- textGrob(
    label = "F",
    x = unit(1, "lines"), 
    y = unit(0, "lines"),
    hjust = 0, vjust = 0,
    gp = gpar(fontsize = 12))
    

pdf('Figure4.pdf',height=2.5*4,width=8.5,family='Helvetica',pointsize=12)
grid.arrange(arrangeGrob(p2,top=A.grob),arrangeGrob(p1,top=B.grob), arrangeGrob(p4,top=C.grob),arrangeGrob(p3,top=D.grob),arrangeGrob(p6,top=E.grob),arrangeGrob(p5,top=F.grob), ncol = 2, 
             layout_matrix = cbind(c(1,3,5),c(2,4,6)))
dev.off()
