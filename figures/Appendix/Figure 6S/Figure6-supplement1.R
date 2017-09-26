require(ggplot2)
require(grid)
require(gridExtra)

load('../../simulations/Persistence/SOR/VEISweepSOR.RData')

# Colour blind palette with grey:
cbPaletteG <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
# Colour blind palette with black:
cbPaletteB <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

source('PowerHerdLevelRecurrenceVES30.R')
source('PowerHerdLevelRecurrenceVES60.R')
source('PowerHerdLevelRecurrenceVES90.R')

p1=(ggplot(df30,aes(x=Herds,y=100*power,linetype=VEI,col=coverage)) + geom_line()  +  scale_x_continuous('Herds') + scale_y_continuous('Power',limits=c(0.0,100)))+ guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Coverage'))+ scale_colour_manual(values=cbPaletteB)
p2=(ggplot(df30,aes(x=Herds,y=100*median,ymin=100*lower,ymax=100*upper,linetype=VEI,col=coverage)) + geom_line()  +  scale_x_continuous('Herds') + scale_y_continuous('% Effect Size',limits=c(0,15)))+ guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Coverage'))+ scale_colour_manual(values=cbPaletteB)

p3=(ggplot(df60,aes(x=Herds,y=100*power,linetype=VEI,col=coverage)) + geom_line()  +  scale_x_continuous('Herds') + scale_y_continuous('Power',limits=c(0.0,100)))+ guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Coverage'))+ scale_colour_manual(values=cbPaletteB)
p4=(ggplot(df60,aes(x=Herds,y=100*median,ymin=100*lower,ymax=100*upper,linetype=VEI,col=coverage)) + geom_line()  +  scale_x_continuous('Herds') + scale_y_continuous('% Effect Size',limits=c(0,15)))+ guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Coverage'))+ scale_colour_manual(values=cbPaletteB)

p5=(ggplot(df90,aes(x=Herds,y=100*power,linetype=VEI,col=coverage)) + geom_line()  +  scale_x_continuous('Herds') + scale_y_continuous('Power',limits=c(0.0,100)))+ guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Coverage'))+ scale_colour_manual(values=cbPaletteB)
p6=(ggplot(df90,aes(x=Herds,y=100*median,ymin=100*lower,ymax=100*upper,linetype=VEI,col=coverage)) + geom_line()  +  scale_x_continuous('Herds') + scale_y_continuous('% Effect Size',limits=c(0,15)))+ guides(linetype = guide_legend(order=1,title=expression(epsilon[I])),colour=guide_legend(order=2,title='Coverage'))+ scale_colour_manual(values=cbPaletteB)


save.image('PaperFigure6-supplement1.RData')

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
    

pdf('Figure6-supplement1.pdf',height=2.5*4,width=8.5,family='Helvetica',pointsize=12)
grid.arrange(arrangeGrob(p2,top=A.grob),arrangeGrob(p1,top=B.grob), arrangeGrob(p4,top=C.grob),arrangeGrob(p3,top=D.grob),arrangeGrob(p6,top=E.grob),arrangeGrob(p5,top=F.grob), ncol = 2, 
             layout_matrix = cbind(c(1,3,5),c(2,4,6)))
dev.off()
