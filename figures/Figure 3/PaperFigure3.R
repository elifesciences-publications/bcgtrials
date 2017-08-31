load('../../simulations/Efficacy/SORI/VEISweepSORI.RData')

# Colour blind palette with grey:
cbPaletteG <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
# Colour blind palette with black:
cbPaletteB <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


source('PowerCompareEfficacyMeasuresWithVI.R')
source('PowerCompareEfficacyMeasuresWithVIVES30.R')
source('PowerCompareEfficacyMeasuresWithVIVES90.R')

save.image('PaperFigure3.RData')

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

 

pdf('Figure3.pdf',height=2.5*4,width=8.5,family='Helvetica',pointsize=12)
grid.arrange(arrangeGrob(p4,top=A.grob),arrangeGrob(p3,top=B.grob), arrangeGrob(p2,top=C.grob),arrangeGrob(p1,top=D.grob),arrangeGrob(p6,top=E.grob),arrangeGrob(p5,top=F.grob), ncol = 2, 
             layout_matrix = cbind(c(1,3,5),c(2,4,6))) 
dev.off()





