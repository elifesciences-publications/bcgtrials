

load('../../simulations/Duration/SOR/DurationSOR.RData')
    
cbPalette1 <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

out_df$duration[out_df$duration==3] = '3 Years'
out_df$duration[out_df$duration==6] = '6 Years'
out_df$duration[out_df$duration==9] = '9 Years'


require(ggplot2)
require(grid)
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
    

p1=ggplot(out_df[out_df$prop==0 & out_df$Vaccinated==0 & is.element(out_df$duration,c('3 Years','6 Years','9 Years')),],aes(Reactor/AtRisk,linetype=retain,fill=retain)) + geom_density(alpha=0.5)+facet_wrap(~as.factor(duration)) + scale_fill_manual(values=cbPalette1) + guides(linetype = guide_legend(order=1,title='Retain Reactors'),fill=guide_legend(order=1,title='Retain Reactors')) + xlab('') + xlim(0,0.5) 

p2=ggplot(out_df[out_df$retain==TRUE & out_df$Vaccinated == 0 & is.element(out_df$duration,c('3 Years','6 Years','9 Years')),],aes(Reactor/AtRisk,linetype=as.factor(prop),fill=as.factor(prop)))+geom_density(alpha=0.5)+facet_wrap(~as.factor(duration)) + scale_fill_manual(values=(cbPalette1[3:4])) + guides(linetype = guide_legend(order=1,title='Coverage (%)'),fill=guide_legend(order=1,title='Coverage (%)')) + xlab('Within-herd prevalance') + xlim(0,0.5)

pdf('FigureR-supplement1.pdf',height=1.5*4,width=8.5,family='Helvetica',pointsize=12)
grid.arrange(arrangeGrob(p1,top=A.grob),arrangeGrob(p2,top=B.grob), ncol = 1, 
             layout_matrix = rbind(c(1),c(2)))
dev.off()
