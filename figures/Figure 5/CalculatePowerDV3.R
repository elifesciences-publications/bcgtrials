require(statmod)
require(ggplot2)
require(directlabels)

samples <- read.csv('../../parameters/samplesSOR.csv')[,2:12]
weight  <- scan('../../parameters/weightSOR.dat')

betaSOR=median(samples[,6])
qSOR=median(samples[,7])

samples <- read.csv('../../parameters/samplesSORI.csv')[,2:13]
weight  <- scan('../../parameters/weightSORI.dat')

betaSORI=median(samples[,6])
qSORI=median(samples[,7])
TO=median(samples[,5])
TR=median(samples[,12])

h=seq(30,400,1)

SORIDempR = c(0.0701609761,0.0001441427,0.8060078782,0.1236870030)

Hm=165

so = 1/TO
sr = 1/TR

# Prediction for transmission experiment

h=seq(30,800,5)
contact_t=seq(1.0,3.0,1)

df <- data.frame()

for(i in h)
{
for(j in contact_t)
{

R = j*betaSOR*i/((i/Hm)^qSOR)

df <- rbind(df,data.frame(h=i,t=j,model='SOR',R=R,Effect=75,Power=power.prop.test(p1=R/(R+2),p2=(0.25*R)/(0.25*R+2), n=i, sig.level=0.025, alternative="one.sided")$power))

df <- rbind(df,data.frame(h=i,t=j,model='SOR',R=R,Effect=50,Power=power.prop.test(p1=R/(R+2),p2=(0.5*R)/(0.5*R+2), n=i, sig.level=0.025, alternative="one.sided")$power))


df <- rbind(df,data.frame(h=i,t=j,model='SOR',R=R,Effect=25,Power=power.prop.test(p1=R/(R+2),p2=(0.75*R)/(0.75*R+2), n=i, sig.level=0.025, alternative="one.sided")$power))

R = (betaSORI*i/((i/Hm)^qSORI)) * (SORIDempR[2]*(j+(so/(sr*(so-sr)))*(exp(-sr*j)-1) - (sr/(so*(so-sr)))*(exp(-so*j)-1)) + SORIDempR[3]*(j+(1/sr)*(exp(-sr*j)-1)) + SORIDempR[4]*j)

df <- rbind(df,data.frame(h=i,t=j,model='SORI',R=R,Effect=75, Power=power.prop.test(p1=R/(R+2),p2=(0.25*R)/(0.25*R+2), n=i, sig.level=0.025, alternative="one.sided")$power))

df <- rbind(df,data.frame(h=i,t=j,model='SORI',R=R,Effect=50, Power=power.prop.test(p1=R/(R+2),p2=(0.5*R)/(0.5*R+2), n=i, sig.level=0.025, alternative="one.sided")$power))

df <- rbind(df,data.frame(h=i,t=j,model='SORI',R=R,Effect=25, Power=power.prop.test(p1=R/(R+2),p2=(0.75*R)/(0.75*R+2), n=i, sig.level=0.025, alternative="one.sided")$power))

}
}

# Colour blind palette with black:
cbPaletteB <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

df$t[df$t==1] = '1 Year'
df$t[df$t==2] = '2 Years'
df$t[df$t==3] = '3 Years'


pdf('Figure 5.pdf',family='Helvetica',points=12,width=10,height=5)
p1=ggplot(df[df$Effect==75,],aes(x=h,y=100*Power,col=as.factor(t),shape=as.factor(t),lty=as.factor(t)))+geom_line()+geom_point(size=0.75)+facet_wrap(~model) + xlab('Group Size') + xlim(c(0,300)) + guides(linetype = guide_legend(order=1,title='Contact Time'),colour=guide_legend(order=1,title='Contact Time')) + scale_colour_manual(values=cbPaletteB) + geom_hline(yintercept=80,col='Grey')+ylab('Power')
print(p1)
dev.off()

pdf('Figure 5 - supplement 1.pdf',family='Helvetica',points=12,width=10,height=5)
p1=ggplot(df[df$Effect==50,],aes(x=h,y=Power*100,col=as.factor(t),shape=as.factor(t),lty=as.factor(t)))+geom_line()+geom_point(size=0.75)+facet_wrap(~model) + xlab('Group Size') + xlim(c(0,500)) + guides(linetype = guide_legend(order=1,title='Contact Time'),colour=guide_legend(order=1,title='Contact Time')) + scale_colour_manual(values=cbPaletteB)+ geom_hline(yintercept=80,col='Grey')+ylab('Power')
print(p1)
dev.off()

pdf('Figure 5 - supplement 2.pdf',family='Helvetica',points=12,width=10,height=5)
p1=ggplot(df[df$Effect==25,],aes(x=h,y=Power*100,shape=as.factor(t),col=as.factor(t),lty=as.factor(t)))+geom_line()+geom_point(size=0.75)+facet_wrap(~model) + xlab('Group Size') + xlim(c(0,800)) + guides(linetype = guide_legend(order=1,title='Contact Time'),colour=guide_legend(order=1,title='Contact Time')) + scale_colour_manual(values=cbPaletteB)+ geom_hline(yintercept=80,col='Grey')+ylab('Power')
print(p1)
dev.off()



