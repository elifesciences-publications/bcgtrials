require(ggplot2)

SORIempR =c(0.0170336134,0.0001009178,0.2425210345,0.7403444343) 
SORIDempR = c(0.0701609761,0.0001441427,0.8060078782,0.1236870030)

df=data.frame(Model=rep('SORI 2012',4),state=c('S','O','R','I'),freq=SORIempR)
df$state=factor(df$state,as.character(df$state))

df=rbind(df,data.frame(Model=rep('SORI 2015',4),state=c('S','O','R','I'),freq=SORIDempR))

pdf('EmpiricalR.pdf',family='Helvetica',pointsize=12)
print(ggplot(df,aes(x=as.factor(state),y=freq,fill=Model))+geom_bar(stat='identity',position='dodge')+xlab('Epidemiological State')+ylab('Frequency'))
dev.off()