load('DarthTablesV2.RData')

PTI1R<-hist(Breakdown_Table$Reactors_at_Start[Breakdown_Table$testinterval==1],breaks=c(seq(0,10),1000),right=FALSE)$counts
PTI2R<-hist(Breakdown_Table$Reactors_at_Start[Breakdown_Table$testinterval==2],breaks=c(seq(0,10),1000),right=FALSE)$counts
PTI4R<-hist(Breakdown_Table$Reactors_at_Start[Breakdown_Table$testinterval==4],breaks=c(seq(0,10),1000),right=FALSE)$counts

Rdistro = rbind(PTI1R/sum(PTI1R),PTI2R/sum(PTI2R),PTI4R/sum(PTI4R))
colnames(Rdistro)<-c('0','1','2','3','4','5','6','7','8','9','10+')

pdf('ReactorDistro.pdf')
barplot(Rdistro,beside=T,legend.text=c('PTI 1', 'PTI 2', 'PTI 4'))
dev.off()