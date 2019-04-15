library(lattice)

setwd("~/Dropbox/Lomas/Nat Geo")
anat<-read.csv('Anatomia.csv')#Datos de anatomía de Camilo
summary(anat)

xyplot(Corteza~Altura|specie,data=anat)
plot(Corteza+Estrato.Algal~Altura,data=anat[anat$specie=='Roccellina',])

par(mfrow=c(1,1),mar=c(5,7,1,1))
plot(Corteza+Estrato.Algal~Altura,data=anat[anat$specie=='Roccellina',],pch=16,cex.lab=1.5,cex.axis=1.5,ylim=c(0,160),ylab='Grosor de capa de almacenamiento de agua
     (corteza + estrato algal)')
#plot(Médula.total~Altura,data=anat[anat$specie=='Roccellina',],pch=16,cex.lab=1.5,ylim=c(0,400),ylab='Grosor (um)')
points(Corteza~Altura,data=anat[anat$specie=='Roccellina',],pch=16,col='blue')
points(Estrato.Algal~Altura,data=anat[anat$specie=='Roccellina',],pch=16,col='green')
summary(lm(Corteza~Altura,data=anat[anat$specie=='Roccellina',]))
summary(lm(Estrato.Algal~Altura,data=anat[anat$specie=='Roccellina',]))
abline(a=128,b=-0.096,col='green',lwd=2,lty=2)
text(450,100,'p = 0.025',col='green',cex=2)
summary(lm(Médula.total~Altura,data=anat[anat$specie=='Roccellina',]))
summary(lm(Corteza+Estrato.Algal~Altura,data=anat[anat$specie=='Roccellina',]))

#Plots ternarios
par(mfrow=c(1,1))
TernaryPlot('Corteza','Algal','Médula')
AddToTernary(points,anat[anat$specie=='Niebla',5:7],col=1,pch=16)
AddToTernary(points,anat[anat$specie=='Roccellina',5:7],col=2,pch=16)
AddToTernary(points,anat[anat$specie=='Heterodermia',5:7],col=3,pch=16)
AddToTernary(points,anat[anat$specie=='Buellia',5:7],col=4,pch=16)

par(mfrow=c(2,2))
TernaryPlot('Corteza','Algal','Médula',main='Niebla')
AddToTernary(points,anat[anat$specie=='Niebla'&anat$Altura<601,5:7], col=2,pch=16)
AddToTernary(points,anat[anat$specie=='Niebla'&anat$Altura>600,5:7], col=4,pch=16)
AddToTernary(points,anat[anat$specie=='Niebla'&anat$Altura>799,5:7], col=1,pch=16)

TernaryPlot('Corteza','Algal','Médula',main='Roccellina')
AddToTernary(points,anat[anat$specie=='Roccellina'&anat$Altura<601,5:7], col=2,pch=16)
AddToTernary(points,anat[anat$specie=='Roccellina'&anat$Altura>600,5:7], col=4,pch=16)
AddToTernary(points,anat[anat$specie=='Roccellina'&anat$Altura>799,5:7], col=1,pch=16)

TernaryPlot('Corteza','Algal','Médula',main='Heterodermia')
AddToTernary(points,anat[anat$specie=='Heterodermia'&anat$Altura<601,5:7], col=2,pch=16)
AddToTernary(points,anat[anat$specie=='Heterodermia'&anat$Altura>600,5:7], col=4,pch=16)
AddToTernary(points,anat[anat$specie=='Heterodermia'&anat$Altura>799,5:7], col=1,pch=16)

TernaryPlot('Corteza','Algal','Médula',main='Buellia')
AddToTernary(points,anat[anat$specie=='Buellia'&anat$Altura<601,5:7], col=2,pch=16)
AddToTernary(points,anat[anat$specie=='Buellia'&anat$Altura>600,5:7], col=4,pch=16)
AddToTernary(points,anat[anat$specie=='Buellia'&anat$Altura>799,5:7], col=1,pch=16)
