#removed setting working directory so paths can be relative for each file read in. Not sure this is the best
#way to do it but allows a user to pull the repo and not have to change working directory
#setwd("~/Dropbox/Lomas/Nat Geo")
areas<-read.csv('data/Images_Wet_Dry.csv')
summary(areas)
#The values for area are not in units of cm for some reason. 
areas$Area_D<-areas$Area_D/10
areas$Area_W<-areas$Area_W/10
areas$Area_Tot<-areas$Area_Tot/10
areas$Water<-(areas$Mass_W-areas$Mass_D)/areas$Mass_D
areas$Swelling<-areas$Area_W/areas$Area_D
areas$A_Ratio<-areas$Area_Tot/areas$Area_W
areas$TMA_V<-areas$Mass_D/areas$Area_D
areas$WHC<-(areas$Mass_W-areas$Mass_D)/areas$Area_W
areas$WHC_surf<-(areas$Mass_W-areas$Mass_D)/areas$Area_Tot

summary(lm(Swelling~Species-1,data=areas))
summary(lm(TMA_V~Species-1,data=areas))
summary(lm(A_Ratio~Species*Elevation,data=areas))
summary(lm(WHC~Species,data=areas))
plot(WHC~Species,data=areas)
par(mfrow=c(1,1),mar=c(5,5,1,1))
plot(WHC~Species,data=areas,xlim=c(1.5,6.5),ylab='Water Holding capacity (g H2O/cm2)',cex.axis=1.25,cex.lab=1.5)
plot(WHC_surf~Species,data=areas,xlim=c(1.5,6.5),ylab='Water Holding capacity (g H2O/cm2)',cex.axis=1.25,cex.lab=1.5)
plot(A_Ratio~Species,data=areas,ylim=c(1,4),xlim=c(1.5,6.5),ylab='Total Area: Projected Area',cex.axis=1.25,cex.lab=1.5)

library(lattice)
rocc<-areas[areas$Species=='Roccellina',]
plot(rocc)
summary(lm(Swelling~Elevation,data=rocc))


summary(lm(WHC~Elevation,data=areas[areas$Species=='Roccellina',]))
summary(lm(TMA_V~Elevation,data=areas[areas$Species=='Roccellina',]))

plot(WHC~Elevation,data=areas[areas$Species=='Roccellina',],ylab='Water Holding capacity (g H2O/cm2)',main='Roccellina cerebriformis',xlab='Elevation (m)',cex.axis=1.5,cex.lab=1.5,ylim=c(0,0.031))
abline(a=2.474e-02,b=-2.849e-05,lty=2)
text(650,0.027,'r2 = 0.50, p « 0.001',cex=1.5)
