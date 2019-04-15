####Overal code for Patache community data, combining various data forms into master files for environment, widely measured traits and community composition###
setwd("~/Dropbox/Lomas/Nat Geo/Combined Data")

##Load libraries
library(smatr,lattice,nlme,labdsv,vegan)
#########Emphasis on Transect B########

##########Taxa####
taxa<-read.csv('Taxa_TransectB.csv')
taxa$Transect<-'B'
taxa$Elevation<-NA #this creates an empty column of data for storing the elevation records from the collection name
taxa$Quadrat<-NA
taxa$Number<-NA

for (i in seq(2,nrow(taxa),1)){ #We are using a for loop to go through the data file row by row, and parse the collection name into Elevation, Quadrat and Number identifiers to align with other datasets
  #first we will parse the "Label" column with the photo file name, to extract plate info:
  taxa[i,10]<-unlist(strsplit(as.character(taxa[i,3]),c("_")))[2] #for example here we remove the second element from the collection name (Elevation) and asign that value to the Elevation column of the data spreadsheet
  taxa[i,11]<-unlist(strsplit(as.character(taxa[i,3]),"_"))[3]
  taxa[i,12]<-unlist(strsplit(as.character(taxa[i,3]),"_"))[4]
}
head(taxa)
summary(taxa)

##############################
#################Cover Master#############
cover<-read.csv('CoverMaster.csv')
summary(cover)
coverB<-cover[cover$Transect=='B',] #just transect B (richer trait data)


#############Contact Angle Data############
angles<-read.csv('ContactAngles_TransectB.csv')
angles$Transect<-'B'
summary(angles)
summary(lme(Angle~Elevation,~1|Number/Replicate/Subreplicate,data=angles,na.action=na.omit))
summary(lme(Time_Absorption~Elevation,~1|Number/Replicate/Subreplicate,data=angles,na.action=na.omit))
ang_rep<-aggregate(angles[,c(7:9)],by=list(Elevation=angles$Elevation,Quadrat=angles$Quadrat,Number=angles$Number,Replicate=angles$Replicate),mean,na.action=na.omit)
ang_num<-aggregate(ang_rep[,c(5:7)],by=list(Elevation=ang_rep$Elevation,Quadrat=ang_rep$Quadrat,Number=ang_rep$Number),mean,na.action=na.omit)
plot(Time_Surface~Elevation,data=ang_num)
sunflowerplot(Time_Surface~Elevation,data=ang_num,ylab='Time to start of uptake (s)',xlab='Elevation (m)',cex.lab=1.5,cex=1.25,cex.axis=1.25)
abline(h=87,lwd=1.5,lty=2)
text(350,90,'Hydrophobic',cex=1.2)

sunflowerplot(Time_Absorption~Elevation,data=ang_num,ylab='Time to droplet absorption (s)',xlab='Elevation (m)',cex.lab=1.5,cex=1.25,cex.axis=1.25)
abline(h=87,lwd=1.5,lty=2)
text(350,90,'Hydrophobic',cex=1.2)

sunflowerplot(Time_Absorption-Time_Surface~Elevation,data=ang_num,ylab='Absorption time(s)',xlab='Elevation (m)',cex.lab=1.5,cex=1.25,cex.axis=1.25)

plot(Angle~Elevation,data=ang_num)

#merge this with cover data to scale by abundance
ang_tax<-merge(ang_num,taxa[,c(1,2,3,8:12)],all=T)


cov_ang<-merge(coverB[,c(1:3,5,6)],ang_tax[,c(4:9,11)]) #merged dataset
sp_trait<-read.csv('Spp.csv')

cov_trait<-merge(cov_ang,sp_trait)

sunflowerplot(Time_Absorption~Elevation,data=cov_ang[,])
sunflowerplot(Time_Absorption~Elevation,data=cov_ang[!cov_ang$Cover=='T',])

##################

########Quadrat Summary Info ###### roughness, aspect, etc
quad<-read.csv('Quadrat_Master.csv')

quad_cov<-merge(quad,cov_ang,all=T)

#################
#########################Species perspective######################
spp_summ<-aggregate(coverB[,c(1)],coverB[,c(1,2,4,6)],length) #number of quadrats for each species/elevation
names(spp_summ)[5]<-'N_Quad_Pres'


sp_persp<-merge(ang_tax,spp_summ)
sp_persp<-merge(sp_persp,sp_trait)

plot(Time_Absorption~jitter(as.numeric(Elevation,1)),data=sp_persp[sp_persp$Form=='Sub-Fruticose',],pch=16)
points(Time_Absorption~jitter(as.numeric(Elevation)),data=sp_persp[sp_persp$Form=='Foliose',],pch=17,col=2)
points(Time_Absorption~jitter(as.numeric(Elevation)),data=sp_persp[sp_persp$Form=='Crustose',],pch=18,col=3)

xyplot(Time_Absorption~as.numeric(Elevation)|Form,data=sp_persp)

library(lattice)#########Microclimate Data##########
ibutt<-read.csv('iButton_combined.csv')
#clim<-read.csv('Microclim_Summary.csv')
ibutt<-aggregate(ibutt[,c(4,5,6,7)],by=list(Site=ibutt$Site,Transect=ibutt$Transect,Hour=ibutt$Hour,Day=ibutt$Day,Month=ibutt$Month,Year=ibutt$Year),mean,na.rm=T)
hourly_medians<-aggregate(ibutt[,c(4,5,6,7)],by=list(Site=ibutt$Site,Hour=ibutt$Hour,Transect=ibutt$Transect),median,na.rm=T)

##### Problem remains with the following aggregates, in that the early data is hourly and the later data is bi-hourly, leading to some asymmetries in values. Only affects Transect A.
saturated99<-aggregate(ibutt[ibutt$RH>99,2],by=list(Site=ibutt[ibutt$RH>99,]$Site,Transect=ibutt[ibutt$RH>99,]$Transect,Month=ibutt[ibutt$RH>99,]$Month),length)
saturated85<-aggregate(ibutt[ibutt$RH>85,2],by=list(Site=ibutt[ibutt$RH>85,]$Site,Transect=ibutt[ibutt$RH>85,]$Transect,Month=ibutt[ibutt$RH>85,]$Month),length)
names(saturated85)[4]<-'Sat85' #Relative humidity for activation of chlorolichens following Lange, Lakatos, etc
saturated0.01kPa<-aggregate(ibutt[ibutt$VPD<0.01,2],by=list(Site=ibutt[ibutt$VPD<0.01,]$Site,Transect=ibutt[ibutt$VPD<0.01,]$Transect,Month=ibutt[ibutt$VPD<0.01,]$Month),length)
names(saturated0.01kPa)[4]<-'Sat0.01kPa'
saturated0.1kPa<-aggregate(ibutt[ibutt$VPD<0.1,2],by=list(Site=ibutt[ibutt$VPD<0.1,]$Site,Transect=ibutt[ibutt$VPD<0.1,]$Transect,Month=ibutt[ibutt$VPD<0.1,]$Month),length)
names(saturated0.1kPa)[4]<-'Sat0.1kPa'
drystress1.0kPa<-aggregate(ibutt[ibutt$VPD>1.0,3],by=list(Site=ibutt[ibutt$VPD>1.0,]$Site,Transect=ibutt[ibutt$VPD>1.0,]$Transect,Month=ibutt[ibutt$VPD>1.0,]$Month),length)
names(drystress1.0kPa)[4]<-'Dry1kPa'
drystress2.0kPa<-aggregate(ibutt[ibutt$VPD>2.0,3],by=list(Site=ibutt[ibutt$VPD>2.0,]$Site,Transect=ibutt[ibutt$VPD>2.0,]$Transect,Month=ibutt[ibutt$VPD>2.0,]$Month),length)
names(drystress2.0kPa)[4]<-'Dry2kPa'
drystress3.0kPa<-aggregate(ibutt[ibutt$VPD>3.0,3],by=list(Site=ibutt[ibutt$VPD>3.0,]$Site,Transect=ibutt[ibutt$VPD>3.0,]$Transect,Month=ibutt[ibutt$VPD>3.0,]$Month),length)
names(drystress3.0kPa)[4]<-'Dry3kPa'
days<-aggregate(ibutt[,9],by=list(Site=ibutt$Site,Transect=ibutt$Transect,Year=ibutt$Year,Day=ibutt$Day),mean)
nday<-aggregate(days[,3],by=list(Site=days$Site,Transect=days$Transect),length) #number of days of data per logger (to get daily means)
#merge
sums<-merge(saturated0.1kPa,drystress1.0kPa,all=T)
sums<-merge(sums,drystress3.0kPa,all=T)
sums<-merge(sums,drystress2.0kPa,all=T)
sums<-merge(sums,saturated0.01kPa,all=T)
sums<-merge(sums,saturated85,all=T)
sums<-merge(sums,nday,all=T)
#and let's trim out all of the dataloggers from A that we have less than a year's worth of data:
sums<-sums[sums$x>365,]

par(mfrow=c(1,1))
plot(2*Dry3kPa*12/x~Site,data=sums[sums$Transect=='B',],ylim=c(0,24),pch=16,cex=2,yaxs='i')
points(2*Sat0.1kPa*12/x~Site,data=sums[sums$Transect=='B',],cex=2,pch=16,col=3)

plot(12*2*Dry3kPa/x~Site,data=sums[sums$Transect=='A',],ylim=c(0,15),pch=16,cex=2,yaxs='i')
points(12*2*Sat0.1kPa/x~Site,data=sums[sums$Transect=='A',],cex=2,pch=16,col=3)

plot(12*2*Sat0.1kPa/x~Site,data=sums[sums$Transect=='B',],cex=2,pch=16,col=3,type='b',ylim=c(0,24),yaxs='i')
lines((24-12*2*Dry1kPa/x)~Site,data=sums[sums$Transect=='B',],cex=2,pch=16,col=2,type='b')

plot(Dry3kPa*2*12/x~Month,data=sums[sums$Transect=='B',],pch=16,ylim=c(0,12),cex=2) #the seasonality is very strong!
points(Sat0.01kPa*24/x~Month,data=sums[sums$Transect=='B',],pch=16,col=3,cex=2)

plot(24-Dry1kPa*2*12/x~Month,data=sums[sums$Transect=='B'&sums$Site==500,],pch=16,ylim=c(0,24),cex=2,type='b',yaxs='i',xaxs='i',ylab='Duration of conditions') #the seasonality is very strong!
points(Sat85*2*12/x~Month,data=sums[sums$Transect=='B'&sums$Site==500,],pch=16,col=3,cex=2,type='b')
x_dry<-c(1:12,12,1)
y_dry<-c(24-sums[sums$Transect=='B'&sums$Site==500,]$Dry1kPa*2*12/sums[sums$Transect=='B'&sums$Site==500,]$x,24,24)
x_sat<-c(1:12,12,1)
y_sat<-c(sums[sums$Transect=='B'&sums$Site==500,]$Sat85*2*12/sums[sums$Transect=='B'&sums$Site==500,]$x,0,0)
polygon(x_dry,y_dry,angle=45,density=8)
text(3,22,'VPD > 1kPa',cex=2)
polygon(x_sat,y_sat,angle=-45,density=8,col=3)
text(7,3,'RH > 85%
     Chlorolichen activation',cex=2)

#This seasonality is drastic enough that we should separate out two key periods: dry (Dec-Feb) and foggy(May-Jul)
sums_dry<-sums[sums$Month==12|sums$Month<3,]
sums_fog<-sums[sums$Month<8&sums$Month>4,]

par(mfrow=c(2,1),mar=c(4,5,1,1),mgp=c(3,1,0))
plot(3*2*Dry2kPa/x~Site,data=sums_dry[sums_dry$Transect=='A',],ylim=c(0,6),pch=16,cex=2,yaxs='i',ylab='Mean Daily Hours',xlab='black: VPD > 2.0 kPa; green: VPD < 0.1 kPa',cex.lab=1.5,cex.axis=1.25)
points(3*2*Sat0.1kPa/x~Site,data=sums_dry[sums_dry$Transect=='A',],cex=2,pch=16,col=3)
text(550,5,'Dry Season',cex=2)
plot(3*2*Dry2kPa/x~Site,data=sums_fog[sums_fog$Transect=='A',],ylim=c(0,6),pch=16,cex=2,yaxs='i',ylab='Mean Daily Hours',xlab='Elevation (m)',cex.lab=1.5,cex.axis=1.25)
points(3*2*Sat0.1kPa/x~Site,data=sums_fog[sums_fog$Transect=='A',],cex=2,pch=16,col=3)
text(550,5,'Fog Season',cex=2)


###################Min/Max data#########
mins<-aggregate(ibutt[,c(7,8,10)],by=list(Site=ibutt[,]$Site,Day=ibutt[,]$Day,Transect=ibutt[,]$Transect,Year=ibutt[,]$Year),min)
maxs<-aggregate(ibutt[,c(7,8,10)],by=list(Site=ibutt[,]$Site,Day=ibutt[,]$Day,Transect=ibutt[,]$Transect,Year=ibutt[,]$Year),max)
meds<-aggregate(ibutt[,c(7,8,10)],by=list(Site=ibutt[,]$Site,Day=ibutt[,]$Day,Transect=ibutt[,]$Transect,Year=ibutt[,]$Year),median)

boxplot(T~Site,data=mins[mins$Transect=='A'&mins$Year==2016,],ylim=c(5,50),ylab='Min/Max Surface Temperature (°C)',xlab='Elevation (m.a.s.l.',notch=T,main='Transect A')
boxplot(T~Site,data=maxs[maxs$Transect=='A'&maxs$Year==2016,],add=T,notch=T)
boxplot(T~Site,data=mins[mins$Transect=='B',],ylim=c(5,50),ylab='Min/Max Surface Temperature (°C)',xlab='Elevation (m.a.s.l.',notch=T,main='Transect B')
boxplot(T~Site,data=maxs[maxs$Transect=='B',],add=T,notch=T)

boxplot(VPD~Site,data=mins[mins$Transect=='A'&mins$Year==2016,],ylim=c(0,6),ylab='Min/Max VPD (kPa)',xlab='Elevation (m.a.s.l.)',notch=T,main='Transect A')
boxplot(VPD~Site,data=maxs[maxs$Transect=='A'&mins$Year==2016,],add=T,notch=T)
boxplot(VPD~Site,data=mins[mins$Transect=='B',],ylim=c(0,10),ylab='Min/Max VPD (kPa)',xlab='Elevation (m.a.s.l.)',notch=T,main='Transect B')
boxplot(VPD~Site,data=maxs[maxs$Transect=='B',],add=T,notch=T)

summary(lm(T~Site,data=mins[mins$Year==2016,]))
summary(lm(T~Site*Transect,data=mins[mins$Year==2016,]))
summary(lm(T~Site,data=meds[meds$Year==2016,]))
summary(lm(T~Site*Transect,data=meds[meds$Year==2016,]))
summary(lm(T~Site,data=maxs[maxs$Year==2016,]))
summary(lm(T~Site*Transect,data=maxs[maxs$Year==2016,]))

summary(lm(VPD~Site,data=mins[mins$Year==2016,]))
summary(lm(VPD~Site*Transect,data=mins[mins$Year==2016,]))
########################


#####Anatomy#####
anat<-read.csv('Anatomia2.csv')


