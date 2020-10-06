####Overall code for Patache community data, combining various data forms into master files for environment, widely measured traits and community composition###
setwd("~/Dropbox/Lomas/Nat Geo/Combined Data")

##Load libraries
library(smatr)
library(lattice)
library(nlme)
library(labdsv)
library(vegan)
#########Emphasis on Transect B########

##########Taxa####
taxa<-read.csv('Taxa_TransectB.csv',stringsAsFactors=T)
taxa$Transect<-'B'
taxa$Elevation<-NA #this creates an empty column of data for storing the elevation records from the collection name
taxa$Quadrat<-NA
taxa$Number<-NA

for (i in seq(2,nrow(taxa),1)){ #We are using a for loop to go through the data file row by row, and parse the collection name into Elevation, Quadrat and Number identifiers to align with other datasets
    taxa[i,10]<-unlist(strsplit(as.character(taxa[i,3]),c("_")))[2] #for example here we remove the second element from the collection name (Elevation) and asign that value to the Elevation column of the data spreadsheet
  taxa[i,11]<-unlist(strsplit(as.character(taxa[i,3]),"_"))[3]
  taxa[i,12]<-unlist(strsplit(as.character(taxa[i,3]),"_"))[4]
}
head(taxa)
summary(taxa)

##############################
#################Cover Master#############
cover<-read.csv('CoverMaster.csv',stringsAsFactors=T)
summary(cover)
coverB<-cover[cover$Transect=='B',] #just transect B (richer trait data)


#############Contact Angle Data############
angles<-read.csv('ContactAngles_TransectB.csv',stringsAsFactors=T)
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


sunflowerplot(Time_Absorption~Elevation,data=cov_ang[,])
sunflowerplot(Time_Absorption~Elevation,data=cov_ang[!cov_ang$Cover=='T',])

sp_trait<-read.csv('Spp.csv',stringsAsFactors=T)

cov_trait<-merge(cov_ang,sp_trait)


##################
isos<-read.csv('Isotopes_TransectB.csv',stringsAsFactors=T)
summary(isos)
iso_tax<-merge(isos,taxa[,c(1,2,3,8:12)],all=T)

cov_trait<-merge(cov_trait,iso_tax,all=T)

summary(lm(d13C~as.numeric(Elevation)*Time_Absorption,cov_trait[,]))

##############Spectral features as traits (from PatacheSpectra.R output)
specs<-read.csv('Spectral_Features.csv',stringsAsFactors=T)
spec_trait<-merge(cov_trait,specs,all=T)


########Quadrat Summary Info ###### roughness, aspect, etc
quad<-read.csv('Quadrat_Master.csv',stringsAsFactors=T)

quad_cov<-merge(quad,cov_ang,all=T)

#################
#########################Species perspective######################
spp_summ<-aggregate(coverB[,c(1)],coverB[,c(1,2,4,6)],length) #number of quadrats for each species/elevation
names(spp_summ)[5]<-'N_Quad_Pres'


sp_persp<-merge(spec_trait,spp_summ)
sp_persp<-merge(sp_persp,sp_trait)


plot(Time_Absorption~jitter(as.numeric(Elevation,1)),data=sp_persp[sp_persp$Form=='Sub-Fruticose',],pch=16)
points(Time_Absorption~jitter(as.numeric(Elevation)),data=sp_persp[sp_persp$Form=='Foliose',],pch=17,col=2)
points(Time_Absorption~jitter(as.numeric(Elevation)),data=sp_persp[sp_persp$Form=='Crustose',],pch=18,col=3)

xyplot(Time_Absorption~as.numeric(Elevation)|Form,data=sp_persp)
sunflowerplot(Time_Absorption~as.numeric(Elevation),data=sp_persp[sp_persp$Photobiont=='Trebouxia',],xlab='Elevation (m above sea level',ylab='Time to absorption (seconds)',main='Trebouxioid lichens')
sunflowerplot(Time_Absorption~as.numeric(Elevation),data=sp_persp[sp_persp$Photobiont=='Trentepohlia',],xlab='Elevation (m above sea level',ylab='Time to absorption (seconds)',main='Trentepohlioid lichens')



library(lattice)#########Microclimate Data##########
ibutt<-read.csv('iButton_combined.csv',stringsAsFactors=T)
#clim<-read.csv('Microclim_Summary.csv')
ibutt<-aggregate(ibutt[,c(4,5,6,7)],by=list(Site=ibutt$Site,Transect=ibutt$Transect,Hour=ibutt$Hour,Day=ibutt$Day,Month=ibutt$Month,Year=ibutt$Year),mean,na.rm=T)
ibutt<-na.omit(ibutt)
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

###Round 2 of variables: aim for fits that better reflect the physiologies of Trebouxioid and Trentepohlioid, focusing on VPD rather than RH. These can be the approx VPD threshold for trebouxioid activation (0.3 kPa, based on 85% at 20ºC) and a conservative guess at one for trentepohlioids (0.7 kPa, based on 65% at 20ºC)
saturated0.3kPa<-aggregate(ibutt[ibutt$VPD<0.3,2],by=list(Site=ibutt[ibutt$VPD<0.3,]$Site,Transect=ibutt[ibutt$VPD<0.3,]$Transect,Month=ibutt[ibutt$VPD<0.3,]$Month),length)
names(saturated0.3kPa)[4]<-'Sat0.3kPa'
saturated0.7kPa<-aggregate(ibutt[ibutt$VPD<0.7,2],by=list(Site=ibutt[ibutt$VPD<0.7,]$Site,Transect=ibutt[ibutt$VPD<0.7,]$Transect,Month=ibutt[ibutt$VPD<0.7,]$Month),length)
names(saturated0.7kPa)[4]<-'Sat0.7kPa'

days<-aggregate(ibutt[,9],by=list(Site=ibutt$Site,Transect=ibutt$Transect,Year=ibutt$Year,Day=ibutt$Day),mean,na.rm=T)
nday<-aggregate(days[,3],by=list(Site=days$Site,Transect=days$Transect),length) #number of days of data per logger (to get daily means)
#merge
sums<-merge(saturated0.1kPa,drystress1.0kPa,all=T)
sums<-merge(sums,drystress3.0kPa,all=T)
sums<-merge(sums,drystress2.0kPa,all=T)
sums<-merge(sums,saturated0.01kPa,all=T)
sums<-merge(sums,saturated85,all=T)
sums<-merge(sums,saturated0.3kPa,all=T)
sums<-merge(sums,saturated0.7kPa,all=T)
sums<-merge(sums,nday,all=T)
#and let's trim out all of the dataloggers from A that we have less than a year's worth of data:
#sums<-sums[sums$x>365,]

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
points(Sat0.3kPa*2*12/x~Month,data=sums[sums$Transect=='B'&sums$Site==500,],pch=16,col=3,cex=2,type='b')
x_dry<-c(1:12,12,1)
y_dry<-c(24-sums[sums$Transect=='B'&sums$Site==500,]$Dry1kPa*2*12/sums[sums$Transect=='B'&sums$Site==500,]$x,24,24)
x_sat<-c(1:12,12,1)
y_sat<-c(sums[sums$Transect=='B'&sums$Site==500,]$Sat0.3kPa*2*12/sums[sums$Transect=='B'&sums$Site==500,]$x,0,0)
polygon(x_dry,y_dry,angle=45,density=8)
text(3,22,'VPD > 1kPa',cex=2)
polygon(x_sat,y_sat,angle=-45,density=8,col=3)
text(7,3,'RH > 85%
     Chlorolichen activation',cex=2)

par(mfrow=c(3,1))
plot(24-Dry1kPa*2*12/x~Month,data=sums[sums$Transect=='B'&sums$Site==800,],pch=16,ylim=c(0,24),cex=2,type='b',yaxs='i',xaxs='i',ylab='Duration of conditions',xlab='',cex.lab=1.5,cex.axis=1.5) #the seasonality is very strong!
points(Sat0.1kPa*2*12/x~Month,data=sums[sums$Transect=='B'&sums$Site==800,],pch=16,col=3,cex=2,type='b')
x_dry<-c(1:12,12,1)
y_dry<-c(24-sums[sums$Transect=='B'&sums$Site==800,]$Dry1kPa*2*12/sums[sums$Transect=='B'&sums$Site==800,]$x,24,24)
x_sat<-c(1:12,12,1)
y_sat<-c(sums[sums$Transect=='B'&sums$Site==800,]$Sat0.1kPa*2*12/sums[sums$Transect=='B'&sums$Site==800,]$x,0,0)
polygon(x_dry,y_dry,angle=45,density=8)
text(3,22,'Rapid Drying',cex=2)
polygon(x_sat,y_sat,angle=-45,density=8,col=3)
text(7,3,'Near Saturation',cex=2)

plot(24-Dry1kPa*2*12/x~Month,data=sums[sums$Transect=='B'&sums$Site==500,],pch=16,ylim=c(0,24),cex=2,type='b',yaxs='i',xaxs='i',ylab='Duration of conditions',xlab='',cex.lab=1.5,cex.axis=1.5) #the seasonality is very strong!
points(Sat0.1kPa*2*12/x~Month,data=sums[sums$Transect=='B'&sums$Site==500,],pch=16,col=3,cex=2,type='b')
x_dry<-c(1:12,12,1)
y_dry<-c(24-sums[sums$Transect=='B'&sums$Site==500,]$Dry1kPa*2*12/sums[sums$Transect=='B'&sums$Site==500,]$x,24,24)
x_sat<-c(1:12,12,1)
y_sat<-c(sums[sums$Transect=='B'&sums$Site==500,]$Sat0.1kPa*2*12/sums[sums$Transect=='B'&sums$Site==500,]$x,0,0)
polygon(x_dry,y_dry,angle=45,density=8)
text(3,22,'Rapid Drying',cex=2)
polygon(x_sat,y_sat,angle=-45,density=8,col=3)
text(7,3,'Near Saturation',cex=2)

plot(24-Dry1kPa*2*12/x~Month,data=sums[sums$Transect=='B'&sums$Site==300,],pch=16,ylim=c(0,24),cex=2,type='b',yaxs='i',xaxs='i',ylab='Duration of conditions',cex.lab=1.5,cex.axis=1.5) #the seasonality is very strong!
points(Sat0.1kPa*2*12/x~Month,data=sums[sums$Transect=='B'&sums$Site==300,],pch=16,col=3,cex=2,type='b')
x_dry<-c(1:12,12,1)
y_dry<-c(24-sums[sums$Transect=='B'&sums$Site==300,]$Dry1kPa*2*12/sums[sums$Transect=='B'&sums$Site==300,]$x,24,24)
x_sat<-c(1:12,12,1)
y_sat<-c(sums[sums$Transect=='B'&sums$Site==300,]$Sat0.1kPa*2*12/sums[sums$Transect=='B'&sums$Site==300,]$x,0,0)
polygon(x_dry,y_dry,angle=45,density=8)
text(3,22,'Rapid Drying',cex=2)
polygon(x_sat,y_sat,angle=-45,density=8,col=3)
text(7,3,'Near Saturation',cex=2)

####
#This seasonality is drastic enough that we should separate out two key periods: dry (Dec-Feb) and foggy(May-Jul)
sums_dry<-sums[sums$Month==12|sums$Month<3,]
sums_fog<-sums[sums$Month<8&sums$Month>4,]

par(mfrow=c(2,1),mar=c(4,5,1,1),mgp=c(3,1,0))
plot(2*12*Dry2kPa/x~Site,data=sums_dry[sums_dry$Transect=='A',],ylim=c(0,6),pch=16,cex=2,yaxs='i',ylab='Mean Daily Hours',xlab='black: VPD > 2.0 kPa; green: VPD < 0.1 kPa',cex.lab=1.5,cex.axis=1.25)
points(2*12*Sat0.1kPa/x~Site,data=sums_dry[sums_dry$Transect=='A',],cex=2,pch=16,col=3)
text(550,5,'Dry Season',cex=2)
plot(2*12*Dry2kPa/x~Site,data=sums_fog[sums_fog$Transect=='A',],ylim=c(0,6),pch=16,cex=2,yaxs='i',ylab='Mean Daily Hours',xlab='Elevation (m)',cex.lab=1.5,cex.axis=1.25)
points(2*12*Sat0.1kPa/x~Site,data=sums_fog[sums_fog$Transect=='A',],cex=2,pch=16,col=3)
text(550,5,'Fog Season',cex=2)

############Aggregate all of these microclimate results into site descriptors for the multivariate analyses#############33

agg_dry<-aggregate(sums_dry[,c(4:12)],by=list(Site=sums_dry$Site,Transect=sums_dry$Transect),sum,na.rm=T)
agg_dry$Sat85_Dry<-24*agg_dry$Sat85/agg_dry$x
agg_dry$Sat0.1kPa_Dry<-24*agg_dry$Sat0.1kPa/agg_dry$x
agg_dry$Dry1kPa_Dry<-24*agg_dry$Dry1kPa/agg_dry$x
agg_dry$Dry3kPa_Dry<-24*agg_dry$Dry3kPa/agg_dry$x
agg_dry$Dry0.3kPa_Dry<-24*agg_dry$Sat0.3kPa/agg_dry$x
agg_dry$Dry0.7kPa_Dry<-24*agg_dry$Sat0.7kPa/agg_dry$x
agg_fog<-aggregate(sums_fog[,c(4:12)],by=list(Site=sums_fog$Site,Transect=sums_fog$Transect),sum,na.rm=T)
agg_fog$Sat85_Fog<-24*agg_fog$Sat85/agg_fog$x
agg_fog$Sat0.1kPa_Fog<-24*agg_fog$Sat0.1kPa/agg_fog$x
agg_fog$Dry1kPa_Fog<-24*agg_fog$Dry1kPa/agg_fog$x
agg_fog$Dry3kPa_Fog<-24*agg_fog$Dry3kPa/agg_fog$x
agg_fog$Sat0.3kPa_Fog<-24*agg_fog$Sat0.3kPa/agg_fog$x
agg_fog$Sat0.7kPa_Fog<-24*agg_fog$Sat0.7kPa/agg_fog$x

combined<-merge(agg_dry[,c(1,2,12:17)],agg_fog[,c(1,2,12:17)],by=c('Site','Transect'),all=T)

#####Alternative perpective: daylight vs nocturnal hydration periods#####
#we'll only bother
#classify data by night/day
ibutt$DayNight<-as.factor('Day')
levels(ibutt$DayNight)<-c('Day','Night')
for (i in seq(1,nrow(ibutt),1)){ 
  if (ibutt[i,3]>19|ibutt[i,3]<7) ibutt[i,11]<-factor('Night')
}
summary(ibutt)
sat85<-aggregate(ibutt[ibutt$RH>85,2],by=list(Site=ibutt[ibutt$RH>85,]$Site,Transect=ibutt[ibutt$RH>85,]$Transect,Month=ibutt[ibutt$RH>85,]$Month,DayNight=ibutt[ibutt$RH>85,]$DayNight),length)
names(sat85)[5]<-'Sat85' #Relative humidity for activation of chlorolichens following Lange, Lakatos, etc
plot(Sat85~Site,data=sat85[sat85$Transect=='B',])
points(Sat85~Site,data=sat85[sat85$Transect=='B'&sat85$DayNight=='Night',],pch=16,col=3)
points(Sat85~Site,data=sat85[sat85$Transect=='B'&sat85$DayNight=='Day',],pch=16,col=2)
summary(lm(Sat85~DayNight*Site,data=sat85[sat85$Transect=='B',]))
###################Min/Max data#########

mins<-aggregate(ibutt[,c(7,8,10)],by=list(Site=ibutt[,]$Site,Day=ibutt[,]$Day,Transect=ibutt[,]$Transect,Year=ibutt[,]$Year),min,na.rm=T)
maxs<-aggregate(ibutt[,c(7,8,10)],by=list(Site=ibutt[,]$Site,Day=ibutt[,]$Day,Transect=ibutt[,]$Transect,Year=ibutt[,]$Year),max,na.rm=T)
meds<-aggregate(ibutt[,c(7,8,10)],by=list(Site=ibutt[,]$Site,Day=ibutt[,]$Day,Transect=ibutt[,]$Transect,Year=ibutt[,]$Year),median,na.rm=T)


min<-aggregate(mins[,c(5,6,7)],by=list(Site=mins[,]$Site,Transect=mins[,]$Transect),mean)
names(min)[c(3,4,5)]<-c('T_min','RH_min','VPD_min')
max<-aggregate(maxs[,c(5,6,7)],by=list(Site=maxs[,]$Site,Transect=maxs[,]$Transect),mean)
names(max)[c(3,4,5)]<-c('T_max','RH_max','VPD_max')
med<-aggregate(meds[,c(5,6,7)],by=list(Site=meds[,]$Site,Transect=meds[,]$Transect),mean)
names(med)[c(3,4,5)]<-c('T_med','RH_med','VPD_med')

min_max<-merge(min,max)
med_max<-merge(min_max,med)

combined<-merge(combined,med_max,all=T)
names(combined)[1]<-'Elevation'
#write.csv(combined,'Summarized_Microclimate.csv')

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

#Combine this in with the species traits:
sp_clim<-merge(sp,combined,all=T)

summary(lm(d13C~Photobiont*as.numeric(Elevation),sp_clim))
#####Anatomy#####
anat<-read.csv('Anatomia2.csv',stringsAsFactors=T)

#####
comb_B<-combined[combined$Transect=='B',]
#since we lack data for 600m, we can try to interpolate it. Let's start with the crude version, the mean of 500 and 700. (Have done this manually and externally instead)
microclim_B<-read.csv('Microclim_Interp_B_new.csv',stringsAsFactors=T)
quadB<-quad[quad$Transect=='B',]

env_B<-merge(quadB[,c(2,3,4,5,6,12)],microclim_B,all=T)
#write.csv(env_B,'Env_matrix_B.csv')


##Isotope data was measured for both quadrat vouchers and intra-specific variation collections. Since these were formatted differently, they need to be brought back together if one is to examine the full data set for those. This has been done manually:
iso_only<-read.csv('allisos.csv',stringsAsFactors=T)
summary(iso_only)

plot(d13C~d15N,data=iso_only[iso_only$Genus=='Roccellina',],ylab=expression(paste(delta^13,'C')),xlab=expression(paste(delta^15,'N')),cex=1.5,cex.lab=1.5,cex.axis=1.25)
points(d13C~d15N,data=iso_only[iso_only$Genus=='Roccellina'&iso_only$Elevation=='300',],pch=16,col='darkred',cex=1.5)
points(d13C~d15N,data=iso_only[iso_only$Genus=='Roccellina'&iso_only$Elevation=='400',],pch=16,col='red',cex=1.5)
points(d13C~d15N,data=iso_only[iso_only$Genus=='Roccellina'&iso_only$Elevation=='500',],pch=16,col='orange',cex=1.5)
points(d13C~d15N,data=iso_only[iso_only$Genus=='Roccellina'&iso_only$Elevation=='600',],pch=16,col='yellow',cex=1.5)
points(d13C~d15N,data=iso_only[iso_only$Genus=='Roccellina'&iso_only$Elevation=='700',],pch=16,col='light blue',cex=1.5)
points(d13C~d15N,data=iso_only[iso_only$Genus=='Roccellina'&iso_only$Elevation=='800',],pch=16,col='blue',cex=1.5)



###############Trait summaries at the quadrat level:
#Version 1: unweighted by cover
QTrait_UnW<-aggregate(spec_trait[,c(8,9,10,21,22,23,24,27:138)],by=list(Transect=spec_trait$Transect,Elevation=spec_trait$Elevation,Quadrat=spec_trait$Quadrat),FUN=mean,na.rm=TRUE)
#lets add in photobiont proportions&frequencies as well:
QTrait_UnW_Treb_Count<-aggregate(spec_trait[spec_trait$Photobiont=='Trebouxia',7],by=list(Transect=spec_trait[spec_trait$Photobiont=='Trebouxia',]$Transect,Elevation=spec_trait[spec_trait$Photobiont=='Trebouxia',]$Elevation,Quadrat=spec_trait[spec_trait$Photobiont=='Trebouxia',]$Quadrat),FUN=length)
QTrait_UnW_Trent_Count<-aggregate(spec_trait[spec_trait$Photobiont=='Trentepohlia',7],by=list(Transect=spec_trait[spec_trait$Photobiont=='Trentepohlia',]$Transect,Elevation=spec_trait[spec_trait$Photobiont=='Trentepohlia',]$Elevation,Quadrat=spec_trait[spec_trait$Photobiont=='Trentepohlia',]$Quadrat),FUN=length)
names(QTrait_UnW_Treb_Count)[4]<-'NSpp_Treb'
names(QTrait_UnW_Trent_Count)[4]<-'NSpp_Trent'
QTrait_UnW_2<-merge(QTrait_UnW,QTrait_UnW_Treb_Count,all=T)
QTrait_UnW_2<-merge(QTrait_UnW_2,QTrait_UnW_Trent_Count,all=T)
QTrait_UnW_2$Prop_Trent<-(QTrait_UnW_2$NSpp_Trent)/(QTrait_UnW_2$NSpp_Treb+QTrait_UnW_2$NSpp_Trent)

#write.csv(QTrait_UnW_2,'QuadratMeanTraits_Unweighted.csv')

##Version 2: weighted by cover
#This obviously takes a few more steps. We should come up with a way of treating Trace values
#The simple version is to eliminate Trace abundance species:
cov_trait_noT<-spec_trait[!spec_trait$Cover=='T',]
cov_trait_noT$Cover<-as.numeric(as.character(cov_trait_noT$Cover)) #then the column can be read as numeric
summary(cov_trait_noT)
#and then to create some weighted columns:
#cov_trait_noT$Time_Surface_W<-cov_trait_noT$Time_Surface*cov_trait_noT$Cover
#kept that first one in as an example, but since adding in spectral indices, need to automate this process rather spelling it out at each step:
cov_trait_noT[,c(8,9,10,21:24,27:138)]<-mutate(cov_trait_noT[,c(8,9,10,21:24,27:138)],cov_trait_noT[,c(8,9,10,21:24,27:138)]*cov_trait_noT$Cover)
summary(cov_trait_noT)

#And now we can aggregate those weighted columns by summing them. Note that since the cover doesn't sum to 100, we need to account for that as well. At this step we will also sum all cover, and then use that to reweight the values
QTrait_W_noT<-aggregate(cov_trait_noT[,c(7,8,9,10,21:24,27:138)],by=list(Transect=cov_trait_noT$Transect,Elevation=cov_trait_noT$Elevation,Quadrat=cov_trait_noT$Quadrat),FUN=sum,na.rm=TRUE)
QTrait_W_noT[,c(5:123)]<-QTrait_W_noT[,c(5:123)]/QTrait_W_noT$Cover
write.csv(QTrait_UnW,'QuadratMeanTraits_Weighted_NoTrace2.csv')



#Figures:
#sp_persp is still quadrat level rather than transect level summary:
sp_trans<-aggregate(sp_persp[,c(15,16,17,21,22,23,24,27:139)],by=list(Form=sp_persp$Form,Photobiont=sp_persp$Photobiont,Coll=sp_persp$Collection2,Family=sp_persp$Family,Order=sp_persp$Order,Species=sp_persp$Species,Elevation=sp_persp$Elevation),mean)
sp_trans$Elevation<-as.numeric(sp_trans$Elevation)
par(mfrow=c(1,1))
plot(Time_Surface~jitter(Elevation,factor=0.5),data=sp_trans,ylab='Time to start of uptake (s)',xlab='Elevation (m)',cex.lab=1.5,cex=1.25,cex.axis=1.25,pch='')
points(Time_Surface~jitter(Elevation),data=sp_trans[sp_trans$Photobiont=='Trebouxia',],col='darkgreen',pch=16)
points(Time_Surface~jitter(Elevation),data=sp_trans[sp_trans$Photobiont=='Trentepohlia',],col='orange',pch=16)

plot(N_Quad_Pres~jitter(Elevation),data=sp_trans[sp_trans$Photobiont=='Trebouxia',],col='darkgreen',pch=16,ylim=c(0,10))
points(N_Quad_Pres~jitter(Elevation),data=sp_trans[sp_trans$Photobiont=='Trentepohlia',],col='orange',pch=16)
  
plot(NDVI~Elevation,data=sp_trans[,])
points(NDVI~Elevation,data=sp_trans[sp_trans$Photobiont=='Trentepohlia',],col='orange',pch=16)
points(NDVI~Elevation,data=sp_trans[sp_trans$Photobiont=='Trebouxia',],col='darkgreen',pch=16)
