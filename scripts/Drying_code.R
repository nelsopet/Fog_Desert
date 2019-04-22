###Process the drying curves from Patache lichens
library(lattice)
#removed setting working directory so paths can be relative for each file read in. Not sure this is the best
#way to do it but allows a user to pull the repo and not have to change working directory
#setwd("~/Dropbox/Lomas/Nat Geo/Combined Data")
dry<-read.csv('data/Drying_Patache_1.csv',na.string='na')
summary(dry)
dry$RWC<-(dry$Mass-dry$DryMass)/(dry$DryMass-dry$EmptyMass)
dry$DryNet<-dry$DryMass-dry$EmptyMass
xyplot(RWC~Elevation|Species,data=dry[dry$Time==0&dry$Treatment=='dry',])

plot(RWC~I(DryMass-EmptyMass),data=dry[dry$Time==0&dry$Treatment=='dry',])#Seems to be an effect of thallus mass, but inverse what I might have expected
xyplot(RWC~I(DryMass-EmptyMass)|Species,data=dry[dry$Time==0&dry$Treatment=='dry',])

library(plyr)
dry_rate<-function(df){ #here we create the function that we want to apply to each unique gridcell
  lm(RWC~log10(Time+1),data=df,na.action=na.omit) #in this case, just a linear model of area by time since plating, which for the moment seems to be ok (see the plot above), but may not hold throughout the experiment
}
models<-dlply(na.exclude(dry), .(Species,Elevation,Transect,Replicate),dry_rate) #this line applies the grow_rate function to every unique combination of Species/Individual/Plate/Cell, ie to each colony separately
#We don't want all of the model output at this point, so we want to trim what we get to the intercept, slope and r-squared value of each one:
rsq <- function(x) summary(x)$r.squared
rates<-ldply(models,function(x) c(coef(x), rsquare=rsq(x)))
names(rates)[5:6]=c('intercept','slope')
summary(rates)
plot(slope~Species,data=rates[!rates$Species=='Blank',])
xyplot(slope~Elevation|Species,data=rates[!rates$Species=='Blank',])

##Bring in additional drying rates for Rocc and Het
#removed setting working directory so paths can be relative for each file read in. Not sure this is the best
#way to do it but allows a user to pull the repo and not have to change working directory
#setwd("~/Dropbox/Lomas/Nat Geo/Data 2016/Drying")
#PRN: Where is RoccHetB_Try1.csv?
dry<-read.csv('data/RoccHetB_Try1.csv')
summary(dry)
dry$AdjMass<-(dry$Mass-dry$Paper)/dry$DryMass
dry$RWC<-(dry$Mass-dry$DryMass)/(dry$DryMass)

plot(AdjMass~Timepoint,data=dry[dry$Treatment=='Drying',])
plot(AdjMass~log(Timepoint),data=dry[dry$Treatment=='Drying',])


dry_rate<-function(df){ #here we create the function that we want to apply to each unique gridcell
  lm(RWC~log10(Timepoint+1),data=df,na.action=na.omit) #Calculate the 1st order decay constant (slope of log-time rate)
}
models<-dlply(na.exclude(dry[dry$Treatment=='Drying',]), .(Species,Transect,Elevation,Replicate),dry_rate)
rsq <- function(x) summary(x)$r.squared
rates2<-ldply(models,function(x) c(coef(x), rsquare=rsq(x)))
names(rates2)[5:6]=c('intercept','slope')
#rates<-rates[rates$rsquare>0.8,]#remove the cases where the r-square fit of the growth model is low (<0.8)
rates2<-unique(merge(rates2,dry[,c(1,2,3,4,7)])) #Bind on Dry Mass, since that likely also affects drying rate
summary(rates2)
summary(lm(slope~Species*DryMass*Elevation,data=rates2))

dry_rates<-merge(rates,rates2,all=T)

summary(lm(slope~Species*Elevation,data=dry_rates[!dry_rates$Species=='Control',]))
xyplot(slope~Elevation|Species,data=dry_rates,xlim=c(200,900))

summary(lm(slope~Species*Elevation,data=dry_rates[!dry_rates$Species=='Control'&dry_rates$Transect=='B',]))

rate_Rocc<-dry_rates[dry_rates$Species=='Roccellina',]

#Another set of Roccellina curves
#removed setting working directory so paths can be relative for each file read in. Not sure this is the best
#way to do it but allows a user to pull the repo and not have to change working directory
#setwd("~/Dropbox/Lomas/Nat Geo/Combined Data")
dry<-read.csv('data/Rocc_Dry_Nov18.csv',na.string='na')

summary(dry)


dry$RWC<-(dry$AdjMass-dry$DryMass)/(dry$DryMass)

plot(AdjMass~Time,data=dry[dry$Treatment=='dry',])
plot(AdjMass~log(Time),data=dry[dry$Treatment=='dry',])


dry_rate<-function(df){ #here we create the function that we want to apply to each unique gridcell
  lm(RWC~log10(Time+1),data=df,na.action=na.omit) #Calculate the 1st order decay constant (slope of log-time rate)
}
models<-dlply(na.exclude(dry[dry$Treatment=='dry',c(1:10,12)]), .(Species,Transect,Elevation,Replicate),dry_rate)
rsq <- function(x) summary(x)$r.squared
rates3<-ldply(models,function(x) c(coef(x), rsquare=rsq(x)))
names(rates3)[5:6]=c('intercept','slope')
rates3<-rates3[rates3$rsquare>0.8,]#remove the cases where the r-square fit of the growth model is low (<0.8)
rates3<-unique(merge(rates3,dry[,c(1,2,3,4,10)])) #Bind on Dry Mass, since that likely also affects drying rate
summary(rates3)
summary(lm(slope~Elevation,data=rates3))

dry_rates<-merge(dry_rates,rates3,all=T)

summary(lm(slope~Species*Elevation,data=dry_rates[!dry_rates$Species=='Control',]))
xyplot(slope~Elevation|Species,data=dry_rates,xlim=c(200,900))

summary(lm(slope~Species*Elevation,data=dry_rates[!dry_rates$Species=='Control'&dry_rates$Transect=='B',]))

rate_Rocc<-dry_rates[dry_rates$Species=='Roccellina',]
summary(lm(slope~Elevation+I(Elevation^2),data=rate_Rocc[rate_Rocc$Transect=='B',]))
plot(slope~Elevation,data=rate_Rocc[rate_Rocc$Transect=='B',])


#Area and mass contrasts (WHC protocol vs Full Soaking)
rocc_WHC<-dry[dry$Treatment=='blot',c(1,2,3,4,9,10,11,12)]
names(rocc_WHC)[c(5,7,8)]<-c('Mass_WHC','Area_WHC','RWC_WHC')
rocc_soak<-dry[dry$Treatment=='full_wet',c(1,2,3,4,9,10,11,12)]
names(rocc_soak)[c(5,7,8)]<-c('Mass_max','Area_max','RWC_max')

rocc<-merge(rocc_WHC,rocc_soak,all=T)
summary(lm(Mass_WHC~Mass_max,data=rocc))
