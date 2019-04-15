####################

setwd("~/Dropbox/Lomas/Nat Geo/Data 2016")
library(labdsv)
library(lattice)
library(polynom)

######1 read in data files, first looks, combine them into master dataset:
rou<-read.csv('Roughnesses.csv') ##Roughness, still needs to be compiled and analysed
summary(rou)
rouVar<-aggregate(rou[,7],by=list(Elevation=rou$Elevation,Quadrat=rou$Quadrat,Transect=rou$Transect),sd)
names(rouVar)<-c('Elevation','Quadrat','Transect','Rou_SD')

cov<-read.csv('CoverMaster.csv') ##Cover data from each quadrat, to be converted into multivariates, but also richness, etc 
summary(cov)
cov_summ<-aggregate(cov[,c(1)],cov[,c(1,2,3)],FUN=length)#Number of species for each quadrat
names(cov_summ)<-c('Elevation','Transect','Quadrat','NSpp')
plot(NSpp~Elevation,data=cov_summ)
boxplot(NSpp~Elevation,data=cov_summ[cov_summ$Transect=='A',],ylab='Quadrat Species Richness',yaxs='r',notch=T,xlab='Elevation(m.a.s.l.)')
boxplot(NSpp~Elevation,data=cov_summ[cov_summ$Transect=='B',],ylab='Quadrat Species Richness',yaxs='r',notch=T,xlab='Elevation(m.a.s.l.)')
spp_summ<-aggregate(cov[,c(1)],cov[,c(1,4)],length) #number of quadrats for each species/elevation

div<-unique(cov[,c(1,4)])
plot(summary(as.factor(div[,1])),pch=16,cex=2.5,ylab='# of Morphospecies',ylim=c(0,40))

spp<-unique(cov$Species)
spinf<-read.csv('Spp.csv')
covs<-merge(cov,spinf,all=T)
covs_summ<-aggregate(covs[,c(1)],covs[,c(1,2,3,13)],FUN=length)#Number of species for each quadrat
plot(x~Elevation,data=covs_summ[covs_summ$Form=='Crustose',],xlim=c(300,850))

trans<-unique(covs[,c(1,2,3)])
transect<-aggregate(trans[,c(2)],trans[,c(2,3)],FUN=length) #Number of species per elevation
plot(x~Elevation,data=transect[transect$Transect=='A',],cex=2,pch=16,cex.lab=2,cex.axis=2,ylab='#Species in transect',xlim=c(300,850),xlab='Elevation (m.a.s.l.)',ylim=c(0,50))
points(x~Elevation,data=transect[transect$Transect=='B',],cex=2,pch=16,col=2)

quad<-unique(covs[,c(2,3,4)])
quads<-aggregate(quad[,3],quad[,c(1,2)],FUN=length) #Counts the number of quadrats for each elevation/transect
names(quads)[3]<-'NQuad'

spp_summ<-aggregate(covs[,c(2)],covs[,c(1,2,3,11,13)],length) #number of quadrats for each species/elevation

gener<-unique(covs[,c(2,3,4,11)])
gen_summ<-aggregate(gener[,c(4)],gener[,c(1,2,4)],length) #number of quadrats for each genus/elevation Interpretation is complicated by the varying number of quadrats/elevation, so need to merge in that info....
gen_sum<-merge(gen_summ,quads,all=T)
summary(lm(x/NQuad~Elevation,data=gen_sum[gen_sum$Genus=='Caloplaca',])) #Contrast an elevation independent genus to...
summary(lm(x/NQuad~Elevation,data=gen_sum[gen_sum$Genus=='Niebla',])) #A highly elevation dependent one
######Generic form to batch produce these####
par(mfrow=c(2,4))
taxa<-c('Buellia','Rhizocarpon','Chrysothrix','Caloplaca','Roccellina','Heterodermia','Niebla','Pentagenella')
for (i in c(1:8)){
species<-taxa[i]
#species<-'Buellia'
lm_sp<-lm(x/NQuad~Elevation+Transect,data=gen_sum[gen_sum$Genus==species,])
summary(lm_sp)
plot(x/NQuad~Elevation,data=gen_sum[gen_sum$Genus==species,],cex=2,pch=16,cex.lab=2,cex.axis=2,ylab='Proportion of quadrats present',xlim=c(300,850),ylim=c(0,1))
points(x/NQuad~Elevation,data=gen_sum[gen_sum$Genus==species&gen_sum$Transect=='B',],pch=16,col=2,cex=2)
text(500,0.99,species,cex=2,font=2)
text(500,0.9,'r2 = ',cex=2)
text(560,0.9,round(summary(lm_sp)$adj.r.squared,digits=2),cex=2)
abline(a=as.numeric(coef(lm_sp)[1]),b=as.numeric(coef(lm_sp)[2]),lty=2)
}

###################### ##############
#text(500,0.8,'p > ',cex=2)
#text(560,0.8,round(as.numeric(summary(lm_sp)[9]),digits=2),cex=2)

forms<-unique(covs[,c(2,3,4,13)])
form_summ<-aggregate(forms[,c(4)],forms[,c(1,2,4)],length) #number of quadrats for each growthform/elevation. Interpretation is complicated by the varying number of quadrats/elevation, so need to merge in that info....
form_sum<-merge(form_summ,quads,all=T)
xyplot(x~Elevation|Form,data=form_sum)
plot(100*x/NQuad~Elevation,data=form_sum[form_sum$Form=='Sub-Fruticose'&form_sum$Transect=='A',],ylab='Percentage Quadrats Occupied',ylim=c(0,105),yaxs='i',xlim=c(250,850),pch=16,cex=2,cex.axis=2,cex.lab=2,type='b')
lines(100*x/NQuad~Elevation,data=form_sum[form_sum$Form=='Crustose'&form_sum$Transect=='A',],pch=16,col=2,cex=2,type='b')
lines(100*x/NQuad~Elevation,data=form_sum[form_sum$Form=='Fruticose'&form_sum$Transect=='A',],pch=16,col=3,cex=2,type='b')
lines(100*x/NQuad~Elevation,data=form_sum[form_sum$Form=='Foliose'&form_sum$Transect=='A',],pch=16,col=4,cex=2,type='b')
lines(100*x/NQuad~Elevation,data=form_sum[form_sum$Form=='Leprose'&form_sum$Transect=='A',],pch=16,col=5,cex=2,type='b')
points(250,100,cex=2,col=2,pch=16)
text(350,100,'Crustose',cex=2)
points(250,92,cex=2,col=5,pch=16)
text(350,92,'Leprose',cex=2)
points(250,84,cex=2,col=1,pch=16)
text(360,84,'Sub-fruticose',cex=2)
points(250,76,cex=2,col=4,pch=16)
text(350,76,'Foliose',cex=2)
points(250,68,cex=2,col=3,pch=16)
text(350,68,'Fruticose',cex=2)

sigcov<-droplevels(covs[!covs$Cover=='T'&!covs$Cover=='',]) #subset with just cover greater than trace
sigcov$Cover<-as.numeric(levels(sigcov$Cover))[sigcov$Cover] #needed to interpret cover values as numbers
xyplot(Cover~jitter(Elevation)|Form,data=sigcov)

sigcov_form<-aggregate(sigcov[,c(5)],sigcov[,c(2,3,4,13)],FUN=sum)#Total cover by growth form
names(sigcov_form)[5]<-'Cover'
bwplot(Cover~Elevation|Form,data=sigcov_form,horizontal = F)
plot(Cover~Elevation,data=sigcov_form)

sigcov_summ<-aggregate(sigcov[,c(1)],sigcov[,c(1,2,3,13)],FUN=length)#Number of species for each quadrat

form_div<-aggregate(forms[,c(3)],forms[,c(1,2,3)],length) #number of growth forms for each quadrat. 
plot(jitter(x)~jitter(Elevation),data=form_div,pch=16,ylab='# of Growth Forms/Quadrat',ylim=c(1,6),xlab='Elevation (m)',cex.lab=2,cex.axis=2)
boxplot(x~Elevation,data=form_div,pch=16,ylab='# of Growth Forms/Quadrat',ylim=c(1,6),xlab='Elevation (m)',cex.lab=2,cex.axis=1.25)


?par(mgp=c(5,1,0),mar=c(8,8,1,1),las=1)
des<-read.csv('Desc.csv') ##Site descriptive data (Inclination and Aspect)
summary(des)
summary(lm(Inclination~as.factor(Elevation),data=des))
summary(lm(Aspect~as.factor(Elevation),data=des))
plot(Aspect~jitter(Elevation,0.5),data=des[des$Transect=='A',],las=1,ylab='Aspect (° from N)',xlab='Elevation (m.a.s.l.)',pch=16,cex.lab=2.3,cex.axis=2,xlim=c(300,850))
points(Aspect~jitter(Elevation,0.5),data=des[des$Transect=='B',],col=2,pch=16)
text(520,140,'Black-Gradient A, Red-Gradient B',cex=1.7,font=2)
plot(Inclination~jitter(Elevation,0.5),data=des[des$Transect=='A',],ylab='Slope (° from horizontal)',xlab='Elevation (m.a.s.l.)',pch=16,cex.lab=2.3,cex.axis=2,xlim=c(300,850),ylim=c(0,60),yaxs='i')
points(Inclination~jitter(Elevation,0.5),data=des[des$Transect=='B',],col=2,pch=16)
text(520,5,'Black-Gradient A, Red-Gradient B',cex=1.7,font=2)

he<-read.csv('Heights.csv') ##Thallus heights
summary(he)
boxplot(Height~Elevation,data=he,notch=T)
xyplot(Height~Elevation|Species,data=he,horiz=F)
summary(lm(Height~Elevation*Species,data=he))
summary(lm(Height~Elevation+I(Elevation^2),data=he[he$Species=='Heterodermia1'&he$Elevation<800,]))
plot(Height~Elevation,data=he[he$Species=='Heterodermia1',],notch=T)

plot(Height~jitter(Elevation,0.5),data=he[he$Species=='Follm',],xlim=c(250,850),pch=16,cex=2,ylim=c(0,4.25),yaxs='i',xlab='Elevation (m)',ylab='Height (cm)')
summary(lm(Height~Elevation,data=he[he$Species=='Follm',]))
abline(a=-0.0688149,b=0.0009978,lty=2)
points(260,4,pch=16,cex=2)
text(400,4,'F. orthoclada: p « 0.001, r2=0.43')
points(Height~jitter(Elevation,0.5),data=he[he$Species=='RocCer',],xlim=c(250,850),pch=17,col=2,cex=2)
summary(lm(Height~Elevation,data=he[he$Species=='RocCer',]))
abline(a=0.673143,b=0.001677,col=2,lty=2)
points(260,3.7,pch=17,cex=2,col=2)
text(400,3.7,'R. cerebriformis: p < 0.005, r2=0.17')
points(Height~jitter(Elevation,0.5),data=he[he$Species=='Pentagenella',],xlim=c(250,850),pch=18,col=3,cex=2)
summary(lm(Height~Elevation,data=he[he$Species=='Pentagenella',]))
summary(lm(Height~Elevation+I(Elevation*Elevation),data=he[he$Species=='Pentagenella',]))
lines(polynomial(coef=c(-1.510e+01,5.391e-02,-4.109e-05)),lty=2,col=3)
#abline(a=0.673143,b=0.001677)
points(260,3.4,pch=18,cex=2,col=3)
text(400,3.4,'P. fragillima: p < 0.001, r2=0.26')

##version for grant proposal:
par(mfrow=c(1,1))
plot(Height~jitter(Elevation,0.5),data=he[he$Species=='Follm',],xlim=c(250,850),pch=16,cex=2,ylim=c(0,3.45),yaxs='i',xlab='Elevation (m)',ylab='Height (cm)',cex.axis=2,cex.lab=2.5)
abline(a=-0.0688149,b=0.0009978,lty=2)
points(260,3.25,pch=16,cex=2)
text(430,3.25,'Follmannia orthoclada',cex=1.75)
text(750,0.25,'p « 0.001, r2=0.43',cex=1.5)
points(Height~jitter(Elevation,0.5),data=he[he$Species=='RocCer',],xlim=c(250,850),pch=17,col=2,cex=2)
summary(lm(Height~Elevation,data=he[he$Species=='RocCer',]))
abline(a=0.673143,b=0.001677,col=2,lty=2)
points(260,3,pch=17,cex=2,col=2)
text(430,3,'Roccellina cerebriformis',cex=1.75)
text(650,2,'p < 0.005, r2=0.17',cex=1.5)


points(Height~jitter(Elevation,0.5),data=he[he$Species=='Niebla',],xlim=c(250,850),pch=1,col=4,cex=2)
summary(lm(Height~Elevation,data=he[he$Species=='Niebla',]))
points(Height~jitter(Elevation,0.5),data=he[he$Species=='HetFoll',],xlim=c(250,850),pch=2,col=5,cex=2)
summary(lm(Height~Elevation,data=he[he$Species=='HetFoll',]))

xyplot(Height~Elevation|Species,data=he,horiz=F)
#Note that for a number of species, there are only a few elevations sampled: this is probably too few to be useful. Also, the Ramalinas were not included in other trait studies. So we can reduce this down to a subset:
he_red<-subset(he,Species != c('Ramalina','Ramalina2','Rhizoplaca','Roccellina2','HetFoll','HetLeu'))
he_red<-subset(he,Species != 'Rhizoplaca'&Species != 'HetLeu'&Species != 'HetFoll'&Species != 'Ramalina'&Species != 'Ramalina2'&Species != 'Rhizoplaca'&Species != 'Roccellina2')
xyplot(Height~Elevation|Species,data=he_red)
summary(lm(Height~Elevation*Species,data=he_red))



rock<-read.csv('Rockcover.csv') ##Rock cover data, source of the summary stats below
rock_st<-read.csv('Rock_Cover_stats.csv') ##Rock cover descriptive statistics. 
summary(rock_st)

quad<-merge(des,rock_st[,c(2,3,5,6,7,8)],all=T)
quad<-merge(quad,cov_summ,all=T)
quad<-merge(quad,rouVar)
quad<-merge(quad,he,all=T)
summary(quad)
#write.csv(quad,'Quadrat_Master.csv')
#################################
###########################Part 2: explore the explanatory power of the quadrat attributes
#quad<-read.csv('Quadrat_Master.csv')

par(mfrow=c(2,2))
plot(R_Sum~Elevation,data=quad)
plot(R_Median~Elevation,data=quad)
plot(R_Max~Elevation,data=quad)
plot(Rou_SD~Elevation,data = quad)

summary(lm(NSpp~Elevation*R_Median*Rou_SD,data=quad))



####################Combine with microclimate data#####
ibutt<-read.csv('iButton_combined.csv')
########Need to figure out ways to handle different sampling frequencies across the period covered.
#The simplest is probably to aggregate by the hour values, although this has the effect of truncating rather than rounding in some cases:
ibutt<-aggregate(ibutt[,c(4,5,6,7)],by=list(Site=ibutt$Site,Hour=ibutt$Hour,Day=ibutt$Day,Month=ibutt$Month),mean)
ndays<-nrow(ibutt[ibutt$Hour==12&ibutt$Site==450,])
ndays750<-nrow(ibutt[ibutt$Hour==12&ibutt$Site==750,])
saturated99<-aggregate(ibutt[ibutt$RH>99,2],by=list(Site=ibutt[ibutt$RH>99,]$Site),length)
saturated97<-aggregate(ibutt[ibutt$RH>97,2],by=list(Site=ibutt[ibutt$RH>97,]$Site),length)
saturated0.01kPa<-aggregate(ibutt[ibutt$VPD<0.01,2],by=list(Site=ibutt[ibutt$VPD<0.01,]$Site),length)
saturated0.1kPa<-aggregate(ibutt[ibutt$VPD<0.1,2],by=list(Site=ibutt[ibutt$VPD<0.1,]$Site),length)
drystress1.0kPa<-aggregate(ibutt[ibutt$VPD>1.0,3],by=list(Site=ibutt[ibutt$VPD>1.0,]$Site),length)
drystress2.0kPa<-aggregate(ibutt[ibutt$VPD>2.0,3],by=list(Site=ibutt[ibutt$VPD>2.0,]$Site),length)
drystress3.0kPa<-aggregate(ibutt[ibutt$VPD>3.0,3],by=list(Site=ibutt[ibutt$VPD>3.0,]$Site),length)

names(saturated0.1kPa)<-c('Elevation','Period')
dry_spp<-merge(saturated0.1kPa,transect)
