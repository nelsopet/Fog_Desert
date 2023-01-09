
#Deprecated first part, retained as indicator of initial data handling
#chaiten<-read.csv('/home/daniel/Desktop/Isotopes/20190117_EA_Stanton_DES_ChaitenVolcanoLichen_processed.csv',skip=1)
# summary(chaiten)
# #some of the N and C values spike way above the calibration linearity range, these would suggest an instrumental problem, and are best trimmed:
# chaiten<-chaiten[chaiten$N.Peak.height<15,]
# #Then we can trim down to the columns and rows of interest
# chaiten<-chaiten[chaiten$X.2=='unknown',c(2,4:8)]
# names(chaiten)<-c('Sample','Mass','.N','d15N','.C','d13C')
# 
# 
# #Parse the sample names in site, replicate and taxon
# chaiten$Site<-NA
# chaiten$Taxon<-NA
# chaiten$Rep<-NA
# for (i in seq(1,nrow(chaiten),1)){ 
#  chaiten[i,7]<-unlist(strsplit(as.character(chaiten[i,1]),"\\.|\\_"))[1] #for example here we remove the first element from the sample name (the site)
#  chaiten[i,8]<-unlist(strsplit(as.character(chaiten[i,1]),"\\.|\\_"))[3]
#  chaiten[i,9]<-unlist(strsplit(as.character(chaiten[i,1]),"\\.|\\_"))[2]
# }
# chaiten$Site<-factor(chaiten$Site)
# chaiten$Taxon<-factor(chaiten$Taxon)
# chaiten$Rep<-factor(chaiten$Rep)

#write.csv(chaiten,'ChaitenIsotopes.csv')

#Manually merged with a second round of isotopes, that had slightly different file naming structure. Easier to do manually.

chaiten<-read.csv('/home/daniel/Dropbox/Chaiten/ChaitenIsotopes.csv',stringsAsFactors = T)

attrib<-read.csv('/home/daniel/Dropbox/Chaiten/ChaitenTaxa.csv',stringsAsFactors = T)
chaiten<-merge(chaiten,attrib,all=T)
summary(chaiten)

#So, focusing on those few taxa with decent replication, it looks like there is some variation between sites in 15N, at least for the fixer
plot(d15N~Site,data=chaiten[chaiten$Taxon=='k',],pch=16,cex=2,ylim=c(-10,0))
points(d15N~Site,data=chaiten[chaiten$Taxon=='c',],pch=16,cex=2,col='lightgreen')
points(d15N~Site,data=chaiten[chaiten$Taxon=='b',],pch=16,cex=2,col='darkgreen')
points(d15N~Site,data=chaiten[chaiten$Taxon=='f',],pch=16,cex=2,col='pink')

#A 2-dimensional view can sometimes help, so if we focus in on Flood1, we can see a trend towards higher 15N (ie greater contribution of recently fixed N) in plants with more N, suggesting that some are benefiting from fixation to boost their N content. That one outlier is a P. berberina, but with 1 rep we can't say much.
plot(d15N~.N,data=chaiten[chaiten$Site=='Flood1',])
points(d15N~.N,data=chaiten[chaiten$Site=='Flood1'&chaiten$Fixer=='yes',],pch=16,col='blue')
points(d15N~.N,data=chaiten[chaiten$Site=='Flood1'&chaiten$Type=='fern',],pch=16,col='lightgreen')
points(d15N~.N,data=chaiten[chaiten$Site=='Flood1'&chaiten$Type=='angio',],pch=16,col='pink')

plot(d15N~.N,data=chaiten[,])
points(d15N~.N,data=chaiten[chaiten$Fixer=='yes',],pch=16,col='blue')
points(d15N~.N,data=chaiten[chaiten$Type=='fern',],pch=16,col='lightgreen')
points(d15N~.N,data=chaiten[chaiten$Type=='angio',],pch=16,col='pink')
points(d15N~.N,data=chaiten[chaiten$Type=='lichen',],pch=16,col='yellow')

plot(d15N~d13C,data=chaiten[,])
points(d15N~d13C,data=chaiten[chaiten$Fixer=='yes',],pch=16,col='blue')
points(d15N~d13C,data=chaiten[chaiten$Type=='fern',],pch=16,col='lightgreen')
points(d15N~d13C,data=chaiten[chaiten$Type=='angio',],pch=16,col='pink')

plot(d13C~Site,data=chaiten[chaiten$Taxon=='k',],pch=16,cex=2,ylim=c(-40,-20))
points(d13C~Site,data=chaiten[chaiten$Taxon=='c',],pch=16,cex=2,col='lightgreen')
points(d13C~Site,data=chaiten[chaiten$Taxon=='b',],pch=16,cex=2,col='darkgreen')
points(d13C~Site,data=chaiten[chaiten$Taxon=='f',],pch=16,cex=2,col='pink')

library(ggplot2)
library(hrbrthemes)
ggplot(chaiten, aes(x=Site, y=d15N, color=Taxon)) +   geom_point(size=6)
