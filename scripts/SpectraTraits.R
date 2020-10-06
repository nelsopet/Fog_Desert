setwd("~/Dropbox/Lomas/Nat Geo/Combined Data")
patache.spec<-read.csv('PatacheSpectra.csv')

#there are waay too many columns to be workable (~2000!), so let's trim down to the ones that might help us calculate some indices
pata.spec<-patache.spec[,c('Field.Name','Replicate','Elevation','X300.373','X532.283','X570.723','X630.393','X670.715','X800.683','X930.070')]
pata.spec$NDVI<-(pata.spec$X800.683-pata.spec$X630.393)/(pata.spec$X800.683+pata.spec$X630.393) #NDVI
pata.spec$PRI<-(pata.spec$X532.283-pata.spec$X570.723)/(pata.spec$X532.283+pata.spec$X570.723) #PRI
pata.spec$UVI<-(pata.spec$X800.683-pata.spec$X300.373)/(pata.spec$X800.683+pata.spec$X300.373) #guess at UV index
pata.spec$wet<-(pata.spec$X930.070-pata.spec$X670.715)/(pata.spec$X930.070+pata.spec$X670.715) #guess at wet response

summary(pata.spec) #some adjustments need to be made to the names manually to align them with the rest of the Patache data:
write.csv(pata.spec,'Patache_Spec_Sum1.csv')

pat.spc<-read.csv('Patache_Spec_Sum.csv')
taxa<-read.csv('Taxa_TransectB.csv')

specs<-merge(pat.spc,taxa,all=T) #these then need to be cleaned for duplicate records
specs<-read.csv('Patache_Spec_Clean.csv',stringsAsFactors = T)


############Alternative approach (retaining more data): see PatacheSpectra.R


