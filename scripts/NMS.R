coverB<-subset(cover, Transect=='B')
coverB$Cover<-as.character(coverB$Cover)
coverB$Cover[coverB$Cover=="T"]<-0.25

coverB_flat<-coverB %>% select(-Collection, -Long_Form, -Photos) %>% spread(Species, Cover)
coverB_flat_forMDS<-select(coverB_flat, -Elevation, -Transect, -Quadrat,-V1)
coverB_flat_forMDS<-sapply(coverB_flat_forMDS, as.numeric) %>% as.data.frame()
coverB_flat_forMDS[is.na(coverB_flat_forMDS)]<-0

coverB_flat_forMDS<-cbind(coverB_flat_forMDS,as.data.frame(rowSums(coverB_flat_forMDS)))
coverB_flat_forMDS<-as.data.frame(coverB_flat_forMDS) 

coverB_flat_forMDS<-rename(coverB_flat_forMDS, tot_cov =`rowSums(coverB_flat_forMDS)`)
coverB_flat_forMDS<-coverB_flat_forMDS %>% subset(tot_cov>0)

plot(metaMDS(coverB_flat_forMDS, distance = "bray"))
