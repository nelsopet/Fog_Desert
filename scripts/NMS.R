
require(tidyverse)
require(vegan)


## Build a regular NMDS using only transect B, which has the complete community and environmental data
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
  coverB_flat_forMDS<-coverB_flat_forMDS %>% subset(tot_cov>0) %>% select(-tot_cov)
  coverB_MDS<-metaMDS(coverB_flat_forMDS, distance = "bray", try=100, trymax = 500)
  plot(coverB_MDS)  

## The NMDS above includes all quadrats, including those which had no lichens.
## Those empty plots need to be removed by calculating and filtering by total cover > 0
  coverB_filt<-subset(cover, Transect=='B') ## 'data.frame':	342 obs. of  8 variables:, 7 factors and 1 integer
  coverB_filt$Cover<-as.character(coverB_filt$Cover) ## same dimensions, one var turned to char
  ##Replace "T" meaning trace cover values iwth a number
  coverB_filt$Cover[coverB_filt$Cover=="T"]<-0.25
  ##Remove superfluous columns
  coverB_filt <- coverB_filt %>% select(-Collection, -Long_Form, -Photos) #'data.frame':	342 obs. of  5 variables:, 1 integer, 3 factors and 1 char
    ##Spread data so columns are species, values are cover and rows are single quadrats
    coverB_filt_flat<-coverB_filt %>% spread(Species, Cover) %>% select(-V1) # 'data.frame':	60 obs. of  92 variables: 1 integer, two factors, 89 chars
    ## Replace NAs with 0's and turn it into a dataframe
    coverB_filt_flat[is.na(coverB_filt_flat)==T]<-0
    ## Make the species cover columns numeric but this makes everything numeric so no need for the row subscript
    coverB_filt_flat$Quadrat<-as.character(coverB_filt_flat$Quadrat)
    coverB_filt_flat<-mapply(as.numeric, coverB_filt_flat) %>% as.data.frame()
    coverB_filt_flat[is.na(coverB_filt_flat)==T]<-0 # 'data.frame':	60 obs. of  92 variables:, all vars numeric
    ## Create rownames that are UIDs
    coverB_filt_flat$UID<-paste(coverB_filt_flat$Elevation, coverB_filt_flat$Quadrat, sep="_")
    rownames(coverB_filt_flat)<-coverB_filt_flat$UID
    coverB_filt_flat<-coverB_filt_flat %>% select(-UID)
      ## Remove columns that aren't species to conform to input requirements for metaMDS
      coverB_filt_flat_forMDS<-select(coverB_filt_flat, -Elevation, -Transect, -Quadrat)
      coverB_filt_flat_forMDS[is.na(coverB_filt_flat_forMDS)==T]<-0
      
      coverB_filt_flat_forMDS$rowsums<-coverB_filt_flat_forMDS[,4:nrow(coverB_filt_flat_forMDS)] %>% rowSums()
      coverB_filt_flat_forMDS<-cbind(coverB_filt_flat_forMDS,as.data.frame(rowSums(coverB_filt_flat_forMDS)))
      coverB_filt_flat_forMDS[is.na(coverB_filt_flat_forMDS)==T]<-0
      
      coverB_filt_flat_forMDS<-as.data.frame(coverB_filt_flat_forMDS) 
      coverB_filt_flat_forMDS<-rename(coverB_filt_flat_forMDS, tot_cov =`rowSums(coverB_filt_flat_forMDS)`)
      coverB_filt_flat_forMDS<-coverB_filt_flat_forMDS %>% subset(rowsums>0) %>% select(-tot_cov)#subset(tot_cov>0)# %>% select(-tot_cov)

coverB_MDS<-metaMDS(coverB_filt_flat_forMDS, distance = "bray", try=100, trymax = 500)

Env_matrix_B$UID<-paste(Env_matrix_B$Elevation, Env_matrix_B$Quadrat, sep="_")
#Need to filter Env_matrix to only have the 49 rows but now has rownames
stuff<-coverB_MDS$points %>% as.data.frame() %>% mutate(UID= rownames(coverB_MDS$points)) 

stuff_env<-Env_matrix_B %>% inner_join(stuff, by="UID")  %>% as.data.frame() # %>% select(VPD_min)

## Now this works and returns the r.sq
  stuff_out1<-ordisurf(coverB_MDS~Inclination, stuff_env, main= paste(colnames(stuff_env[5]))) %>% summary() %>% select(r.sq)
  stuff_out1$r.sq
  dev.off()
  
  
##automate ordisurf for all columns
## Need to unlist the output of each ordisurf ... str(stuff_out) resultsing a list of 54
## Error indicates input dimensions don't match
  vars<-colnames(stuff_env[5:24])  #%>% as.data.frame() 
  paste(vars[1])
  paste("stuff$",vars[2],sep="")   
  
  tst_hill<-ordisurf(eval(parse(text=paste("coverB_MDS~",vars[1],sep=""))), stuff_env)
  tst_hill_stats<-summary(tst_hill) 
  tst_hill_stats$r.sq
  dev.off()

  tst_hill<-ordisurf(eval(parse(text=paste("coverB_MDS~",vars[4],sep=""))), stuff_env)
  tst_hill_stats<-summary(tst_hill) 
  tst_hill_stats$r.sq
  dev.off()
  
  
    
  make_mods<-function (x) {paste("coverB_MDS~",vars[x],sep="") %>% as.formula()}
  var_mods<-lapply(1:length(vars),make_mods)
  
  eval(var_mods[1])
  
    ordi_fit<-function(x) {ordisurf(eval(parse(text=paste("coverB_MDS~",vars[x],sep=""))), stuff_env, main= paste(vars[x]))}
    
    length(stuff_env)
    
    ## Why does this produce a vector of numbers read as characters?
    #vars<-as.character(vars)
    pdf("Patache_Hilltop_Plots.pdf")
    lapply(1:length(vars),ordi_fit)
    dev.off()
    
    ordi_stats<-function(x)
    {out1<-ordisurf(eval(parse(text=paste("coverB_MDS~",vars[x],sep=""))), stuff_env) %>% summary()
    return(round(out1$r.sq,2))}
    
    ordi_stats_out<-lapply(1:length(vars),ordi_stats)
    ordi_stats_out<-cbind(vars,ordi_stats_out)
    colnames(ordi_stats_out)<-c("Env Variable","R squared")
    write.csv(ordi_stats_out, "Patached_TransectB_CommunityEnv_NMS_GAM_fit.csv")
######### Second matrix ... Daniel recommeded starting looking at microclimate, specifically
######### sat85_dry, Tmed, VPDmed ... where are these columns?
## Find files with data from same quadrats/elevations on transect B
## anat, he, he_red, quad, rou, rouVar all are filterable by transect
## angles, dry only has data for transect B
## unclear which transect rock and rock_st came from
