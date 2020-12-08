#install.packages('devtools')
#Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")

devtools::install_github('phytomosaic/foggy')
devtools::install_github('phytomosaic/ecole')

require(tidyverse)
require(vegan)
require(Metrics)
require(hsdar)
## Read in cover data
cover<-read.csv('data/CoverMaster.csv')
  cover$UID<-paste(cover$Elevation, cover$Quadrat, sep="_")
    dim(cover) #[1] 1522    9

# read in species-level traits
traits_sp<-read.csv('data/Traits_TransB.csv')
  View(traits_sp)
    dim(traits_sp) #[1] 356  26
      traits_sp %>% 
        dplyr::select(Species, Photobiont) %>%
          unique() %>%
            write.csv("traits_species_photobiont_for_review.csv")
colnames(traits_sp)


##where is the Env data read in?
Env_matrix_B<-read.csv('data/Env_matrix_B.csv')

## Read in trait data
#trait<-read.csv('data/QuadratMeanTraits_Weighted_NoTrace.csv')
trait<-read.csv('data/QuadratMeanTraits_Unweighted.csv')
trait$UID<-paste(trait$Elevation, trait$Quadrat, sep="_")
View(trait)
trait<-trait[,-1:-4]
anti_join(cover, trait, by="UID") %>% nrow()
nrow(Env_matrix_B)
# 69 species with no trait data
anti_join(cover, trait, by="Species") %>% select(Species) %>% unique() %>% nrow()

#Read in spectra by quadrat
spectra<-read.csv('data/Spectral_Features.csv')
View(spectra)
dim(spectra) # [1] 130 115


## Build a regular NMDS using only transect B, which has the complete community and environmental data
  coverB<-subset(cover, Transect=='B')
    coverB$Cover<-as.character(coverB$Cover)
      coverB$Cover[coverB$Cover=="T"]<-0.25
      
      #Check to see if all the taxa in coverB, the subset of the cover data used, have traits to filter by
      biont = c("Trebouxia","Trentepohlia")
      traits_sp_filt<-traits_sp %>% dplyr::select(Species, Photobiont) %>% unique() %>% as.data.frame()
      coverB_Treb<-coverB %>% 
        left_join(traits_sp_filt,by="Species", keep=F) %>% 
        subset(Photobiont %in% biont) %>%
        subset(Photobiont == "Trebouxia") #%>% View()
      #()#only 11 rows don't match #[1] 11
      coverB_Trent<-coverB %>% 
        left_join(traits_sp_filt,by="Species", keep=F) %>% #dim()
        subset(Photobiont %in% biont) %>%
        subset(Photobiont == "Trentepohlia")
      
      
  coverB_flat<-coverB %>% dplyr::select(-Collection, -Long_Form, -Photos) %>% spread(Species, Cover) #%>% dim()
  coverB_flat_forMDS<-dplyr::select(coverB_flat, -Elevation, -Transect, -Quadrat,-V1) #%>% colnames()
    coverB_flat_forMDS<-sapply(coverB_flat_forMDS, as.numeric) %>% as.data.frame()
      coverB_flat_forMDS[is.na(coverB_flat_forMDS)]<-0
        coverB_flat_forMDS<-cbind(coverB_flat_forMDS,as.data.frame(rowSums(coverB_flat_forMDS)))
          coverB_flat_forMDS<-as.data.frame(coverB_flat_forMDS) 
            coverB_flat_forMDS<-rename(coverB_flat_forMDS, tot_cov =`rowSums(coverB_flat_forMDS)`)
              coverB_flat_forMDS<-coverB_flat_forMDS %>% subset(tot_cov>0) %>% dplyr::select(-tot_cov)
  
  #Trebouxia matrix
              #LEFT OFF HERE
   coverB_Treb_flat<- coverB_Treb %>% #dplyr::select(-Collection, -Long_Form, -Photos, -Photobiont) %>% spread(Species, Cover)
     #rename(Elevation = Elevation.x, Transect = Transect.x, Quadrat = Quadrat.x, Cover = Cover.x) %>% #View()
      dplyr::select(Elevation, Transect, Quadrat, UID, Species, Cover) %>% #group_by(Cover) %>% tally()
        pivot_wider(-Elevation, names_from = Species, values_from = Cover) #%>% View() 
   #coverB_Treb_flat_forMDS<-dplyr::select(coverB_Treb_flat, -Transect, -Quadrat, -Elevation) %>% colnames()
   coverB_Treb_flat_forMDS[is.na(coverB_Treb_flat_forMDS)]<-0
   coverB_Treb_flat_forMDS<-sapply(coverB_Treb_flat_forMDS, as.numeric) %>% as.data.frame()
   coverB_Treb_flat_forMDS<-cbind(coverB_Treb_flat_forMDS,as.data.frame(rowSums(coverB_Treb_flat_forMDS)))
   coverB_Treb_flat_forMDS<-as.data.frame(coverB_Treb_flat_forMDS) 
   coverB_Treb_flat_forMDS<-rename(coverB_Treb_flat_forMDS, tot_cov =`rowSums(coverB_Treb_flat_forMDS)`)
   coverB_Treb_flat_forMDS<-coverB_Treb_flat_forMDS %>% subset(tot_cov>0) %>% dplyr::select(-tot_cov)
        coverB_Treb_flat_forMDS %>% group_by(tot_cov) %>% tally()
    
    #Trentepohlia matrix
   coverB_Trent_flat<-coverB_Trent %>% #colnames()
     #rename(Elevation = Elevation.x, Transect = Transect.x, Quadrat = Quadrat.x, Cover = Cover.x) %>% #View()
     dplyr::select(Elevation, Transect, Quadrat, UID, Species, Cover) %>% #head()
     pivot_wider( names_from = Species, values_from = Cover) #%>% View()
   coverB_Trent_flat_forMDS<-dplyr::select(coverB_Trent_flat, -Elevation, -Transect, -Quadrat)
   coverB_Trent_flat_forMDS<-sapply(coverB_Trent_flat_forMDS, as.numeric) %>% as.data.frame()
   coverB_Trent_flat_forMDS[is.na(coverB_Trent_flat_forMDS)]<-0
   coverB_Trent_flat_forMDS<-cbind(coverB_Trent_flat_forMDS,as.data.frame(rowSums(coverB_Trent_flat_forMDS)))
   coverB_Trent_flat_forMDS<-as.data.frame(coverB_Trent_flat_forMDS) 
   coverB_Trent_flat_forMDS<-rename(coverB_Trent_flat_forMDS, tot_cov =`rowSums(coverB_Trent_flat_forMDS)`)
   coverB_Trent_flat_forMDS<-coverB_Trent_flat_forMDS %>% subset(tot_cov>0) %>% dplyr::select(-tot_cov)
   
              
              
  coverB_MDS<-metaMDS(coverB_flat_forMDS, distance = "bray", try=500, trymax = 1000)
    plot(coverB_MDS)  
      dev.off()
  coverB_Treb_MDS<-metaMDS(coverB_Treb_flat_forMDS, distance = "bray", try=500, trymax = 1000)
    plot(coverB_Treb_MDS)  
       dev.off()
  coverB_Trent_MDS<-metaMDS(coverB_Trent_flat_forMDS, distance = "bray", try=500, trymax = 1000)
    plot(coverB_Trent_MDS)  
      dev.off()
   #### LEFT OFF HERE     
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

## Run NMDS with relative distance (Bray) measure
coverB_MDS<-metaMDS(coverB_filt_flat_forMDS, distance = "bray", try=100, trymax = 500)
plot(coverB_MDS)
# Make the a column in the NMDS output that is the quadrat UID to be able to join with the env data 
#Need to filter Env_matrix to only have the 49 rows but now has rownames
stuff<-coverB_MDS$points %>% as.data.frame() %>% mutate(UID= rownames(coverB_MDS$points)) 

## Make the Env data have a column that makes a unique identifier 
Env_matrix_B$UID<-paste(Env_matrix_B$Elevation, Env_matrix_B$Quadrat, sep="_")

## Join the NMDS data with the 
stuff_env<-Env_matrix_B %>% inner_join(stuff, by="UID")  %>% as.data.frame() # %>% select(VPD_min)
  #Allrows are in both dasetss
  #stuff_env %>% anti_join(trait, by = "UID") %>% nrow()
    #[1] 0
## Add the traits to the env matrix
stuff_env<-stuff_env %>% 
  inner_join(trait, by="UID")  %>% 
  as.data.frame() %>% 
  select(UID, X.1, Transect, Quadrat, Elevation, everything())
colnames(stuff_env)
  #Write out summary of environmental data
  stuff_env %>% View()
## Make a list of env variable names to use in functions later
vars<-colnames(stuff_env[4:length(stuff_env)])  #%>% as.data.frame() 
colnames(stuff_env)
##PASS: Unit test of using the GAM model of Env var ~ NMS Axis 1 + Axis 2 and returns the r.sq
  stuff_out1<-ordisurf(coverB_MDS~Inclination, stuff_env, main= paste(colnames(stuff_env[5]))) %>% summary() #%>% select(r.sq)
  #stuff_out1$r.sq
  #dev.off()
  
  
## automate ordisurf for all columns
## Need to unlist the output of each ordisurf ... str(stuff_out) resultsing a list of 54
  ##Unit test for hilstop plot and  grabbing R2 PASSES
  tst_hill<-ordisurf(eval(parse(text=paste("coverB_MDS~",vars[2],sep=""))), stuff_env)
  tst_hill_stats<-summary(tst_hill) 
  tst_hill_stats$r.sq
  #dev.off()

  #tst_hill<-ordisurf(eval(parse(text=paste("coverB_MDS~",vars[4],sep=""))), stuff_env, )
  #tst_hill_stats<-summary(tst_hill) 
  #tst_hill_stats$r.sq
  #str(tst_hill_stats)
  #dev.off()
  
  
  ##Function to make models ... not sure if it even still get used though
  make_mods<-function (x) {paste("coverB_MDS~",vars[x],sep="") %>% as.formula()}
  var_mods<-lapply(1:length(vars),make_mods)
    eval(var_mods[1])
    ##Unit test for making linear biplot vector for elevation PASSSES
    tst1<-envfit(eval(parse(text=paste("coverB_MDS~",vars[2],sep=""))), stuff_env, main= paste(vars[2]))
    plot(tst1)
    dev.off()
    #tst2<-elev_hilltop<-ordisurf(eval(parse(text=paste("coverB_MDS~",vars[2],sep=""))), stuff_env, main= paste(vars[2]),labcex=0, col='black')
    #elev_hilltop
    #ordisurf(eval(parse(text=paste("coverB_MDS~",vars[4],sep=""))), stuff_env, main= paste(vars[4]),labcex=0, add=T)
    #text(max(stuff$MDS1)*0.9, max(stuff$MDS2)*0.95, paste("R2",ordi_stats_out[2,2])) 
    #dev.off()
    #tst2
##Make R2, give them sensible names and write out as CSV
      ordi_stats<-function(x)
      {
      out1<-ordisurf(eval(parse(text=paste("coverB_MDS~",vars[x],sep=""))), stuff_env) %>% summary()
      return(round(out1$r.sq,2))
      return(round(out1$se,2))
      }
        ordi_stats_out<-lapply(1:length(vars),ordi_stats) %>% unlist()
          ordi_stats_out<-cbind(vars,ordi_stats_out) %>% as.data.frame()
            colnames(ordi_stats_out)<-c("Env Variable","R squared")
              ##Now write.csv doesn't work!?
              write_csv(ordi_stats_out, "Patached_TransectB_CommunityEnv_NMS_GAM_fit.csv")
    
#Make hilltop plots with R2 
    ordi_fit<-function(x) 
            {
      ##jpeg(paste("Patache_Hilltop",vars[x],".jpeg"))
      ordisurf(eval(parse(text=paste("coverB_MDS~",vars[2],sep=""))), stuff_env, main= "",labcex=0, col='black')  
        ordisurf(eval(parse(text=paste("coverB_MDS~",vars[x],sep=""))), stuff_env, main= paste(vars[x]),labcex=1, add=T)
          text(max(stuff$MDS1)*0.5, max(stuff$MDS2)*0.95, paste("R2=",ordi_stats_out[x,2]), cex=2)
           #elev<-envfit(eval(parse(text=paste("coverB_MDS~",vars[2],sep=""))), stuff_env, main= paste(vars[2]))
           # plot(elev, cex=1.5)
      #dev.off()
          
            }
    
    pdf("Patache_Hilltop_Plots.pdf")
    lapply(1:length(vars),ordi_fit)
    dev.off()
    
######### Second matrix ... Daniel recommeded starting looking at microclimate, specifically
######### sat85_dry, Tmed, VPDmed ... where are these columns?
## Find files with data from same quadrats/elevations on transect B
## anat, he, he_red, quad, rou, rouVar all are filterable by transect
## angles, dry only has data for transect B
## unclear which transect rock and rock_st came from
    stuff_env[,-1:-4] %>% cor(use="pairwise.complete.obs") %>% cov2cor() %>% view()