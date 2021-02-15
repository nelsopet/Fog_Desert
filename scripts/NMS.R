##### Required packages
#install.packages('devtools')
#Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")

devtools::install_github('phytomosaic/foggy')
devtools::install_github('phytomosaic/ecole')

require(tidyverse)
require(vegan)
require(Metrics)
require(hsdar)
##### Read in data
## Read in cover data
cover<-read.csv('data/CoverMaster.csv')
  cover$UID<-paste(cover$Elevation, cover$Quadrat, sep="_")
    dim(cover) #[1] 1522    9
      colnames(cover)
      
# read in species-level traits, including absorption time and angles
traits_sp<-read.csv('data/Traits_TransB.csv')
  View(traits_sp)
    dim(traits_sp) #[1] 356  26
      traits_sp %>% 
        dplyr::select(Species, Photobiont) %>%
          unique() %>%
            write.csv("traits_species_photobiont_for_review.csv")

##Read in the environment data, including climate derivatives
Env_matrix_B<-read.csv('data/Env_matrix_B.csv')
  ## Make the Env data have a column that makes a unique identifier 
  Env_matrix_B$UID<-paste(Env_matrix_B$Elevation, Env_matrix_B$Quadrat, sep="_")
    dim(Env_matrix_B) #[1] 60 30
    
# USE THIS TRAIT MATRIX: Read in trait data, including all the spectra, absorption time, number of trent, treb photobionts
#trait<-read.csv('data/QuadratMeanTraits_Weighted_NoTrace.csv')
trait<-read.csv('data/QuadratMeanTraits_Unweighted.csv')
  trait$UID<-paste(trait$Elevation, trait$Quadrat, sep="_")
    trait<-trait[,-1:-4]
      dim(trait) #[1]  49 123 
        colnames(trait)
      
#Read in spectra by quadrat, including band ratios
spectra<-read.csv('data/Spectral_Features.csv')
  dim(spectra) #[1] 130 115 
  
#Read in substrate data, including quadrat-level rockiness measurements
  #QUESTION: Joint this the traits?
substrate<-read.csv('data/Substrate_PatacheB.csv',header=T)
  substrate$UID<-paste(substrate$Elevation, substrate$Quadrat, sep="_")
    dim(substrate) #[1] 60 10 

##### NMS 
# Build a regular NMDS using only transect B, which has the complete community and environmental data
  coverB<-subset(cover, Transect=='B')
    coverB$Cover<-as.character(coverB$Cover)
      coverB$Cover[coverB$Cover=="T"]<-0.25
      
    #Check to see if all the taxa in coverB, the subset of the cover data used, have traits to filter by
    biont = c("Trebouxia","Trentepohlia")
    
    traits_sp_filt<-traits_sp %>% 
      dplyr::select(Species, Photobiont) %>% 
        unique() %>% 
          as.data.frame()
    
    coverB_Treb<-coverB %>% 
      left_join(traits_sp_filt,by="Species", keep=F) %>% 
        subset(Photobiont %in% biont) %>%
          subset(Photobiont == "Trebouxia") #%>% View()
      #()#only 11 rows don't match #[1] 11
    
    coverB_Trent<-coverB %>% 
      left_join(traits_sp_filt,by="Species", keep=F) %>% #dim()
        subset(Photobiont %in% biont) %>%
          subset(Photobiont == "Trentepohlia")

    coverB_flat_forMDS<-
    coverB %>% 
    dplyr::filter(Species!="") %>% #nrow()
    dplyr::select(-Collection, -Long_Form, -Photos) %>% 
    pivot_wider(names_from=Species, values_from=Cover) %>% #colnames()
      dplyr::select(-Elevation, -Transect, -Quadrat) %>%
      mutate(across(c(-UID),as.numeric)) %>%
      replace(is.na(.),0) %>%
      rowwise() %>% 
      mutate(tot_cov=sum(across(c(-UID)))) %>% 
      dplyr::filter(tot_cov>0) %>%
      dplyr::select(-tot_cov,-UID)
    
  
      #OLD coverB_flat<-
      #OLD coverB %>% 
      #OLD dplyr::filter(Species!="") %>% #nrow()
      #OLD dplyr::select(-Collection, -Long_Form, -Photos) %>% 
      #OLD spread(Species, Cover) #%>% str() #'data.frame':	60 obs. of  94 variables:
      #OLD coverB_flat_forMDS<-dplyr::select(coverB_flat, -Elevation, -Transect, -Quadrat,-V1)
      #OLD: coverB_flat_forMDS<-sapply(coverB_flat_forMDS, as.numeric) %>% as.data.frame()
      #OLD coverB_flat_forMDS[is.na(coverB_flat_forMDS)]<-0
      #OLD coverB_flat_forMDS<-cbind(coverB_flat_forMDS,as.data.frame(rowSums(coverB_flat_forMDS)))
      #OLDcoverB_flat_forMDS<-as.data.frame(coverB_flat_forMDS) 
      #OLD coverB_flat_forMDS<-rename(coverB_flat_forMDS, tot_cov =`rowSums(coverB_flat_forMDS)`)
      #OLD coverB_flat_forMDS<-coverB_flat_forMDS %>% subset(tot_cov>0) %>% dplyr::select(-tot_cov)
  
  #Trebouxia matrix
   coverB_Treb_flat_forMDS<- 
     coverB_Treb %>% #colnames()
      dplyr::select(Elevation, Transect, Quadrat, UID, Species, Cover) %>% #colnames()
      #OLD pivot_wider(-Elevation, names_from = Species, values_from = Cover) #UID still here
     dplyr::filter(Species!="") %>% #nrow()
     #dplyr::select(-Collection, -Long_Form, -Photos) %>% 
     pivot_wider(names_from=Species, values_from=Cover) %>% #colnames()
     dplyr::select(-Elevation, -Transect, -Quadrat) %>% #colnames()
     mutate(across(c(-UID),as.numeric)) %>% #colnames()
     replace(is.na(.),0) %>% #dim()
     rowwise() %>% 
     mutate(tot_cov=sum(across(c(-UID)))) %>% 
     dplyr::filter(tot_cov>0) %>%
     dplyr::select(-tot_cov,-UID)
   
   names_nmds_Treb<-coverB_Treb %>% #colnames()
     dplyr::select(Elevation, Transect, Quadrat, UID, Species, Cover) %>% #colnames()
     #OLD pivot_wider(-Elevation, names_from = Species, values_from = Cover) #UID still here
     dplyr::filter(Species!="") %>% #nrow()
     #dplyr::select(-Collection, -Long_Form, -Photos) %>% 
     pivot_wider(names_from=Species, values_from=Cover) %>% #colnames()
     dplyr::select(-Elevation, -Transect, -Quadrat) %>% #colnames()
     mutate(across(c(-UID),as.numeric)) %>% #colnames()
     replace(is.na(.),0) %>% #dim()
     rowwise() %>% 
     mutate(tot_cov=sum(across(c(-UID)))) %>% 
     dplyr::filter(tot_cov>0) %>%
     dplyr::select(UID)
            
   #OLD coverB_Treb_flat_forMDS<-sapply(coverB_Treb_flat_forMDS, as.numeric) %>% as.data.frame() ##
   #OLD  coverB_Treb_flat_forMDS[is.na(coverB_Treb_flat_forMDS)==T]<-0
   #OLD    coverB_Treb_flat_forMDS<-cbind(coverB_Treb_flat_forMDS,as.data.frame(rowSums(coverB_Treb_flat_forMDS)))
   #OLD      coverB_Treb_flat_forMDS<-as.data.frame(coverB_Treb_flat_forMDS) 
   #OLD        coverB_Treb_flat_forMDS<-rename(coverB_Treb_flat_forMDS, tot_cov =`rowSums(coverB_Treb_flat_forMDS)`)
   #OLD          coverB_Treb_flat_forMDS<-coverB_Treb_flat_forMDS %>% subset(tot_cov>0) %>% dplyr::select(-tot_cov)
            
  #Trentepohlia matrix
   coverB_Trent_flat_forMDS<- 
     coverB_Trent %>% #colnames()
     dplyr::select(Elevation, Transect, Quadrat, UID, Species, Cover) %>% #colnames()
     #OLD pivot_wider(-Elevation, names_from = Species, values_from = Cover) #UID still here
     dplyr::filter(Species!="") %>% #nrow()
     #dplyr::select(-Collection, -Long_Form, -Photos) %>% 
     pivot_wider(names_from=Species, values_from=Cover) %>% #colnames()
     dplyr::select(-Elevation, -Transect, -Quadrat) %>% #colnames()
     mutate(across(c(-UID),as.numeric)) %>% #colnames()
     replace(is.na(.),0) %>% #dim()
     rowwise() %>% 
     mutate(tot_cov=sum(across(c(-UID)))) %>% 
     dplyr::filter(tot_cov>0) %>%
     dplyr::select(-tot_cov,-UID)
   
   
   names_nmds_Trent<-coverB_Trent %>% #colnames()
     dplyr::select(Elevation, Transect, Quadrat, UID, Species, Cover) %>% #colnames()
     #OLD pivot_wider(-Elevation, names_from = Species, values_from = Cover) #UID still here
     dplyr::filter(Species!="") %>% #nrow()
     #dplyr::select(-Collection, -Long_Form, -Photos) %>% 
     pivot_wider(names_from=Species, values_from=Cover) %>% #colnames()
     dplyr::select(-Elevation, -Transect, -Quadrat) %>% #colnames()
     mutate(across(c(-UID),as.numeric)) %>% #colnames()
     replace(is.na(.),0) %>% #dim()
     rowwise() %>% 
     mutate(tot_cov=sum(across(c(-UID)))) %>% 
     dplyr::filter(tot_cov>0) %>%
     dplyr::select(UID)
   
   
   
   #OLD coverB_Trent_flat<-coverB_Trent %>% #colnames()
   #OLD #rename(Elevation = Elevation.x, Transect = Transect.x, Quadrat = Quadrat.x, Cover = Cover.x) %>% #View()
   #OLD dplyr::select(Elevation, Transect, Quadrat, UID, Species, Cover) %>% #head()
   #OLD  pivot_wider( names_from = Species, values_from = Cover) #%>% View()
   #OLD    coverB_Trent_flat_forMDS<-dplyr::select(coverB_Trent_flat, -Elevation, -Transect, -Quadrat) #%>% str() #tibble [43 × 26] (S3: tbl_df/tbl/data.frame)
   #OLD      colnames(coverB_Trent_flat_forMDS)
   #OLD      #Creates NAs unnecessarily for the UID column
   #OLD      tst<-coverB_Trent_flat_forMDS$UID %>% 
   #OLD      cbind(sapply(coverB_Trent_flat_forMDS[,2:25], as.numeric)) #%>% dim()
   #OLD      #coverB_Trent_flat_forMDS<-tst
   #OLD        
   #OLD      #coverB_Trent_flat_forMDS<-sapply(coverB_Trent_flat_forMDS, as.numeric) %>% as.data.frame()
   #OLD      coverB_Trent_flat_forMDS[is.na(coverB_Trent_flat_forMDS)]<-0
   #OLD        tst[is.na(tst)]<-0
   #OLD          tst<-as.data.frame(tst)
   #OLD            rowSums(tst)
   #OLD              str(tst)
   #OLD      rowSums(sapply(coverB_Trent_flat_forMDS[,2:25], as.numeric))
   #OLD        coverB_Trent_flat_forMDS<-cbind(coverB_Trent_flat_forMDS,as.data.frame(rowSums(coverB_Trent_flat_forMDS)))
   #OLD          coverB_Trent_flat_forMDS<-as.data.frame(coverB_Trent_flat_forMDS) 
   #OLD            coverB_Trent_flat_forMDS<-rename(coverB_Trent_flat_forMDS, tot_cov =`rowSums(coverB_Trent_flat_forMDS)`)
   #OLD              coverB_Trent_flat_forMDS<-coverB_Trent_flat_forMDS %>% subset(tot_cov>0) %>% dplyr::select(-tot_cov)
                    
                  
                  
  #Run NMDS with Bray-Curtis distance measure with more than necessary number of min/max tries to avoid locale optima
  coverB_MDS<-metaMDS(coverB_flat_forMDS, distance = "bray", try=500, trymax = 1000)
    plot(coverB_MDS)  
      dev.off()
      
      
      #Shephard plots and stress vs dim from Stackoverflow post
      #https://stackoverflow.com/questions/49223740/cumulative-variance-explained-for-nmds-in-r/49224573
      
      #par(mfrow = c(1,2), mar = c(3.5,3.5,3,1), mgp = c(2, 0.6, 0), cex = 0.8, las = 1)
      spear <- round(cor(vegdist(coverB_flat_forMDS, method = "bray"), dist(coverB_MDS$points), method = "spearman"),3)
      plot(vegdist(coverB_flat_forMDS, method = "bray"), dist(coverB_MDS$points), main = "Shepard diagram of coverB_MDS", 
           xlab = "True Bray-Curtis distance", ylab = "Distance in the reduced space")
      mtext(line = 0.1, text = paste0("Spearman correlation = ", spear), cex = 0.7)
      
      n = 10
      stress <- vector(length = n)
      for (i in 1:n) {
        stress[i] <- metaMDS(coverB_flat_forMDS, distance = "bray", k = i)$stress
      }
      names(stress) <- paste0(1:n, "Dim")
      # x11(width = 10/2.54, height = 7/2.54)
      par(mar = c(3.5,3.5,1,1), mgp = c(2, 0.6, 0), cex = 0.8, las = 2)
      barplot(stress, ylab = "stress")
      
  coverB_Treb_MDS<-metaMDS(coverB_Treb_flat_forMDS, distance = "bray", try=500, trymax = 1000)
    plot(coverB_Treb_MDS)  
       dev.off()
       
       #Shephard plots and stress vs dim from Stackoverflow post
       #https://stackoverflow.com/questions/49223740/cumulative-variance-explained-for-nmds-in-r/49224573
       
       spear_Treb <- round(cor(vegdist(coverB_Treb_flat_forMDS, method = "bray"), dist(coverB_Treb_MDS$points), method = "spearman"),3)
       plot(vegdist(coverB_Treb_flat_forMDS, method = "bray"), dist(coverB_Treb_MDS$points), main = "Shepard diagram of coverB_MDS", 
            xlab = "True Bray-Curtis distance", ylab = "Distance in the reduced space")
       mtext(line = 0.1, text = paste0("Spearman correlation = ", spear_Treb), cex = 0.7)
       
       
       n = 10
       stress <- vector(length = n)
       for (i in 1:n) {
         stress[i] <- metaMDS(coverB_Treb_flat_forMDS, distance = "bray", k = i)$stress
       }
       names(stress) <- paste0(1:n, "Dim")
       # x11(width = 10/2.54, height = 7/2.54)
       par(mar = c(3.5,3.5,1,1), mgp = c(2, 0.6, 0), cex = 0.8, las = 2)
       barplot(stress, ylab = "stress")
       
  coverB_Trent_MDS<-metaMDS(coverB_Trent_flat_forMDS, distance = "bray", try=500, trymax = 1000)
    plot(coverB_Trent_MDS)  
      dev.off()
      
      #Shephard plots and stress vs dim from Stackoverflow post
      #https://stackoverflow.com/questions/49223740/cumulative-variance-explained-for-nmds-in-r/49224573
      
      spear_Trent <- round(cor(vegdist(coverB_Trent_flat_forMDS, method = "bray"), dist(coverB_Trent_MDS$points), method = "spearman"),3)
      plot(vegdist(coverB_Trent_flat_forMDS, method = "bray"), dist(coverB_Trent_MDS$points), main = "Shepard diagram of coverB_MDS", 
           xlab = "True Bray-Curtis distance", ylab = "Distance in the reduced space")
      mtext(line = 0.1, text = paste0("Spearman correlation = ", spear_Trent), cex = 0.7)
      
      
      n = 10
      stress <- vector(length = n)
      for (i in 1:n) {
        stress[i] <- metaMDS(coverB_Trent_flat_forMDS, distance = "bray", k = i)$stress
      }
      names(stress) <- paste0(1:n, "Dim")
      # x11(width = 10/2.54, height = 7/2.54)
      par(mar = c(3.5,3.5,1,1), mgp = c(2, 0.6, 0), cex = 0.8, las = 2)
      barplot(stress, ylab = "stress")

#OLD ## The NMDS above includes all quadrats, including those which had no lichens.
#OLD       #QUESTION: Is this still true? I filtero out quadrats with cover<0
#OLD ## Those empty plots need to be removed by calculating and filtering by total cover > 0
#OLD   coverB_filt<-subset(cover, Transect=='B') ## 'data.frame':	342 obs. of  8 variables:, 7 factors and 1 integer
#OLD     coverB_filt$Cover<-as.character(coverB_filt$Cover) ## same dimensions, one var turned to char
#OLD       ##Replace "T" meaning trace cover values iwth a number
#OLD         coverB_filt$Cover[coverB_filt$Cover=="T"]<-0.25
#OLD           ##Remove superfluous columns
#OLD             coverB_filt <- coverB_filt %>% 
#OLD               dplyr::select(-Collection, -Long_Form, -Photos) #'data.frame':	342 obs. of  5 variables:, 1 integer, 3 factors and 1 char
#OLD     ##Spread data so columns are species, values are cover and rows are single quadrats
#OLD     coverB_filt_flat<-coverB_filt %>% 
#OLD       spread(Species, Cover) %>% 
#OLD         dplyr::select(-V1) # 'data.frame':	60 obs. of  92 variables: 1 integer, two factors, 89 chars
#OLD           ## Replace NAs with 0's and turn it into a dataframe
#OLD           coverB_filt_flat[is.na(coverB_filt_flat)==T]<-0
#OLD             ## Make the species cover columns numeric but this makes everything numeric so no need for the row subscript
#OLD             coverB_filt_flat$Quadrat<-as.character(coverB_filt_flat$Quadrat)
#OLD               names<-coverB_filt_flat$UID
#OLD             
#OLD                coverB_filt_flat<-lapply(coverB_filt_flat[-1:-4],as.numeric) %>% as.data.frame()
#OLD               
#OLD                 coverB_filt_flat[is.na(coverB_filt_flat)==T]<-0 # 'data.frame':	60 obs. of  92 variables:, all vars numeric
#OLD                   ## Create rownames that are UIDs
#OLD                     rownames(coverB_filt_flat)<-names
#OLD                 
#OLD         ## Remove columns that aren't species to conform to input requirements for metaMDS
#OLD       coverB_filt_flat_forMDS<-coverB_filt_flat
#OLD         coverB_filt_flat_forMDS[is.na(coverB_filt_flat_forMDS)==T]<-0
#OLD       #numeric
#OLD           coverB_filt_flat_forMDS$rowsums<-rowSums(coverB_filt_flat_forMDS)
#OLD               coverB_filt_flat_forMDS[is.na(coverB_filt_flat_forMDS)==T]<-0
#OLD                 #coverB_filt_flat_forMDS<-mapply(as.numeric, coverB_filt_flat_forMDS) %>% as.data.frame()
#OLD                   #coverB_filt_flat_forMDS[is.na(coverB_filt_flat_forMDS)==T]<-0
#OLD                     coverB_filt_flat_forMDS<-rename(coverB_filt_flat_forMDS, tot_cov =rowsums)
#OLD                       coverB_filt_flat_forMDS<-coverB_filt_flat_forMDS %>% subset(tot_cov>0) %>% dplyr::select(-tot_cov)#subset(tot_cov>0)# %>% select(-tot_cov)
#OLD                         names_nmds<-rownames(coverB_filt_flat_forMDS)
#OLD ## Run NMDS with relative distance (Bray) measure
#OLD coverB_MDS<-metaMDS(coverB_filt_flat_forMDS, distance = "bray", try=100, trymax = 500)
#OLD plot(coverB_MDS)

# Make the a column in the NMDS output that is the quadrat UID to be able to join with the env data 
#Need to filter Env_matrix to only have the 49 rows but now has rownames
stuff<-coverB_MDS$points %>% as.data.frame() %>% mutate(UID= names_nmds)  %>% replace(is.na(.),0)
stuff_Trent<-coverB_Trent_MDS$points %>% as.data.frame() %>% mutate(UID= names_nmds_Trent$UID) %>% replace(is.na(.),0)
stuff_Treb<-coverB_Treb_MDS$points %>% as.data.frame() %>% mutate(UID= names_nmds_Treb$UID) %>% replace(is.na(.),0)
  stuff_Treb$UID<-as.character(stuff_Treb$UID)
## Join the NMDS data with the 
  str(Env_matrix_B)
  str(stuff_Treb)
stuff_env<-Env_matrix_B %>% inner_join(stuff, by="UID")  %>% as.data.frame() # %>% select(VPD_min)
stuff_env_Trent<-Env_matrix_B %>% inner_join(stuff_Trent, by="UID")  %>% as.data.frame() # %>% select(VPD_min)
stuff_env_Treb<-Env_matrix_B %>% inner_join(stuff_Treb, by="UID")  %>% as.data.frame() # %>% select(VPD_min)

## Add the traits to the env matrix
stuff_env<-
  stuff_env %>% 
  dplyr::select(-X, -X.1) %>% #colnames()
    inner_join(trait, by="UID", keep=FALSE)  %>% #colnames()
    #inner_join(trait, by=c("UID","Elevation","Transect","Quadrat"), keep=F)  %>% #colnames()
    #inner_join(substrate, by="UID", keep=FALSE) %>% colnames()
    inner_join(substrate, by=c("UID","Elevation","Quadrat"), keep=FALSE) %>%
     as.data.frame() %>% 
    dplyr::select(UID, Transect, Quadrat, Elevation, -Total, everything())

stuff_env_Trent<-stuff_env_Trent %>% 
  dplyr::select(-X, -X.1) %>% #colnames()
  inner_join(trait, by="UID", keep=FALSE)  %>% #colnames()
  #inner_join(trait, by=c("UID","Elevation","Transect","Quadrat"), keep=F)  %>% #colnames()
  #inner_join(substrate, by="UID", keep=FALSE) %>% colnames()
  inner_join(substrate, by=c("UID","Elevation","Quadrat"), keep=FALSE) %>%
  as.data.frame() %>% 
  dplyr::select(UID, Transect, Quadrat, Elevation, -Total, everything())

stuff_env_Treb<-stuff_env_Treb %>% 
  dplyr::select(-X, -X.1) %>% #colnames()
  inner_join(trait, by="UID", keep=FALSE)  %>% #colnames()
  #inner_join(trait, by=c("UID","Elevation","Transect","Quadrat"), keep=F)  %>% #colnames()
  #inner_join(substrate, by="UID", keep=FALSE) %>% colnames()
  inner_join(substrate, by=c("UID","Elevation","Quadrat"), keep=FALSE) %>%
  as.data.frame() %>% 
  dplyr::select(UID, Transect, Quadrat, Elevation, -Total, everything())


## Make a list of env variable names to use in functions later
vars<-colnames(stuff_env[5:length(stuff_env)-1])  #%>% #as.data.frame() %>% dplyr::select(-Total)
vars_Treb<-colnames(stuff_env_Treb[5:length(stuff_env_Treb)-1])#  %>% replace(is.na(.),0) %>% select(-Total)
vars_Trent<-colnames(stuff_env_Trent[5:length(stuff_env_Trent)-1])#  %>% replace(is.na(.),0) %>% select(-Total)

##PASS: Unit test of using the GAM model of Env var ~ NMS Axis 1 + Axis 2 and returns the r.sq
  stuff_out1<-ordisurf(coverB_MDS~Inclination, stuff_env, main= paste(colnames(stuff_env[1]))) %>% summary() #%>% select(r.sq)
  stuff_out1_Trent<-ordisurf(coverB_Trent_MDS~Inclination, stuff_env_Trent, main= paste(colnames(stuff_env_Trent[5]))) %>% summary() #%>% select(r.sq)
  stuff_out1_Treb<-ordisurf(coverB_Treb_MDS~Inclination, stuff_env_Treb, main= paste(colnames(stuff_env_Treb[5]))) %>% summary() #%>% select(r.sq)
  
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
      #out1<-ordisurf(eval(parse(text=paste(MDS, "~",vars[x],sep=""))), stuff_env) %>% summary()
      out1<-ordisurf(eval(parse(text=paste("coverB_MDS~",vars[x],sep=""))), stuff_env) %>% summary()
      return(round(out1$r.sq,2))
      return(round(out1$se,2))
      }
        #ordi_stats_out<-lapply(1:length(vars),ordi_stats(MDS=coverB_MDS)) %>% unlist()
        ordi_stats_out<-lapply(1:length(vars),ordi_stats) %>% unlist()
        
          ordi_stats_out<-cbind(vars,ordi_stats_out) %>% as.data.frame()
            colnames(ordi_stats_out)<-c("Env Variable","R squared")
              ##Now write.csv doesn't work!?
              write_csv(ordi_stats_out, "Patached_TransectB_CommunityEnv_NMS_GAM_fit.csv")
   
       ##Stats for Treb 
        ordi_stats_Treb<-function(x)
        {
          #out1<-ordisurf(eval(parse(text=paste(MDS, "~",vars[x],sep=""))), stuff_env) %>% summary()
          out1<-ordisurf(eval(parse(text=paste("coverB_Treb_MDS~",vars[x],sep=""))), stuff_env_Treb) %>% summary()
          return(round(out1$r.sq,2))
          return(round(out1$se,2))
        }
        #ordi_stats_out<-lapply(1:length(vars),ordi_stats(MDS=coverB_MDS)) %>% unlist()
        ordi_stats_out_Treb<-lapply(1:length(vars_Treb),ordi_stats_Treb) %>% unlist()
        
        ordi_stats_out_Treb<-cbind(vars,ordi_stats_out_Treb) %>% as.data.frame()
        colnames(ordi_stats_out_Treb)<-c("Env Variable","R squared")
        ##Now write.csv doesn't work!?
        write_csv(ordi_stats_out_Treb, "Patached_TransectB_CommunityEnv_NMS_GAM_fit_Treb.csv")
       
        ##Stats for Trent 
        ordi_stats_Trent<-function(x)
        {
          #out1<-ordisurf(eval(parse(text=paste(MDS, "~",vars[x],sep=""))), stuff_env) %>% summary()
          out1<-ordisurf(eval(parse(text=paste("coverB_Trent_MDS~",vars[x],sep=""))), stuff_env_Trent) %>% summary()
          return(round(out1$r.sq,2))
          return(round(out1$se,2))
        }
        #ordi_stats_out<-lapply(1:length(vars),ordi_stats(MDS=coverB_MDS)) %>% unlist()
        ordi_stats_out_Trent<-lapply(1:length(vars_Trent),ordi_stats_Trent) %>% unlist()
        
        ordi_stats_out_Trent<-cbind(vars,ordi_stats_out_Trent) %>% as.data.frame()
        colnames(ordi_stats_out_Trent)<-c("Env Variable","R squared")
        ##Now write.csv doesn't work!?
        write_csv(ordi_stats_out_Trent, "Patached_TransectB_CommunityEnv_NMS_GAM_fit_Trent.csv")
        
              
              
#Make hilltop plots with R2 
    ordi_fit<-function(x) 
            {
      ##jpeg(paste("Patache_Hilltop",vars[x],".jpeg"))
      ordisurf(eval(parse(text=paste("coverB_MDS~",vars[1],sep=""))), stuff_env, main= "",labcex=0, col='black')  
        ordisurf(eval(parse(text=paste("coverB_MDS~",vars[x],sep=""))), stuff_env, main= paste(vars[x]),labcex=1, add=T)
          text(max(stuff$MDS1)*0.5, max(stuff$MDS2)*0.95, paste("R2=",ordi_stats_out[x,2]), cex=2)
           #elev<-envfit(eval(parse(text=paste("coverB_MDS~",vars[2],sep=""))), stuff_env, main= paste(vars[2]))
           # plot(elev, cex=1.5)
      #dev.off()
          
            }
    
    pdf("Patache_Hilltop_Plots.pdf")
    lapply(1:length(vars),ordi_fit)
    dev.off()
    
    
    
    ######Hilltop plots for Trebouxioud lichens
    
    ordi_fit_Treb<-function(x) 
    {
      ##jpeg(paste("Patache_Hilltop",vars_Treb[x],".jpeg"))
      ordisurf(eval(parse(text=paste("coverB_Treb_MDS~",vars_Treb[1],sep=""))), stuff_env_Treb, main= "",labcex=0, col='black')  
      ordisurf(eval(parse(text=paste("coverB_Treb_MDS~",vars_Treb[x],sep=""))), stuff_env_Treb, main= paste(vars_Treb[x]),labcex=1, add=T)
      text(max(stuff_Treb$MDS1)*0.5, max(stuff_Treb$MDS2)*0.95, paste("R2=",ordi_stats_out_Treb[x,2]), cex=2)
      #elev<-envfit(eval(parse(text=paste("coverB_Treb_MDS~",vars_Treb[2],sep=""))), stuff_env_Treb, main= paste(vars_Treb[2]))
      # plot(elev, cex=1.5)
      #dev.off()
      
    }
    
    pdf("Patache_Hilltop_Plots_Treb.pdf")
    lapply(1:length(vars_Treb),ordi_fit_Treb)
    dev.off()
    
    
    
    
    ######Hilltop plots for Trentepohlioid lichens
    
    ordi_fit_Trent<-function(x) 
    {
      ##jpeg(paste("Patache_Hilltop",vars_Trent[x],".jpeg"))
      ordisurf(eval(parse(text=paste("coverB_Trent_MDS~",vars_Trent[1],sep=""))), stuff_env_Trent, main= "",labcex=0, col='black')  
      ordisurf(eval(parse(text=paste("coverB_Trent_MDS~",vars_Trent[x],sep=""))), stuff_env_Trent, main= paste(vars_Trent[x]),labcex=1, add=T)
      text(max(stuff_Trent$MDS1)*0.5, max(stuff_Trent$MDS2)*0.95, paste("R2=",ordi_stats_out_Trent[x,2]), cex=2)
      #elev<-envfit(eval(parse(text=paste("coverB_Trent_MDS~",vars_Trent[2],sep=""))), stuff_env_Trent, main= paste(vars_Trent[2]))
      # plot(elev, cex=1.5)
      #dev.off()
      
    }
    
    pdf("Patache_Hilltop_Plots_Trent.pdf")
    lapply(1:length(vars_Trent),ordi_fit_Trent)
    dev.off()
    
    
######### Second matrix ... Daniel recommeded starting looking at microclimate, specifically
######### sat85_dry, Tmed, VPDmed ... where are these columns?
## Find files with data from same quadrats/elevations on transect B
## anat, he, he_red, quad, rou, rouVar all are filterable by transect
## angles, dry only has data for transect B
## unclear which transect rock and rock_st came from
    stuff_env[,-1:-4] %>% cor(use="pairwise.complete.obs") %>% cov2cor() %>% view()