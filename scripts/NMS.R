##### Required packages
#install.packages('devtools')
#Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")

#devtools::install_github('phytomosaic/foggy')
#devtools::install_github('phytomosaic/ecole')

require(tidyverse)
require(vegan)
require(Metrics)
require(hsdar)
set.seed(1234)
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
  
## automate ordisurf for all columns
## Need to unlist the output of each ordisurf ... str(stuff_out) resultsing a list of 54
  ##Unit test for hilstop plot and  grabbing R2 PASSES
  tst_hill<-ordisurf(eval(parse(text=paste("coverB_MDS~",vars[2],sep=""))), stuff_env)
  tst_hill_stats<-summary(tst_hill) 
  tst_hill_stats$r.sq
  #dev.off()

  ##Function to make models ... not sure if it even still get used though
  make_mods<-function (x) {paste("coverB_MDS~",vars[x],sep="") %>% as.formula()}
  var_mods<-lapply(1:length(vars),make_mods)
    eval(var_mods[1])
    ##Unit test for making linear biplot vector for elevation PASSSES
    tst1<-envfit(eval(parse(text=paste("coverB_MDS~",vars[2],sep=""))), stuff_env, main= paste(vars[2]))
    plot(tst1)
    dev.off()
 
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
              #write_csv(ordi_stats_out, "Patached_TransectB_CommunityEnv_NMS_GAM_fit.csv")
   
       ##Stats for Treb 
              ##Unit test for hilstop plot and  grabbing R2 
              ordisurf(eval(parse(text=paste("coverB_Treb_MDS~",vars_Treb[1],sep=""))), stuff_env_Treb)  
              ordisurf(eval(parse(text=paste("coverB_Treb_MDS~",vars_Treb[8],sep=""))), stuff_env_Treb)
              tst_hill<-ordisurf(eval(parse(text=paste("coverB_Treb_MDS~",vars[8],sep=""))), stuff_env_Treb)
              tst_hill_stats<-summary(tst_hill) 
              tst_hill_stats$r.sq
              dev.off()
              
              
              tst<-lm(stuff_env_Treb$Dry3kPa_Dry~stuff_Treb$MDS2+stuff_Treb$MDS1)
              summary(tst)
              
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
        #write_csv(ordi_stats_out_Treb, "Patached_TransectB_CommunityEnv_NMS_GAM_fit_Treb.csv")
       
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
        #write_csv(ordi_stats_out_Trent, "Patached_TransectB_CommunityEnv_NMS_GAM_fit_Trent.csv")
        
              
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
    
###############Combine stats from full, Trent and Treb ordinations 
reject = c("MDS1","MD2", ".C", "d13C", ".N") %>% as.data.frame()
colnames(reject) <- "Env Variable"
ordi_stats_out %>% 
      rename(All_R2 = 'R squared') %>%
        inner_join(ordi_stats_out_Treb, by="Env Variable") %>% 
          rename(Treb_R2 = 'R squared') %>%
            inner_join(ordi_stats_out_Trent, by="Env Variable") %>%
              rename(Trent_R2 = 'R squared') %>%
                anti_join(reject, by="Env Variable") %>%
                  arrange(desc(All_R2), desc(Treb_R2), desc(Trent_R2)) %>%
                    dplyr::filter(All_R2 > 0.3) %>%
                     write_csv("Patache_TransectB_CommunityEnv_NMS_GAM.csv")
    