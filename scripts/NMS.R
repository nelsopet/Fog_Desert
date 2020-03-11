
require(tidyverse)
require(vegan)


## Two objects need to be read in
cover<-read.csv('data/CoverMaster.csv')
#cover$UID<-paste(cover$Elevation, cover$Quadrat, sep="_")

#str(cover)
#'data.frame':	1522 obs. of  8 variables:
#  $ Elevation : int  450 450 450 450 450 450 450 450 450 450 ...
#$ Transect  : Factor w/ 3 levels "","A","B": 2 2 2 2 2 2 2 2 2 2 ...
#$ Quadrat   : Factor w/ 36 levels "","1","10","10b",..: 2 20 23 25 27 27 27 29 31 33 ...
#$ Species   : Factor w/ 158 levels "","AcarosporaBull",..: 1 38 1 1 38 114 152 38 1 1 ...
#$ Cover     : Factor w/ 17 levels "","1","10","11",..: 1 8 1 1 8 15 17 17 1 1 ...
#$ Collection: Factor w/ 204 levels "","no coll","T450_A2.1",..: 1 3 1 1 3 4 5 3 1 1 ...
#$ Long_Form : Factor w/ 94 levels "","Acarospora bullate",..: 1 20 1 1 1 1 1 1 1 1 ...
#$ Photos    : Factor w/ 2 levels "","DSC3979": 1 1 1 1 1 1 1 1 1 1 ...

##where is the Env data read in?
Env_matrix_B<-read.csv('data/Env_matrix_B.csv')

## Read in trait data
trait<-read.csv('data/QuadratMeanTraits_Weighted_NoTrace.csv')
trait$UID<-paste(trait$Elevation, trait$Quadrat, sep="_")
trait<-trait[,-1:-4]
  #'data.frame':	49 obs. of  11 variables:
  #  $ X              : int  1 2 3 4 5 6 7 8 9 10 ...
  #$ Transect       : Factor w/ 1 level "B": 1 1 1 1 1 1 1 1 1 1 ...
  #$ Elevation      : int  300 400 500 600 700 800 400 500 600 700 ...
  #$ Quadrat        : int  1 1 1 1 1 1 10 10 10 10 ...
  #$ Time_Surface   : num  1 1.16 13.67 10.18 25.83 ...
  #$ Time_Absorption: num  7.17 23.89 26.58 32.23 65.74 ...
  #$ Angle          : num  130 97.5 103.3 93.8 92.8 ...
  #$ .N             : num  NA 1.18 0.95 1.76 1.37 ...
  #$ d15N           : num  NA 1.71 -1.24 -0.215 -0.529 ...
  #$ .C             : num  NA 40.7 28.8 28.8 30.5 ...
  #$ d13C           : num  NA -24.4 -28.6 -21.3 -24 ... 
#test to make sure all species in cover are in trait df
#976 rows without 
anti_join(cover, trait, by="UID") %>% nrow()
nrow(Env_matrix_B)
# 69 species with no trait data
anti_join(cover, trait, by="Species") %>% select(Species) %>% unique() %>% nrow()
#2              Caloplaca2
#6             Roccellina1
#22          ArthoniaRocc1
#23         LightBrownRocc
#25            Roccellina2
#28               Buellia3
#30              Ramalina2
#36               Buellia2
#37     YellowRhizocarpon2
#38            BuelliaBlue
#40             Caloplaca3
#41             GrayCrust1
#44               Buellia4
#55             BrownLumpy
#56   Caloplacapergracilis
#57                Niebla3
#58              WhiteFoam
#59             Rhizoplaca
#60          Sclerophyton1
#61         Candelariella1
#64    Roccellinaterricola
#65  TerricolousWhiteCrust
#68               Buellia1
#70         DarkBrownSmear
#74            Amigdalaria
#75   Caloplacataltalensis
#78              Lepraria1
#80              Xanthoria
#83    YellowConvexBuellia
#88               Ramalina
#94    DarkGrayRhizocarpon
#108        WhiteRhizoPara
#110        Candelariella3
#123   TinyGreyRhizocarpon
#135  ClumpyMustardBuellia
#137    BrownConvexBuellia
#138 Caloplacasubhervidela
#140          Chrysothrix2
#150        Candelariella2
#152         Heterodermia1
#208               Niebla1
#219            BlackRhizo
#231          Chrysothrix1
#245         Heterodermia2
#246               Niebla4
#254          GrayLepraria
#260        BrownAspicilia
#308         GreenLepraria
#339               Niebla2
#360          GrayIsidiate
#425         WhiteBuellia2
#583         Cyanobacteria
#596      GraySorAspicilia
#607     RoccellaFruticose
#611         BrownIsidiate
#620              BrownSor
#623            Caloplaca1
#650         Sclerophyton2
#652          YellowCrust1
#684           BrownCrust1
#760            Aspicilia1
#826       PaleRhizocarpon
#833         WhiteSorAreol
#839          Pyrenodesmia
#854       RoccellinaCereb
#860    GraySorRhizocarpon
#876          Nieblavagans
#879            Endolithic

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
  dev.off()
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
  select(UID, X, Transect, Quadrat, Elevation, everything())
colnames(stuff_env)

## Make a list of env variable names to use in functions later
vars<-colnames(stuff_env[4:length(stuff_env)])  #%>% as.data.frame() 




##PASS: Unit test of using the GAM model of Env var ~ NMS Axis 1 + Axis 2 and returns the r.sq
  #stuff_out1<-ordisurf(coverB_MDS~Inclination, stuff_env, main= paste(colnames(stuff_env[5]))) %>% summary() %>% select(r.sq)
  #stuff_out1$r.sq
  #dev.off()
  
  
## automate ordisurf for all columns
## Need to unlist the output of each ordisurf ... str(stuff_out) resultsing a list of 54
## Error indicates input dimensions don't match

  tst_hill<-ordisurf(eval(parse(text=paste("coverB_MDS~",vars[25],sep=""))), stuff_env)
  tst_hill_stats<-summary(tst_hill) 
  tst_hill_stats$r.sq
  dev.off()

  tst_hill<-ordisurf(eval(parse(text=paste("coverB_MDS~",vars[4],sep=""))), stuff_env)
  tst_hill_stats<-summary(tst_hill) 
  tst_hill_stats$r.sq
  str(tst_hill_stats)
  dev.off()
  
  
    
  make_mods<-function (x) {paste("coverB_MDS~",vars[x],sep="") %>% as.formula()}
  var_mods<-lapply(1:length(vars),make_mods)
  
  eval(var_mods[1])
  
    tst1<-envfit(eval(parse(text=paste("coverB_MDS~",vars[2],sep=""))), stuff_env, main= paste(vars[2]))
    tst
    dev.off()
    tst<-ordisurf(eval(parse(text=paste("coverB_MDS~",vars[2],sep=""))), stuff_env, main= paste(vars[2]))
    
    ordi_fit<-function(x) {ordisurf(eval(parse(text=paste("coverB_MDS~",vars[x],sep=""))), stuff_env, main= paste(vars[x]))
      elev<-envfit(eval(parse(text=paste("coverB_MDS~",vars[2],sep=""))), stuff_env, main= paste(vars[2]))
      plot(elev)}
    
    
    length(stuff_env)
    
    ## Why does this produce a vector of numbers read as characters?
    #vars<-as.character(vars)
    pdf("Patache_Hilltop_Plots.pdf")
    lapply(1:length(vars),ordi_fit)
    dev.off()
    
    ordi_stats<-function(x)
    {
    out1<-ordisurf(eval(parse(text=paste("coverB_MDS~",vars[x],sep=""))), stuff_env) %>% summary()
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
