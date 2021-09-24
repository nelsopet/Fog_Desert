######################################################################
#
#  NMS and hilltop plots
#
#    Peter Nelson, w changes from Rob Smith, phytomosaic@gmail.com, 16 Sep 2021
#
##      GNU General Public License, Version 3.0    ###################


### preamble
rm(list=ls())
set.seed(1234)
require(tidyverse)
require(vegan)
require(mgcv)
devtools::install_github("https://github.com/phytomosaic/ecole")


### define functions
# quick function for proportion of no-share sample units
`noshare` <- function(x) {
  z <- vegan::no.shared(x)
  length(z[z == TRUE])/length(z)
}
# step-down dimensionality selection in NMS, for varying 'k' dimensions
`scree_nms` <- function(D, k=9, ...) {
  stress <- rep(NA, k)
  for (i in 1:k) {
    cat('calculating', i, 'of', k, 'dimensions...\n')
    stress[i] <- metaMDS(D, k=i, trace=0, ...)$stress
  }
  plot(1:k, stress, main='', xlab='Dimension', ylab='Stress',
       ylim=c(0, max(stress)*1.05), pch=16, las=1, bty='l')
  lines(1:k, stress)
  abline(0.20, 0, col='red', lty = 2)
  data.matrix(stress)
}


### read in data
# cover abundance
cover     <- read.csv('data/CoverMaster.csv')
cover$UID <- paste(cover$Elevation, cover$Quadrat, sep='_')
# species-level traits, including absorption time and angles
traits_sp <- read.csv('data/Traits_TransB.csv')
# environment data, including climate derivatives
env       <- read.csv('data/Env_matrix_B.csv')
env$UID   <- paste(env$Elevation, env$Quadrat, sep='_')
# trait data, including spectra, absorption time, photobionts
# trait <- read.csv('data/QuadratMeanTraits_Weighted_NoTrace.csv')
trait     <- read.csv('data/QuadratMeanTraits_Unweighted.csv')
trait$UID <- paste(trait$Elevation, trait$Quadrat, sep='_')
trait     <- trait[,-1:-4]
# spectra by quadrat, including band ratios
spectra   <- read.csv('data/Spectral_Features.csv')
# substrate data, including quadrat-level rockiness measurements
#   QUESTION: Join this to traits?
substrate     <- read.csv('data/Substrate_PatacheB.csv',header=T)
substrate$UID <- paste(substrate$Elevation, substrate$Quadrat, sep='_')


### examine
str(cover)     # 1522 obs. of  9 variables
str(traits_sp) # 356 obs. of  26 variables
str(env)       # 60 obs. of  30 variables
str(trait)     # 49 obs. of  127 variables
str(spectra)   # 49 obs. of  127 variables
str(substrate) # 60 obs. of  10 variables


### --- data pre-processing (RJS assumes this works as intended, not checked)
# use only transect B, which has complete community and environmental data
coverB       <- subset(cover, Transect=='B')
coverB$Cover <- as.character(coverB$Cover)
coverB$Cover[coverB$Cover=='T'] <- 0.25
biont        <- c('Trebouxia','Trentepohlia')
table(traits_sp$Photobiont, useNA='always')
traits_sp_filt <- traits_sp %>% dplyr::select(Species, Photobiont) %>% 
  unique() %>% as.data.frame()
# all matrix
coverB_flat_forMDS <-
  coverB %>% 
  dplyr::filter(Species != '') %>%
  dplyr::select(-Collection, -Long_Form, -Photos) %>% 
  pivot_wider(names_from=Species, values_from=Cover) %>%
  dplyr::select(-Elevation, -Transect, -Quadrat) %>%
  mutate(across(c(-UID),as.numeric)) %>%
  replace(is.na(.),0) %>%
  rowwise() %>% 
  mutate(tot_cov=sum(across(c(-UID)))) %>% 
  dplyr::filter(tot_cov>0) %>%
  dplyr::select(-tot_cov)
dimnames(coverB_flat_forMDS)
coverB_flat_forMDS <- as.data.frame(coverB_flat_forMDS)
dimnames(coverB_flat_forMDS)[[1]] <- coverB_flat_forMDS$UID
coverB_flat_forMDS$UID <- NULL
# Trebouxia matrix
coverB_Treb <- coverB %>%  left_join(traits_sp_filt, by='Species', keep=F) %>% 
  subset(Photobiont == 'Trebouxia')
coverB_Treb_flat_forMDS <- 
  coverB_Treb %>%
  dplyr::select(Elevation, Transect, Quadrat, UID, Species, Cover) %>%
  dplyr::filter(Species!='') %>%
  pivot_wider(names_from=Species, values_from=Cover) %>%
  dplyr::select(-Elevation, -Transect, -Quadrat) %>%
  mutate(across(c(-UID),as.numeric)) %>%
  replace(is.na(.),0) %>%
  rowwise() %>% 
  mutate(tot_cov=sum(across(c(-UID)))) %>% 
  dplyr::filter(tot_cov>0) %>%
  dplyr::select(-tot_cov,-UID)
names_nmds_Treb <- 
  coverB_Treb %>% 
  dplyr::select(Elevation, Transect, Quadrat, UID, Species, Cover) %>%
  dplyr::filter(Species!='') %>%
  pivot_wider(names_from=Species, values_from=Cover) %>%
  dplyr::select(-Elevation, -Transect, -Quadrat) %>%
  mutate(across(c(-UID),as.numeric)) %>%
  replace(is.na(.),0) %>%
  rowwise() %>% 
  mutate(tot_cov=sum(across(c(-UID)))) %>% 
  dplyr::filter(tot_cov>0) %>%
  dplyr::select(UID)
# Trentepohlia matrix
coverB_Trent <- coverB %>% left_join(traits_sp_filt, by='Species', keep=F) %>%
  subset(Photobiont == 'Trentepohlia')
coverB_Trent_flat_forMDS <- 
  coverB_Trent %>%
  dplyr::select(Elevation, Transect, Quadrat, UID, Species, Cover) %>%
  dplyr::filter(Species!='') %>%
  pivot_wider(names_from=Species, values_from=Cover) %>%
  dplyr::select(-Elevation, -Transect, -Quadrat) %>%
  mutate(across(c(-UID),as.numeric)) %>%
  replace(is.na(.),0) %>%
  rowwise() %>% 
  mutate(tot_cov=sum(across(c(-UID)))) %>% 
  dplyr::filter(tot_cov>0) %>%
  dplyr::select(-tot_cov,-UID)
names_nmds_Trent <-
  coverB_Trent %>%
  dplyr::select(Elevation, Transect, Quadrat, UID, Species, Cover) %>%
  dplyr::filter(Species!='') %>%
  pivot_wider(names_from=Species, values_from=Cover) %>%
  dplyr::select(-Elevation, -Transect, -Quadrat) %>%
  mutate(across(c(-UID),as.numeric)) %>%
  replace(is.na(.),0) %>%
  rowwise() %>% 
  mutate(tot_cov=sum(across(c(-UID)))) %>% 
  dplyr::filter(tot_cov>0) %>%
  dplyr::select(UID)


### double-standardization by species max then site totals
coverB_flat_forMDS        <- wisconsin(coverB_flat_forMDS)
coverB_Treb_flat_forMDS   <- wisconsin(coverB_Treb_flat_forMDS)
coverB_Trent_flat_forMDS  <- wisconsin(coverB_Trent_flat_forMDS)
# # remove rare species?
# table(colSums(coverB_flat_forMDS > 0)) # 38 of 89 taxa are singletons
# table(colSums(coverB_Treb_flat_forMDS > 0)) # 22 of 60 taxa are singletons
# table(colSums(coverB_Trent_flat_forMDS > 0)) # 13 of 25 taxa are singletons


### dissimilarity matrices
noshare(coverB_flat_forMDS)       # 0.346 no-share, advise stepacross
noshare(coverB_Treb_flat_forMDS)  # 0.480 no-share, advise stepacross
noshare(coverB_Trent_flat_forMDS) # 0.454 no-share, advise stepacross
D_all   <- stepacross(vegdist(coverB_flat_forMDS, 'bray'))
D_treb  <- stepacross(vegdist(coverB_Treb_flat_forMDS, 'bray'))
D_trent <- stepacross(vegdist(coverB_Trent_flat_forMDS, 'bray'))


### NMS runs
par(mfrow=c(1,3))
# NMS all
scree_nms(D_all, k=5, try=100)             # dimensionality selection
coverB_MDS <- metaMDS(D_all, k=2, try=500, trymax=501, auto=F, noshare=F)
stressplot(coverB_MDS)                     # shepard plot
plot(coverB_MDS, type='t', cex=0.6)        # NMS plot
# NMS Trebouxia
scree_nms(D_treb, k=5, try=100)            # dimensionality selection
coverB_Treb_MDS <- metaMDS(D_treb, k=2, try=500, trymax=501, auto=F, noshare=F)
stressplot(coverB_Treb_MDS)                # shepard plot
plot(coverB_Treb_MDS, type='t', cex=0.6)   # NMS plot
# NMS Trentepohlia
scree_nms(D_trent, k=5, try=100)           # dimensionality selection
coverB_Trent_MDS <- metaMDS(D_trent, k=2, try=500, trymax=501, auto=F, noshare=F)
stressplot(coverB_Trent_MDS)               # shepard plot
plot(coverB_Trent_MDS, type='t', cex=0.6)  # CAUTION! many ties, and wonky Axis 1...


### join NMS scores/env/traits
names_nmds  <- dimnames(coverB_flat_forMDS)[[1]]
scr       <- coverB_MDS$points %>% as.data.frame() %>% 
  mutate(UID=names_nmds) %>% replace(is.na(.),0)
scr_Trent <- coverB_Trent_MDS$points %>% as.data.frame() %>% 
  mutate(UID=names_nmds_Trent$UID) %>% replace(is.na(.),0)
scr_Treb  <- coverB_Treb_MDS$points %>% as.data.frame() %>% 
  mutate(UID=names_nmds_Treb$UID) %>% replace(is.na(.),0)
# add NMS scores to env
scr_env       <- env %>% inner_join(scr, by='UID') %>% as.data.frame() 
scr_env_Trent <- env %>% inner_join(scr_Trent, by='UID')  %>% as.data.frame() 
scr_env_Treb  <- env %>% inner_join(scr_Treb, by='UID')  %>% as.data.frame() 
# add traits to the env matrix
scr_env <-
  scr_env %>% 
  dplyr::select(-X, -X.1) %>%
  inner_join(trait, by='UID', keep=FALSE)  %>%
  inner_join(substrate, by=c('UID','Elevation','Quadrat'), keep=FALSE) %>%
  as.data.frame() %>% 
  dplyr::select(UID, Transect, Quadrat, Elevation, -Total, everything())
scr_env_Treb <- scr_env_Treb %>% 
  dplyr::select(-X, -X.1) %>%
  inner_join(trait, by='UID', keep=FALSE)  %>%
  inner_join(substrate, by=c('UID','Elevation','Quadrat'), keep=FALSE) %>%
  as.data.frame() %>% 
  dplyr::select(UID, Transect, Quadrat, Elevation, -Total, everything())
scr_env_Trent <- scr_env_Trent %>% 
  dplyr::select(-X, -X.1) %>%
  inner_join(trait, by='UID', keep=FALSE)  %>%
  inner_join(substrate, by=c('UID','Elevation','Quadrat'), keep=FALSE) %>%
  as.data.frame() %>% 
  dplyr::select(UID, Transect, Quadrat, Elevation, -Total, everything())
# env variable names to use later
vars <- colnames(scr_env[5:length(scr_env)-1]) 
vars[vars=='UID'] <- NULL
# cleanup
rm(substrate, cover, traits_sp, env, trait, spectra,
   names_nmds, names_nmds_Treb, names_nmds_Trent,
   D_all, D_treb, D_trent)


### CAUTION!!! incomplete replication within Elevations...
cbind(xtabs(~ Elevation + Quadrat, scr_env), Total=table(scr_env$Elevation))
cbind(xtabs(~ Elevation + Quadrat, scr_env_Treb), Total=table(scr_env_Treb$Elevation))
cbind(xtabs(~ Elevation + Quadrat, scr_env_Trent), Total=table(scr_env_Trent$Elevation))


### plot NMS by elevation
`plot_nms` <- function(d, ...) {
  plot(d$MDS1, d$MDS2, xlab='NMS 1', ylab='NMS 2',
       pch=as.numeric(factor(d$Elevation))+14, 
       col=as.numeric(factor(d$Elevation)), 
       bty='L', las=1, asp=1, main='', ...)
}
par(mfrow=c(1,3), pty='s')
plot_nms(scr_env)
plot_nms(scr_env_Treb)
plot_nms(scr_env_Trent)
# cant do mixed model because Elevations occupy diff regions of ordination space


### fit GAM models (aka hilltop plots...)
`fit_gam` <- function(y, d, saveplot=FALSE, ...) {
  o  <- ordisurf(scores(d[,c('MDS1','MDS2')]), d[,y], plot=F)
  r2 <- round(summary(o)$r.sq,3)
  if(isTRUE(saveplot)) { png(paste0('./fig/fig_surface_',y,'.png'), 
                             hei=10, wid=10, unit='in', res=350) }
  plot(d$MDS1, d$MDS2, xlab='NMS 1', ylab='NMS 2',
       pch=21, cex=1.0, col='#00000099' , #bg=ecole::colvec(d$Elevation), 
       bty='L', las=1, asp=1, main=paste0(y,', R2=',r2), ...)
  contour(o$grid$x, o$grid$y, o$grid$z, col='#000000FF', add=T)
  if(isTRUE(saveplot)) dev.off()
  return(o)
}


### for testing only: reduce number of candidate responses
length(vars) # 155 candidate responses (!)
vars <- vars[c(1,2,3,4, 5,7,8,9, 28,29,30,37, 51,141,151,153)]


### fit the GAMs
png("output/Patache_Hilltops.png", height = 5000, width = 5000, res = 350)
par(mfrow=c(4,4), mar=c(3,3,1,0), oma=c(0,0,0,0), pty='s')
m_all   <- lapply(vars, function(i) fit_gam(i, d=scr_env))
names(m_all)   <- vars
dev.off()

png("output/Patache_Hilltops_Treb.png", height = 5000, width = 5000, res = 350)
par(mfrow=c(4,4), mar=c(3,3,1,0), oma=c(0,0,0,0), pty='s')
m_treb  <- lapply(vars, function(i) fit_gam(i, d=scr_env_Treb))
names(m_treb)  <- vars
dev.off()

png("output/Patache_Hilltops_Trent.png", height = 5000, width = 5000, res = 350)
par(mfrow=c(4,4), mar=c(3,3,1,0), oma=c(0,0,0,0), pty='s')
m_trent <- lapply(vars, function(i) fit_gam(i, d=scr_env_Trent))
names(m_trent) <- vars
dev.off()

### GAM goodness-of-fit
`get_r2` <- function(mod_lst, savecsv=FALSE, prefix='all') {
  r2        <- sapply(mod_lst, function(i) summary(i)$r.sq)
  names(r2) <- names(mod_lst)
  if(isTRUE(savecsv)) {
    write.csv(t(t(r2)), file=paste0('./', prefix, '_PatacheB_nmsgamfit.csv'))
  }
  return(r2)
}
r2_all   <- get_r2(m_all)
r2_treb  <- get_r2(m_treb)
r2_trent <- get_r2(m_trent)
(r2 <- data.frame(r2_all, r2_treb, r2_trent))
### view R2 across photobiont types
`plot_r2_matrix` <- function(r2, ...) {
  m <- round(as.matrix(r2),2)
  m <- m[NROW(m):1,]
  par(pty='s', mfrow=c(1,1))
  image(1:ncol(m), 1:nrow(m), t(m), col=viridis::viridis(99), axes=F, ann=F)
  axis(1, 1:ncol(m), colnames(m), las=1)
  axis(2, 1:nrow(m), rownames(m), las=1)
  for (x in 1:ncol(m)) { for (y in 1:nrow(m)) { text(x,y,m[y,x]) }}
}
plot_r2_matrix(r2)


####    END    ####