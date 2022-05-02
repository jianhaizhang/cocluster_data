# Read validation results.
auc.arab.rt.sc9 <- readRDS('validate_opt_res/auc.dat.fct.fil1.arab.rt.sc9.pca.knn.rds')
auc.arab.rt.sc51 <- readRDS('validate_opt_res/auc.dat.fct.fil1.arab.rt.sc51.pca.knn.rds') 
auc.mus.brain <- readRDS('validate_opt_res/auc.dat.fct.fil1.mus.brain.pca.knn.rds')
auc.mus.kdn <- readRDS('validate_opt_res/auc.dat.fct.fil1.mus.kdn.pca.knn.rds')
auc.lis <- list(sc9=auc.arab.rt.sc9, sc51=auc.arab.rt.sc51, brain=auc.mus.brain, kdn=auc.mus.kdn)

# Cutoff to filter weak settings.
min.true <- 300; min.total <- 500; min.auc <- 0.5 

# Filter weak settings. 
spd.lis <- lapply(auc.lis, function(x) {
  df0 <- subset(x, true >= min.true & total >= min.total & auc >= min.auc)
  df0$spd.set 
}) 

auc.lis <- lapply(auc.lis, function(x) {
  x$auc[x$true < min.true | x$total < min.total | x$auc < min.auc] <- 0
  x  
})

tab <- table(unlist(spd.lis))
# Optimal spd.sets overlapping across validation results.  
spd.opt <- names(tab[tab==length(auc.lis)]); spd.opt 

# Rank optimal spd.opts by mean AUCs.
spd.auc.all <- do.call(rbind, auc.lis)[, c('spd.set', 'auc')]
spd.auc.ag <- aggregate(spd.auc.all$auc, list(spd.auc.all$spd.set), FUN=mean)
colnames(spd.auc.ag) <- c('spd.set', 'auc.mean') 
spd.auc.ag[order(-spd.auc.ag$auc.mean), ]


# Read results of random settings.
rdn.auc.arab.rt.sc9 <- readRDS('validate_opt_res/auc.dat.cpm.fil3.arab.rt.sc9.umap.rdn.rds')
rdn.auc.arab.rt.sc51 <- readRDS('validate_opt_res/auc.dat.cpm.fil3.arab.rt.sc51.umap.rdn.rds')
rdn.auc.mus.brain <- readRDS('validate_opt_res/auc.dat.cpm.fil3.mus.brain.umap.rdn.rds')
rdn.auc.mus.kdn <- readRDS('validate_opt_res/auc.dat.cpm.fil3.mus.kdn.umap.rdn.rds')

rdn.auc.lis <- list(sc9=rdn.auc.arab.rt.sc9, sc51=rdn.auc.arab.rt.sc51, brain=rdn.auc.mus.brain, kdn=rdn.auc.mus.kdn)

rdn.auc.lis <- lapply(rdn.auc.lis, function(x) {
  x$auc[x$true < min.true | x$total < min.total | x$auc < min.auc] <- 0
  x
})  


