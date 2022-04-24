# Read validation results.
auc.arab.rt.sc9 <- readRDS('validate_opt_res/auc.dat.fct.fil1.arab.rt.sc9.pca.knn.rds')
auc.arab.rt.sc51 <- readRDS('validate_opt_res/auc.dat.fct.fil1.arab.rt.sc51.pca.knn.rds')
auc.mus.brain <- readRDS('validate_opt_res/auc.dat.fct.fil1.mus.brain.pca.knn.rds')
auc.mus.kdn <- readRDS('validate_opt_res/auc.dat.fct.fil1.mus.kdn.pca.knn.rds')

auc.lis <- list(auc.arab.rt.sc9, auc.arab.rt.sc51, auc.mus.brain, auc.mus.kdn)
# Cutoff to filter weak settings.
min.true <- 300; min.total <- 500; min.auc <- 0.5
# Filter weak settings.
spd.lis <- lapply(auc.lis, function(x) { 
  df0 <- subset(x, true >= min.true & total >= min.total & auc >= min.auc) 
  df0$spd.set
})

tab <- table(unlist(spd.lis))
# Optimal spd.sets overlapping across validation results.
spd.opt <- names(tab[tab==length(auc.lis)])


