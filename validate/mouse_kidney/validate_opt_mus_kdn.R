# Loading packages.
library(spatialHeatmap); library(BiocParallel); library(data.table)
source('function/bulk_dat.R')
source('function/scell_dat.R')
source('function/df_match.R')

# Obtain reproducible results.
set.seed(10)

# Data paths.
blk.pa.mus.kdn <- '~/bigdata/single_cell/data/mouse_kidney/result_PRJNA389326_PRJNA438336'
blk.pa.mus.glom <- '~/bigdata/single_cell/data/mouse_kidney/PRJNA435940/result'
sc.pa.mus.kdn <- '~/bigdata/single_cell/data/mouse_kidney'

# Matching table between bulk and single cells.
df.match.mus.kdn <- df_match_mus_kdn()

# Read data.
# Bulk tissues of mouse kidney.
blk.mus.kdn <- blk_dat_mus(pa=file.path(blk.pa.mus.kdn, 'countDFeByg.xls'))
blk.mus.kdn[1:3, ]
blk.mus.glom <- blk_dat_mus(pa=file.path(blk.pa.mus.glom, 'countDFeByg.xls'))
blk.mus.glom[1:3, ]

inter <- intersect(rownames(blk.mus.kdn), rownames(blk.mus.glom))
blk.mus.kdn <- cbind(blk.mus.glom[inter, ], blk.mus.kdn[inter, ])
blk.mus.kdn[1:3, ]

# Single cells of mouse kidney.
sc.mus.kdn <- sc_dat_mus_kdn(file.path(sc.pa.mus.kdn, 'GSE107585_Mouse_kidney_single_cell_datamatrix.txt'))
sc.mus.kdn[1:3, 1:5]

# Inital filtering before normalization.
blk.mus.kdn.init <- filter_data(data=blk.mus.kdn, pOA=c(0.05, 5), CV=c(0.05, 100)); dim(blk.mus.kdn.init)

# Reduce replicates, since too many reps/bulk impact coclustering.
blk.mus.kdn.init <- reduce_rep(dat=blk.mus.kdn.init, n=3)

mus.kdn.lis.init <- filter_cell(lis=list(sc.mus=sc.mus.kdn), bulk=blk.mus.kdn.init, gen.rm=NULL, min.cnt=1, p.in.cell=0.01, p.in.gen=0.05)

# Validate optimal settings.

# Optimal settings.
df.spd.opt <- read.table('df_spd_opt.txt', header=TRUE, row.name=1, sep='\t')
df.spd.opt[1:3, ]

# Normalization by sum.factor.
norm.fct.mus.kdn <- norm_multi(dat.lis=mus.kdn.lis.init, cpm=FALSE)
saveRDS(norm.fct.mus.kdn, file='validate_opt_res/norm.fct.mus.kdn.rds')
norm.fct.mus.kdn <- readRDS('validate_opt_res/norm.fct.mus.kdn.rds')

# Secondary filtering 1 (fil1). In the optimization, "fil1", "fil2", and "fil3" are equally optimal, so only "fil1" is used. 
blk.fct.fil1.mus.kdn <- filter_data(data=norm.fct.mus.kdn$bulk, pOA=c(0.1, 1), CV=c(0.1, 100)); dim(blk.fct.fil1.mus.kdn)
dat.fct.fil1.mus.kdn <- filter_cell(lis=list(sc.mus=norm.fct.mus.kdn$sc.mus), bulk=blk.fct.fil1.mus.kdn, gen.rm=NULL, min.cnt=1, p.in.cell=0.1, p.in.gen=0.01)

saveRDS(dat.fct.fil1.mus.kdn, file='validate_opt_res/dat.fct.fil1.mus.kdn.rds')
dat.fct.fil1.mus.kdn <- readRDS('validate_opt_res/dat.fct.fil1.mus.kdn.rds')

# Coclustering.
opt.mus.kdn <- cocluster(bulk=dat.fct.fil1.mus.kdn$bulk, cell=dat.fct.fil1.mus.kdn$sc.mus, df.match=df.match.mus.kdn, df.para=df.spd.opt, sc.dim.min=10, max.dim=50, sim=0.2, sim.p=0.8, dim=12, graph.meth='knn', dimred='PCA', sim.meth='spearman', return.all=FALSE, multi.core.par=MulticoreParam(workers=5), file='validate_opt_res/auc.dat.fct.fil1.mus.kdn.pca.knn')

opt.mus.kdn[1:3, ]
opt.mus.kdn <- readRDS('validate_opt_res/auc.dat.fct.fil1.mus.kdn.pca.knn.rds')

# Set auc = 0 for weak assignments.
opt.mus.kdn$auc[opt.mus.kdn$true < 300 | opt.mus.kdn$true < 500 | opt.mus.kdn$auc < 0.5] <- 0

# As comparison with optimal settings, random combinations of sub-optimal settings are generated and tested.

# Random settings. In the optimization, "fil4" is sub-optimal and should be used in the random settings. But it filters out too many genes and too less genes remain, so "fil3" is used instead.
par.rdn.mus.kdn <- random_para(fil.set=c('fil3'), norm='cpm', dimred='UMAP', graph.meth=c('knn', 'snn'), sim=round(seq(0.2, 0.8, by=0.1), 1), sim.p=round(seq(0.2, 0.8,by=0.1), 1), dim=seq(5, 40, by=1), df.spd.opt=df.spd.opt)
par.rdn.mus.kdn[1:3, ]

# Normalization by sum.factor + cpm.
norm.cpm.mus.kdn <- norm_multi(dat.lis=mus.kdn.lis.init, cpm=TRUE)
saveRDS(norm.cpm.mus.kdn, file='validate_opt_res/norm.cpm.mus.kdn.rds')
norm.cpm.mus.kdn <- readRDS('validate_opt_res/norm.cpm.mus.kdn.rds')

# Secondary filtering 3 (fil3). 
blk.cpm.fil3.mus.kdn <- filter_data(data=norm.cpm.mus.kdn$bulk, pOA=c(0.3, 1), CV=c(0.3, 100)); dim(blk.cpm.fil3.mus.kdn)
dat.cpm.fil3.mus.kdn <- filter_cell(lis=list(sc.mus=norm.cpm.mus.kdn$sc.mus), bulk=blk.cpm.fil3.mus.kdn, gen.rm=NULL, min.cnt=1, p.in.cell=0.3, p.in.gen=0.1)

saveRDS(dat.cpm.fil3.mus.kdn, file='validate_opt_res/dat.cpm.fil3.mus.kdn.rds')

# Coclustering.
rdn.mus.kdn <- cocluster(bulk=dat.cpm.fil3.mus.kdn$bulk, cell=dat.cpm.fil3.mus.kdn$sc.mus, df.match=df.match.mus.kdn, df.para=par.rdn.mus.kdn, sc.dim.min=10, max.dim=50, sim=0.2, sim.p=0.8, dim=12, graph.meth='knn', dimred='PCA', sim.meth='spearman', return.all=FALSE, multi.core.par=MulticoreParam(workers=5), file='validate_opt_res/auc.dat.cpm.fil3.mus.kdn.umap.rdn')

rdn.mus.kdn[1:3, ]
rdn.mus.kdn <- readRDS('validate_opt_res/auc.dat.cpm.fil3.mus.kdn.umap.rdn.rds')

# Set auc = 0 for weak assignments.
rdn.mus.kdn$auc[rdn.mus.kdn$true < 300 | rdn.mus.kdn$true < 500 | rdn.mus.kdn$auc < 0.5] <- 0

