# Loading packages.
library(spatialHeatmap); library(BiocParallel); library(data.table)
source('function/bulk_dat.R')
source('function/scell_dat.R')
source('function/df_match.R')

# Obtain reproducible results.
set.seed(10)

# Data paths.
sc.pa.mus<- '~/bigdata/single_cell/data/mouse_brain/single_cell'
blk.pa.mus533 <- '~/bigdata/single_cell/data/mouse_brain/PRJNA725533/result'

# Matching table between bulk and single cells.
df.match.mus.brain <- df_match_mus533()

# Read data.

# Bulk tissues of mouse brain.
blk.mus.brain <- blk_dat_mus(pa=file.path(blk.pa.mus533, 'countDFeByg.xls'))
blk.mus.brain[1:3, ]

# Single cells of mouse kidney.
sc.mus.brain <- sc_dat_mus_brain(sc.pa=file.path(sc.pa.mus, 'GSE147747_expr_raw_counts_table.tsv'), meta.pa=file.path(sc.pa.mus, 'GSE147747_meta_table.tsv'))
sc.mus.brain[1:3, 1:5]

# Inital filtering before normalization.
blk.mus.brain.init <- filter_data(data=blk.mus.brain, pOA=c(0.05, 5), CV=c(0.05, 100)); dim(blk.mus.brain)

mus.brain.lis.init <- filter_cell(lis=list(sc.mus=sc.mus.brain), bulk=blk.mus.brain.init, gen.rm=NULL, min.cnt=1, p.in.cell=0.01, p.in.gen=0.05)

# Validate optimal settings.

# Optimal settings.
df.spd.opt <- read.table('df_spd_opt.txt', header=TRUE, row.name=1, sep='\t')
df.spd.opt[1:3, ]

# Normalization by sum.factor.
norm.fct.mus.brain <- norm_multi(dat.lis=mus.brain.lis.init, cpm=FALSE)
saveRDS(norm.fct.mus.brain, file='validate_opt_res/norm.fct.mus.brain.rds')
norm.fct.mus.brain <- readRDS('validate_opt_res/norm.fct.mus.brain.rds')

# Secondary filtering 1 (fil1). In the optimization, "fil1", "fil2", and "fil3" are equally optimal, so only "fil1" is used.
blk.fct.fil1.mus.brain <- filter_data(data=norm.fct.mus.brain$bulk, pOA=c(0.1, 1), CV=c(0.1, 100)); dim(blk.fct.fil1.mus.brain)

dat.fct.fil1.mus.brain <- filter_cell(lis=list(sc.mus=norm.fct.mus.brain$sc.mus), bulk=blk.fct.fil1.mus.brain, gen.rm=NULL, min.cnt=1, p.in.cell=0.1, p.in.gen=0.01)

saveRDS(dat.fct.fil1.mus.brain, file='validate_opt_res/dat.fct.fil1.mus.brain.rds')
dat.fct.fil1.mus.brain <- readRDS('validate_opt_res/dat.fct.fil1.mus.brain.rds')

# Coclustering.
opt.mus.brain <- cocluster(bulk=dat.fct.fil1.mus.brain$bulk, cell=dat.fct.fil1.mus.brain$sc.mus, df.match=df.match.mus.brain, df.para=df.spd.opt, sc.dim.min=10, max.dim=50, sim=0.2, sim.p=0.8, dim=12, graph.meth='knn', dimred='PCA', sim.meth='spearman', return.all=FALSE, multi.core.par=MulticoreParam(workers=5), file='validate_opt_res/auc.dat.fct.fil1.mus.brain.pca.knn')

opt.mus.brain <- readRDS('validate_opt_res/auc.dat.fct.fil1.mus.brain.pca.knn.rds')

# Set auc = 0 for weak assignments.
opt.mus.brain$auc[opt.mus.brain$true < 300 | opt.mus.brain$true < 500 | opt.mus.brain$auc < 0.5] <- 0

# As comparison with optimal settings, random combinations of sub-optimal settings are generated and tested.

# Random settings. In the optimization, "fil4" is sub-optimal and should be used in the random settings. But it filters out too many genes and too less genes remain, so "fil3" is used instead.
par.rdn.mus.brain <- random_para(fil.set=c('fil3'), norm='cpm', dimred='UMAP', graph.meth=c('knn', 'snn'), sim=round(seq(0.2, 0.8, by=0.1), 1), sim.p=round(seq(0.2, 0.8,by=0.1), 1), dim=seq(5, 40, by=1), df.spd.opt=df.spd.opt)

par.rdn.mus.brain[1:3, ]

# Normalization by sum.factor + cpm.
norm.cpm.mus.brain <- norm_multi(dat.lis=mus.brain.lis.init, cpm=TRUE)
saveRDS(norm.cpm.mus.brain, file='validate_opt_res/norm.cpm.mus.brain.rds')
norm.cpm.mus.brain <- readRDS('validate_opt_res/norm.cpm.mus.brain.rds')

# Secondary filtering 3 (fil3).
blk.cpm.fil3.mus.brain <- filter_data(data=norm.cpm.mus.brain$bulk, pOA=c(0.3, 1), CV=c(0.3, 100)); dim(blk.cpm.fil3.mus.brain)

dat.cpm.fil3.mus.brain <- filter_cell(lis=list(sc.mus=norm.cpm.mus.brain$sc.mus), bulk=blk.cpm.fil3.mus.brain, gen.rm=NULL, min.cnt=1, p.in.cell=0.3, p.in.gen=0.1)

saveRDS(dat.cpm.fil3.mus.brain, file='validate_opt_res/dat.cpm.fil3.mus.brain.rds')

# Coclustering.
rdn.mus.brain <- cocluster(bulk=dat.cpm.fil3.mus.brain$bulk, cell=dat.cpm.fil3.mus.brain$sc.mus, df.match=df.match.mus.brain, df.para=par.rdn.mus.brain, sc.dim.min=10, max.dim=50, sim=0.2, sim.p=0.8, dim=12, graph.meth='knn', dimred='PCA', sim.meth='spearman', return.all=FALSE, multi.core.par=MulticoreParam(workers=5), file='validate_opt_res/auc.dat.cpm.fil3.mus.brain.umap.rdn')

rdn.mus.brain <- readRDS('validate_opt_res/auc.dat.cpm.fil3.mus.brain.umap.rdn.rds')

# Set auc = 0 for weak assignments.
rdn.mus.brain$auc[rdn.mus.brain$true < 300 | rdn.mus.brain$true < 500 | rdn.mus.brain$auc < 0.5] <- 0









