# Validate optimal parameter settings with 1 bulk and 2 single cell data sets of Arabidopsis thaliana (Arabidopsis) root.

# Read 1 bulk and 2 single cell data sets.   

# Read Arabidopsis root bulk data.
# RNASeq_counts.csv, ici_metadata.csv, cell_type_marker_bulk.txt are downloaded at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152766.

source('../../function/bulk_dat.R')
# blk.arab.rt
blk.arab.rt <- bulk_dat(rds='blk_data.rds', cnt.path='RNASeq_counts.csv', meta.path='ici_metadata.csv', marker.blk.path='cell_type_marker_bulk.txt')
blk.arab.rt[1:3, 1:5]

# Read Arabidopsis root single cell data.
# Download GSM4625994_sc_9_at_COPILOT.rds, GSM4626006_sc_51_COPILOT.rds at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152766.

source('../../function/scell_dat.R'); source('../../function/df_match.R')

# Matching between bulk and cells.
df.match.arab <- df_match()
# sc.arab.rt9
sc.arab.rt9 <- scell_dat('GSM4625994_sc_9_at_COPILOT.rds', df.match=df.match) # Read data.                       
as.matrix(sc.arab.rt9[1:3, 1:5])
# sc.arab.rt51
sc.arab.rt51 <- scell_dat('GSM4626006_sc_51_COPILOT.rds', df.match=df.match) # Read data.                        
as.matrix(sc.arab.rt51[1:3, 1:5])


# Loading packages.
library(spatialHeatmap); library(data.table); library(BiocParallel)
source('function/bulk_dat.R')
source('function/scell_dat.R')
source('function/df_match.R')

# Obtain reproducible results.
set.seed(10)

# Inital filtering before normalization.
blk.arab.rt.init <- filter_data(data=blk.arab.rt, pOA=c(0.05, 5), CV=c(0.05, 100)); dim(blk.arab.rt.init)

arab.rt.lis.init <- filter_cell(lis=list(sc9=sc.arab.rt9, sc51=sc.arab.rt51), bulk=blk.arab.rt.init, gen.rm='^ATCG|^ATCG', min.cnt=1, p.in.cell=0.01, p.in.gen=0.05)

# Validate optimal settings with sc9.

# Optimal settings.
df.spd.opt <- read.table('../df_spd_opt.txt', header=TRUE, row.name=1, sep='\t')

# Normalization by sum.factor.
norm.fct.arab.rt <- norm_multi(dat.lis=arab.rt.lis.init, cpm=FALSE)
saveRDS(norm.fct.arab.rt, file='../validate_opt_res/norm.fct.arab.rt.rds')
norm.fct.arab.rt <- readRDS('../validate_opt_res/norm.fct.arab.rt.rds')

# Secondary filtering 1 (fil1). In the optimization, "fil1", "fil2", and "fil3" are equally optimal, so only "fil1" is used.
blk.fct.fil1.arab.rt <- filter_data(data=norm.fct.arab.rt$bulk, pOA=c(0.1, 1), CV=c(0.1, 100)); dim(blk.fct.fil1.arab.rt)

dat.fct.fil1.arab.rt <- filter_cell(lis=norm.fct.arab.rt[c('sc9', 'sc51')], bulk=blk.fct.fil1.arab.rt, gen.rm='^ATCG|^ATCG', min.cnt=1, p.in.cell=0.1, p.in.gen=0.01)

saveRDS(dat.fct.fil1.arab.rt, file='../validate_opt_res/dat.fct.fil1.arab.rt.rds')
dat.fct.fil1.arab.rt <- readRDS('../validate_opt_res/dat.fct.fil1.arab.rt.rds')

# Coclustering: bulk + sc9.
opt.arab.rt.sc9 <- cocluster(bulk=dat.fct.fil1.arab.rt$bulk, cell=dat.fct.fil1.arab.rt$sc9, df.match=df.match.arab, df.para=df.spd.opt, sc.dim.min=10, max.dim=50, sim=0.2, sim.p=0.8, dim=12, graph.meth='knn', dimred='PCA', sim.meth='spearman', return.all=FALSE, multi.core.par=MulticoreParam(workers=5), file='../validate_opt_res/auc.dat.fct.fil1.arab.rt.sc9.pca.knn')

opt.arab.rt.sc9 <- readRDS('../validate_opt_res/auc.dat.fct.fil1.arab.rt.sc9.pca.knn.rds')

# Set auc = 0 for weak assignments.
opt.arab.rt.sc9$auc[opt.arab.rt.sc9$true < 300 | opt.arab.rt.sc9$true < 500 | opt.arab.rt.sc9$auc < 0.5 ] <- 0

# As comparison with optimal settings, random combinations of sub-optimal settings are generated and tested.

# Random settings. In the optimization, "fil4" is sub-optimal and should be used in the random settings. But it filters out too many genes and too less genes remain, so "fil3" is used instead. 
par.rdn.arab.rt.sc9 <- random_para(fil.set=c('fil3'), norm='cpm', dimred='UMAP', graph.meth=c('knn', 'snn'), sim=round(seq(0.2, 0.8, by=0.1), 1), sim.p=round(seq(0.2, 0.8,by=0.1), 1), dim=seq(5, 40, by=1), df.spd.opt=df.spd.opt)

par.rdn.arab.rt.sc9[1:3,  ]

# Normalization by sum.factor + cpm.
norm.cpm.arab.rt <- norm_multi(dat.lis=arab.rt.lis.init, cpm=TRUE)
saveRDS(norm.cpm.arab.rt, file='../validate_opt_res/norm.cpm.arab.rt.rds')
norm.cpm.arab.rt <- readRDS('../validate_opt_res/norm.cpm.arab.rt.rds')

# Secondary filtering 3 (fil3).
blk.cpm.fil3.arab.rt <- filter_data(data=norm.cpm.arab.rt$bulk, pOA=c(0.3, 1), CV=c(0.3, 100)); dim(blk.cpm.fil3.arab.rt)

dat.cpm.fil3.arab.rt <- filter_cell(lis=norm.cpm.arab.rt[c('sc9', 'sc51')], bulk=blk.cpm.fil3.arab.rt, gen.rm='^ATCG|^ATCG', min.cnt=1, p.in.cell=0.3, p.in.gen=0.1)

saveRDS(dat.cpm.fil3.arab.rt, file='../validate_opt_res/dat.cpm.fil3.arab.rt.rds')
dat.cpm.fil3.arab.rt <- readRDS('../validate_opt_res/dat.cpm.fil3.arab.rt.rds')

# Coclustering bulk + sc9.
rdn.arab.rt.sc9 <- cocluster(bulk=dat.cpm.fil3.arab.rt$bulk, cell=dat.cpm.fil3.arab.rt$sc9, df.match=df.match.arab, df.para=par.rdn.arab.rt.sc9, sc.dim.min=10, max.dim=50, sim=0.2, sim.p=0.8, dim=12, graph.meth='knn', dimred='PCA', sim.meth='spearman', return.all=FALSE, multi.core.par=MulticoreParam(workers=5), file='../validate_opt_res/auc.dat.cpm.fil3.arab.rt.sc9.umap.rdn')

rdn.arab.rt.sc9 <- readRDS('../validate_opt_res/auc.dat.cpm.fil3.arab.rt.sc9.umap.rdn.rds')

# Set auc = 0 for weak assignments.
rdn.arab.rt.sc9$auc[rdn.arab.rt.sc9$true < 300 | rdn.arab.rt.sc9$true < 500 | rdn.arab.rt.sc9$auc < 0.5] <- 0


# Validate optimal settings with sc51.

# Coclustering: bulk + sc51
opt.arab.rt.sc51 <- cocluster(bulk=dat.fct.fil1.arab.rt$bulk, cell=dat.fct.fil1.arab.rt$sc51, df.match=df.match.arab, df.para=df.spd.opt, sc.dim.min=10, max.dim=50, sim=0.2, sim.p=0.8, dim=12, graph.meth='knn', dimred='PCA', sim.meth='spearman', return.all=FALSE, multi.core.par=MulticoreParam(workers=5), file='../validate_opt_res/auc.dat.fct.fil1.arab.rt.sc51.pca.knn')

opt.arab.rt.sc51[1:3, ]
opt.arab.rt.sc51 <- readRDS('../validate_opt_res/auc.dat.fct.fil1.arab.rt.sc51.pca.knn.rds')

# Set auc = 0 for weak assignments.
opt.arab.rt.sc51$auc[opt.arab.rt.sc51$true < 300 | opt.arab.rt.sc51$true < 500 | opt.arab.rt.sc51$auc < 0.5 ] <- 0

# Sub-optimal/random settings. In the optimization, "fil4" is sub-optimal and should be used in the random settings. But it filters out too many genes and too less genes remain, so "fil3" is used instead. 
par.rdn.arab.rt.sc51 <- random_para(fil.set=c('fil3'), norm='cpm', dimred='UMAP', graph.meth=c('knn', 'snn'), sim=round(seq(0.2, 0.8, by=0.1), 1), sim.p=round(seq(0.2, 0.8,by=0.1), 1), dim=seq(5, 40, by=1), df.spd.opt=df.spd.opt)

par.rdn.arab.rt.sc51[1:3,  ]

# Coclustering bulk + sc51.
rdn.arab.rt.sc51 <- cocluster(bulk=dat.cpm.fil3.arab.rt$bulk, cell=dat.cpm.fil3.arab.rt$sc51, df.match=df.match.arab, df.para=par.rdn.arab.rt.sc51, sc.dim.min=10, max.dim=50, sim=0.2, sim.p=0.8, dim=12, graph.meth='knn', dimred='PCA', sim.meth='spearman', return.all=FALSE, multi.core.par=MulticoreParam(workers=5), file='../validate_opt_res/auc.dat.cpm.fil3.arab.rt.sc51.umap.rdn')

rdn.arab.rt.sc51[1:3, ]
rdn.arab.rt.sc51 <- readRDS('../validate_opt_res/auc.dat.cpm.fil3.arab.rt.sc51.umap.rdn.rds')

# Set auc = 0 for weak assignments.
rdn.arab.rt.sc51$auc[rdn.arab.rt.sc51$true < 300 | rdn.arab.rt.sc51$true < 500 | rdn.arab.rt.sc51$auc < 0.5] <- 0









