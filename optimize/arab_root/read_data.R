# This script illustrates how to read 1 bulk and 5 single cell data sets of Arabidopsis thaliana (Arabidopsis) root for coclustering optimization.   

# Read Arabidopsis root bulk data.
# RNASeq_counts.csv, ici_metadata.csv, cell_type_marker_bulk.txt are downloaded at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152766.

source('../../function/bulk_dat.R')
# blk.arab.rt
blk.all <- bulk_dat(rds='blk_data.rds', cnt.path='RNASeq_counts.csv', meta.path='ici_metadata.csv', marker.blk.path='cell_type_marker_bulk.txt')
blk.all[1:3, 1:5]

# Read Arabidopsis root single cell data.
# Download GSM4625995_sc_10_at_COPILOT.rds, GSM4625996_sc_11_COPILOT.rds, GSM4625997_sc_12_COPILOT.rds, GSM4626001_sc_30_COPILOT.rds, GSM4626002_sc_31_COPILOT.rds at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152766.

source('../../function/scell_dat.R'); source('../../function/df_match.R')

# Matching between bulk and cells.
df.match <- df_match()

# sc.arab.rt10
sc.all10 <- scell_dat('GSM4625995_sc_10_at_COPILOT.rds', df.match=df.match) # Read data.
as.matrix(sc.all10[1:3, 1:5])

# sc.arab.rt11
sc.all11 <- scell_dat('/GSM4625996_sc_11_COPILOT.rds', df.match=df.match) # Read data.
as.matrix(sc.all11[1:3, 1:5])

# sc.arab.rt12
sc.all12 <- scell_dat('/GSM4625997_sc_12_COPILOT.rds', df.match=df.match) # Read data.
as.matrix(sc.all12[1:3, 1:5])

# sc.arab.rt30
sc.all30 <- scell_dat('/GSM4626001_sc_30_COPILOT.rds', df.match=df.match) # Read data.
as.matrix(sc.all30[1:3, 1:5])

# sc.arab.rt31
sc.all31 <- scell_dat('/GSM4626002_sc_31_COPILOT.rds', df.match=df.match) # Read data.
as.matrix(sc.all31[1:3, 1:5])




