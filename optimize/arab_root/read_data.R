# This script illustrates how to read 1 bulk and 5 single cell data sets of Arabidopsis thaliana (Arabidopsis) root for coclustering optimization.   

# Read Arabidopsis root bulk data.
# RNASeq_counts.csv, ici_metadata.csv, cell_type_marker_bulk.txt are downloaded at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152766.

source('../../function/bulk_dat.R')
# RNASeq_counts.csv: bulk data. ici_metadata.csv, cell_type_marker_bulk.txt: annotation information of the bulk tissues. blk_data.rds: file to save the bulk data with annotation labels included.  
blk.arab.rt <- bulk_dat(cnt.path='RNASeq_counts.csv', meta.path='ici_metadata.csv', marker.blk.path='cell_type_marker_bulk.txt', rds='blk_data.rds')
blk.arab.rt[1:3, 1:5]

# Read Arabidopsis root single cell data.
# Download GSM4625995_sc_10_at_COPILOT.rds, GSM4625996_sc_11_COPILOT.rds, GSM4625997_sc_12_COPILOT.rds, GSM4626001_sc_30_COPILOT.rds, GSM4626002_sc_31_COPILOT.rds at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152766.

source('../../function/scell_dat.R'); source('../../function/df_match.R')

# Matching between bulk tissues and single cells.
df.match <- df_match()

# sc.arab.rt10
sc.arab.rt10 <- scell_dat('GSM4625995_sc_10_at_COPILOT.rds', df.match=df.match) 
as.matrix(sc.arab.rt10[1:3, 1:5])

# sc.arab.rt11
sc.arab.rt11 <- scell_dat('GSM4625996_sc_11_COPILOT.rds', df.match=df.match) # Read data.
as.matrix(sc.arab.rt11[1:3, 1:5])

# sc.arab.rt12
sc.arab.rt12 <- scell_dat('GSM4625997_sc_12_COPILOT.rds', df.match=df.match) # Read data.
as.matrix(sc.arab.rt12[1:3, 1:5])

# sc.arab.rt30
sc.arab.rt30 <- scell_dat('GSM4626001_sc_30_COPILOT.rds', df.match=df.match) # Read data.
as.matrix(sc.arab.rt30[1:3, 1:5])

# sc.arab.rt31
sc.arab.rt31 <- scell_dat('GSM4626002_sc_31_COPILOT.rds', df.match=df.match) # Read data.
as.matrix(sc.arab.rt31[1:3, 1:5])




