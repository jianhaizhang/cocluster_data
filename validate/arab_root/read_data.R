# This script illustrates how to read 1 bulk and 2 single cell data sets of Arabidopsis thaliana (Arabidopsis) root for validating optimal parameter settings.  

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
df.match <- df_match()
# sc.arab.rt9
sc.arab.rt9 <- scell_dat('GSM4625994_sc_9_at_COPILOT.rds', df.match=df.match) # Read data.                       
as.matrix(sc.arab.rt9[1:3, 1:5])
# sc.arab.rt51
sc.arab.rt51 <- scell_dat('GSM4626006_sc_51_COPILOT.rds', df.match=df.match) # Read data.                        
as.matrix(sc.arab.rt51[1:3, 1:5])



