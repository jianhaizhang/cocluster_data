# This script illustrates how to read 1 bulk and 1 single cell data sets of mouse brain for validating optimal parameter settings.  

# Mouse brain bulk data.
# The fastq files are downloaded at https://www.ncbi.nlm.nih.gov/bioproject/PRJNA725533 using sratoolkit (3.0.0) and trimmed using fastp (https://github.com/OpenGene/fastp). Then raw counts are obtained using systemPipeR (2.1.12).  

# Read mouse brain bulk data.
source('../../function/bulk_dat.R')
# blk.mus.brain
blk.mus.brain <- blk_dat_mus('bulk_mouse_brain.xls')
blk.mus.brain[1:3, ]

# Download GSE147747_expr_raw_counts_table.tsv, GSE147747_meta_table.tsv at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147747. 

# Read mouse brain single cell data.
source('../../function/scell_dat.R')
# sc.mus.brain
sc.mus.brain <- sc_dat_mus_brain(sc.pa= 'GSE147747_expr_raw_counts_table.tsv', meta.pa=   
'GSE147747_meta_table.tsv') 
sc.mus.brain[1:3, 1:5]

# Matching table.  
source('../../function/df_match.R')
match.mus.brain <- df_match_mus533()



