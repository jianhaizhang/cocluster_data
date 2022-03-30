# This script illustrates how to read 2 bulk and 1 single cell data sets of mouse kidney for validating optimal parameter settings.  

# Mouse kidney bulk data.
# The fastq file of PTS2 are downloaded at https://www.ncbi.nlm.nih.gov/bioproject/PRJNA438336. The fastq files of CCDs, cTALs are downloaded at https://www.ncbi.nlm.nih.gov/bioproject/PRJNA389326. The fastq file of glomeruli is downloaded at https://www.ncbi.nlm.nih.gov/bioproject/PRJNA435940. All fastq files are downloaded using sratoolkit (3.0.0) and trimmed using fastp (https://github.com/OpenGene/fastp). Then raw counts are obtained using systemPipeR (2.1.12).  

# Read mouse kidney bulk data.
source('../../function/bulk_dat.R')
# PTS2, CCDs, cTALs.
blk.all.mus.pct <- blk_dat_mus(pa='bulk_mouse_kidney_p.c.t.xls') 
blk.all.mus.pct[1:3, ]
# Glomeruli.
blk.all.mus.glom <- blk_dat_mus(pa='bulk_mouse_kidney_glom.xls')
blk.all.mus.glom[1:3, ] 
inter <- intersect(rownames(blk.all.mus.pct), rownames(blk.all.mus.glom)) 
# blk.mus.kdn
blk.all.mus.kdn <- cbind(blk.all.mus.glom[inter, ], blk.all.mus.pct[inter, ]) 
blk.all.mus.kdn[1:3, ] 


# Download GSE107585_Mouse_kidney_single_cell_datamatrix.txt at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107585. 

# Read mouse kidney single cell data.
# sc.mus.kdn
sc.mus.kdn <- sc_dat_mus_kdn('GSE107585_Mouse_kidney_single_cell_datamatrix.txt') 
sc.mus.kdn[1:3, ]


