df_match <- function() {
  SVGBulk <- c('NONHAIR', 'COLU', 'COLU', 'COLU', 'COLU', 'COLU', 'CORT', 'CORT', 'ENDO', 'LRC', 'LRC', 'LRC', 'LRC', 'STELE', 'STELE', 'STELE', 'STELE', 'STELE', 'STELE', 'STELE', 'QC', 'QC', 'QC', 'HAIR', 'STELE', 'STELE', 'STELE', 'STELE')
  
  cell <- c('atricho', 'colu.dist.colu', 'colu.dist.lrc', 'colu.proxi.colu', 'colu.proxi.lrc', 'colu', 'cortex', 'cortex.dist.lrc', 'endo', 'lat.rt.cap.dist.colu', 'lat.rt.cap.dist.lrc', 'lat.rt.cap.proxi.lrc', 'lat.rt.cap', 'metaphlo.comp.cell', 'metaxylem', 'phlo.po.per', 'procam', 'protophlo', 'protoxylem', 'protoxylem.dist.lrc', 'quies.cent', 'put.quies.cent', 'stem.niche', 'tricho', 'xylem.po.per', 'xylem', 'per', 'phlo')

  trueBulk <-c('NONHAIR,LRC_NONHAIR', 'COLU', 'COLU', 'COLU', 'COLU', 'COLU', 'CORT', 'CORT', 'ENDO,ENDO_QC', 'LRC_NONHAIR', 'LRC_NONHAIR', 'LRC_NONHAIR', 'LRC_NONHAIR', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'ENDO_QC,QC', 'ENDO_QC,QC', 'ENDO_QC,QC', 'HAIR', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE', 'XYLEM,PERI,PHLM_COMP,PHLOEM,STELE')

  df.match <- data.frame(SVGBulk=SVGBulk, cell=cell, trueBulk=trueBulk)
  return(df.match)
}
# Mouse brain bulk 533.
df_match_mus533 <- function() {
  SVGBulk <- c('cerebellum', 'hippocampus', 'cerebral.cortex', 'hippocampus', 'hypothalamus')
  cell <- c('cere', 'hipp', 'isocort', 'retrohipp', 'hypotha')
  trueBulk <- c('CERE', 'HIPP', 'CERE.CORTEX', 'HIPP', 'HYPOTHA')
  df.match <- data.frame(SVGBulk=SVGBulk, cell=cell, trueBulk=trueBulk)
  return(df.match)
}
# Mouse kidney.
df_match_mus_kdn <- function() {
  SVGBulk <- c('proximal.tubule', 'distal.convoluted.tubule', 'collecting.duct', 'loop.of.henle', 'loop.of.henle', 'collecting.duct', 'collecting.duct', 'glomeruli')
  cell <- c('proxi.tub', 'distal.con.tub', 'col.duct.prin.cell', 'loop.hen', 'endo', 'col.duct.inter.cell', 'col.duct.trans.cell', 'podo')
  trueBulk <- c('PTS1,PTS2,PTS3', 'DCT', 'CCD,OMCD,IMCD', 'DTL1,DTL2,DTL3,ATL,MTAL,CTAL', 'DTL1,DTL2,DTL3,ATL,MTAL,CTAL', 'CCD,OMCD,IMCD', 'CCD,OMCD,IMCD', 'GLOM')
  df.match <- data.frame(SVGBulk=SVGBulk, cell=cell, trueBulk=trueBulk)
  return(df.match)
}

# Mouse brain bulk 988
df_match_mus988 <- function() {
  SVGBulk <- c('hippocampus', 'hippocampus', 'prefrontal.cortex')
  cell <- c('hipp', 'retrohipp', 'isocort')
  trueBulk <- c('HIPP', 'HIPP', 'PRECORT')
  df.match <- data.frame(SVGBulk=SVGBulk, cell=cell, trueBulk=trueBulk)
  return(df.match)
}
