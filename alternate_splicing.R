#### Alternative splicing analysis
## SKR hybrid data
## from Bethany's thesis

setwd('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/Bim_Bam/Brain_bam/')

library(BiocManager)

# BiocManager::install('GenomicFeatures')
# BiocManager::install('txdbmaker')

library(GenomicFeatures)
library(txdbmaker)
library(Rsamtools)
library(GenomicAlignments)
library(DEXSeq)

txdb = makeTxDbFromGFF('GCF_016920845.1_GAculeatus_UGA_version5_genomic.gtf.gz')

flattenedAnnotation = exonicParts(txdb, 
                                  linked.to.single.gene.only = T)
names(flattenedAnnotation) = sprintf('%s:E%0.3d', 
                                     flattenedAnnotation$gene_id, 
                                     flattenedAnnotation$exonic_part)

bam_files = dir('.', 'bam$')

bam_files = BamFileList(bam_files)

se = summarizeOverlaps(flattenedAnnotation, 
                       BamFileList(bam_files), 
                       singleEnd = F, 
                       fragments = F, 
                       ignore.strand = T)



colData(se)$ecotype = factor(c('SKRC','SKRC','SKRC','SKRC','SKRC','SKRC','SKRC','SKRC','SKRC','SKRC',
                               'SKRC','SKRC','SKRC','SKRC','SKRC','SKRC','SKRW','SKRW','SKRW','SKRW',
                               'SKRW','SKRW','SKRW','SKRW','SKRW','SKRW','SKRW','SKRW','SKRW','SKRW',
                               'SKRW','SKRW','SKRHYB','SKRHYB','SKRHYB','SKRHYB','SKRHYB','SKRHYB','SKRHYB',
                               'SKRHYB','SKRHYB','SKRHYB','SKRHYB','SKRHYB','SKRHYB','SKRHYB','SKRHYB'))

colData(se)$temp = factor(c('12','12','12','12','12','12','12','12','18','18',
                               '18','18','18','18','18','18','12','12','12','12',
                               '12','12','12','12','18','18','18','18','18','18',
                               '18','18','18','18','18','18','18','18','18',
                               '18','12','12','12','12','12','12','12'))

colData(se)$libtype = factor(c('paired-end','paired-end','paired-end','paired-end','paired-end','paired-end','paired-end','paired-end','paired-end','paired-end',
                            'paired-end','paired-end','paired-end','paired-end','paired-end','paired-end','paired-end','paired-end','paired-end','paired-end',
                            'paired-end','paired-end','paired-end','paired-end','paired-end','paired-end','paired-end','paired-end','paired-end','paired-end',
                            'paired-end','paired-end','paired-end','paired-end','paired-end','paired-end','paired-end','paired-end','paired-end',
                            'paired-end','paired-end','paired-end','paired-end','paired-end','paired-end','paired-end','paired-end'))


colData(se)$individual = factor(c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20',
                               '21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39',
                               '40','41','42','43','45','46','47','48'))

colData(se)$family = factor(c('006','008','073','073','073','082','082','082','008','008',
                            '008','073','073','073','082','082','025','025','025',
                            '029','034','034','034','046','025','025','034','034','034',
                            '046','046','046', '017','017','017','017','068','068','068',
                            '068','017','017','017','068','068','068','068'))



dxd = DEXSeqDataSetFromSE(se, design = ~ exon + ecotype:exon + temp:exon + ecotype:temp:exon)
