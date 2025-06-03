#### Alternative splicing analysis
## SKR hybrid data
## from Bethany's thesis

setwd('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/')

library(BiocManager)

# BiocManager::install('GenomicFeatures')
# BiocManager::install('txdbmaker')

library(GenomicFeatures)
library(txdbmaker)

txdb = makeTxDbFromGFF('GCF_016920845.1_GAculeatus_UGA_version5_genomic.gtf.gz')

flattenedAnnotation = exonicParts(txdb, 
                                  linked.to.single.gene.only = T)
names(flattenedAnnotation) = sprintf('%s:E%0.3d', 
                                     flattenedAnnotation$gene_id, 
                                     flattenedAnnotation$exonic_part)
