##############################
## Bethany Gene expression
## SKR HYBRIDS
##
## Matt Brachmann (PhDMattyB)
##
## 22.05.2025
##
##############################


setwd('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/')

library(tidyverse)
library(edgeR)

brain_exp = read_tsv('brain_gene_read_counts_table_all_final.tsv')
liver_exp = read_tsv('liver_gene_read_counts_table_all_final.tsv')

brain_metadata = names(brain_exp) %>% 
  as_tibble() %>% 
  slice(-1) %>% 
  separate(col = value, 
           into = c('ecotype', 
                    'temp', 
                    'family', 
                    'sample', 
                    'tissue'), 
           sep = '_') %>% 
  separate(col = ecotype, 
           into = c('sample_num', 
                    'ecotype'), 
           sep = '-')


brain_metadata %>% 
  group_by(ecotype, 
           temp) %>% 
  summarize(n = n())


## filter low expression edgeR

brain_dge_list = DGEList(brain_exp)
brain_norm = calcNormFactors(brain_dge_list)

mm = model.matrix(~0 + ecotype*temp, 
                  data = brain_metadata)

brain_keep = filterByExpr(brain_norm, 
                          min.count = 10,
                          mm)
sum(brain_keep) # number of genes retai
brain_keep = brain_norm[brain_keep,]

brain_voom = voom(brain_keep, mm, plot = T)
