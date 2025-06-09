##############################
## gene network analysis
##
## Matt Brachmann (PhDMattyB)
##
## 09.06.2025
##
##############################


setwd('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/')

library(tidyverse)
library(WGCNA)

brain_exp = read_tsv('brain_gene_read_counts_table_all_final.tsv')
brain_limma = read_tsv('brain_limma_gene_list.txt')

brain_count_limma = inner_join(brain_exp, 
           brain_limma, 
           by = 'GeneID')


liver_exp = read_tsv('liver_gene_read_counts_table_all_final.tsv')
liver_limma = read_tsv("Liver_limma_gene_list.txt")

liver_count_limma = inner_join(liver_exp, 
                               liver_limma, 
                               by = 'GeneID')


metadata = names(brain_exp) %>% 
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
           sep = '-') %>% 
  unite(col = ecotemp, 
        c('ecotype',
          'temp'),
        sep = '_',
        remove = F)
