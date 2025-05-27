##############################
## DEG plasticity between temperatures
##
## Matt Brachmann (PhDMattyB)
##
## 27.05.2025
##
##############################

setwd('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/')

library(tidyverse)
library(edgeR)

theme_set(theme_bw())

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

## Limma and edger analyses were done in another script
## This is a follow up script from
## the SKR_hybrid_analysis script to look at temperature
## and ecotype effects more effectively

brain_common_genes = read_csv('Brain_DEG_common_edger_limma.csv')

liver_common_genes = read_csv('Liver_DEG_Common_edger_limma.csv')
