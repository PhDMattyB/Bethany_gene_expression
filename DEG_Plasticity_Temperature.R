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


liver_common_genes = read_csv('Liver_DEG_Common_edger_limma.csv')


# Brain temp plast --------------------------------------------------------
brain_common_genes = read_csv('Brain_DEG_common_edger_limma.csv')
# brain_exp = read_tsv('brain_gene_read_counts_table_all_final.tsv')
brain_limma_all = read_csv('Brain_LIMMA_model_results_all.csv')

brain_significant = brain_limma_all %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')

brain_neutral = brain_limma_all %>% 
  filter(adj.P.Val > 0.05) %>% 
  mutate(status = 'Neutral')

brain_limma_all = bind_rows(brain_significant, 
                            brain_neutral)

ggplot()+
  geom_point(data = brain_neutral, 
             aes(x = eco18, 
                 y = eco12, 
                 col = status))+
  geom_point(data = brain_significant, 
             aes(x = eco18, 
                 y = eco12, 
                 col = status))


ggplot(data = brain_limma_all, 
       aes(x = plast_amb, 
           y = adj.P.Val, 
           col = status))+
  geom_point()+
  scale_y_reverse()

ggplot(data = brain_limma_all, 
       aes(x = plast_geo, 
           y = adj.P.Val, 
           col = status))+
  geom_point()+
  scale_y_reverse()

ggplot(data = brain_limma_all, 
       aes(x = plast_hyb, 
           y = adj.P.Val, 
           col = status))+
  geom_point()+
  scale_y_reverse()
