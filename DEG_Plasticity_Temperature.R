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
brain_limma_all = read_csv('Brain_LIMMA_model_results_all.csv')

## All genes 

brain_eco12_neutral = read_csv('Brain_eco_div_12.csv') %>% 
  filter(adj.P.Val > 0.05) %>% 
  mutate(status = 'Neutral')
brain_eco18_neutral = read_csv('Brain_eco_div_18.csv') %>% 
  filter(adj.P.Val > 0.05)%>% 
  mutate(status = 'Neutral')
brain_plast_amb_neutral = read_csv('Brain_ambient_plastic.csv') %>% 
  filter(adj.P.Val > 0.05)%>% 
  mutate(status = 'Neutral')
brain_plast_geo_neutral = read_csv('Brain_geothermal_plastic.csv') %>% 
  filter(adj.P.Val > 0.05)%>% 
  mutate(status = 'Neutral')
brain_plast_hyb_neutral = read_csv('Brain_hybrid_plastic.csv') %>% 
  filter(adj.P.Val > 0.05)%>% 
  mutate(status = 'Neutral')

## significantly differentially expressed genes
brain_eco12 = read_csv('Brain_eco_div_12_significant.csv') %>% 
  mutate(status = 'Outlier')
brain_eco18 = read_csv('Brain_eco_div_18_significant.csv')%>% 
  mutate(status = 'Outlier')
brain_plast_amb = read_csv('Brain_ambient_plastic_significant.csv')%>% 
  mutate(status = 'Outlier')
brain_plast_geo = read_csv('Brain_geothermal_plastic_significant.csv')%>% 
  mutate(status = 'Outlier')
brain_plast_hyb = read_csv('Brain_hybrid_plastic_significant.csv')%>% 
  mutate(status = 'Outlier')


brain_eco12_clean = bind_rows(brain_eco12, 
                              brain_eco12_neutral)

brain_eco18_clean = bind_rows(brain_eco18, 
                              brain_eco18_neutral)

brain_plast_amb_clean = bind_rows(brain_plast_amb, 
                              brain_plast_amb_neutral)

brain_plast_geo_clean = bind_rows(brain_plast_geo, 
                                  brain_plast_geo_neutral)

brain_plast_hyb_clean = bind_rows(brain_plast_hyb, 
                                  brain_plast_hyb_neutral)


ggplot(data = brain_eco12_clean, 
       aes(x = logFC, 
           y = adj.P.Val, 
           col = status))+
  geom_point()+
  scale_y_reverse()

ggplot(data = brain_eco18_clean, 
       aes(x = logFC, 
           y = adj.P.Val, 
           col = status))+
  geom_point()+
  scale_y_reverse()

ggplot(data = brain_plast_amb_clean, 
       aes(x = logFC, 
           y = adj.P.Val, 
           col = status))+
  geom_point()+
  scale_y_reverse()

ggplot(data = brain_plast_geo_clean, 
       aes(x = logFC, 
           y = adj.P.Val, 
           col = status))+
  geom_point()+
  scale_y_reverse()


ggplot(data = brain_plast_hyb_clean, 
       aes(x = logFC, 
           y = adj.P.Val, 
           col = status))+
  geom_point()+
  scale_y_reverse()
