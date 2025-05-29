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

##  brain significantly differentially expressed genes
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


# brain overlap differential exppression ----------------------------------------



inner_join(brain_eco12, 
           brain_eco18, 
           by = 'GeneID') 
## only 2 genes are differentially divergent between ecotypes
## across the two temperatures. 

inner_join(brain_plast_amb, 
           brain_plast_geo, 
           by = 'GeneID') %>% 
  select(GeneID) %>% 
  write_tsv('Brain_plasticity_amb_vs_geo.txt', 
            col_names = F)
## 47 genes are plastic in both the ambient and geothermal ecotypes

inner_join(brain_plast_amb, 
           brain_plast_hyb, 
           by = 'GeneID')%>% 
  select(GeneID) %>% 
  write_tsv('Brain_plasticity_amb_vs_hyb.txt', 
            col_names = F)
## 43 genes are plastic in both the amient and hybridl ecotypes

inner_join(brain_plast_geo, 
           brain_plast_hyb, 
           by = 'GeneID')%>% 
  select(GeneID) %>% 
  write_tsv('Brain_plasticity_geo_vs_hyb.txt', 
            col_names = F)
## 85 genes are plastic in both the geothermal and hybrid ecotypes

inner_join(brain_plast_amb, 
           brain_plast_geo, 
           by = 'GeneID') %>% 
  inner_join(.,
             brain_plast_hyb, 
             by = 'GeneID')%>% 
  select(GeneID) %>% 
  write_tsv('Brain_plasticity_amb_vs_geo_vs_hyb.txt', 
            col_names = F)

## 38 that are plastic between all three ecotypes

# brain quick volcano plots -----------------------------------------------------


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



# stickleback v5 annotation -----------------------------------------------

gene_annotation = read_tsv('~/Parsons_Postdoc/Stickleback_Genomic/Stickleback_Annotation_features/stickleback_v5_ensembl_genes.gff3.gz', 
                           col_names = F, 
                           skip = 1) %>% 
  # filter(X3 %in% c('gene', 
  #                  'exon', 
  #                  'CDS')) %>% 
  group_by(X1) %>% 
  arrange(X4, 
          X5) %>% 
  ## arrange each gene by its start and end points on each chromosome
  mutate(mid = X4 + (X5-X4)/2) %>% 
  dplyr::select(X1, 
                X3:X5, 
                X9:mid) %>% 
  rename(chromosome = X1, 
         feature = X3, 
         start = X4, 
         end = X5, 
         gene_id = X9, 
         position = mid) %>% 
  na.omit() %>% 
  separate(col = gene_id, 
                    into = c('ensemble_id', 
                             'gene_name', 
                             'parent_code', 
                             'gene_name2'), 
                    sep = ';') %>%
  separate(col = gene_name, 
           into = c('Garbage', 
                    'gene_name'), 
           sep = '=') %>% 
  dplyr::select(chromosome, 
                position, 
                start, ,
                end, 
                feature,
                gene_name) %>% 
  na.omit()

gene_annotation %>% 
  select(chromosome, 
         gene_id) %>% 
  pull(gene_id)


gene_annotation %>% 
  ungroup() %>% 
  select(gene_name) %>% 
  write_tsv('exonID.txt')
