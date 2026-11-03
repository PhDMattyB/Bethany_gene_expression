library(tidyverse)
library(data.table)
library(qvalue)

setwd('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/')

## extract all of the stickle genome annotation data
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
  dplyr::rename(chromosome = X1, 
         feature = X3, 
         start = X4, 
         end = X5, 
         gene_id = X9, 
         position = mid) %>% 
  na.omit()


gene_metadata = gene_annotation %>% 
  # filter(feature == 'gene') %>% 
  dplyr::select(position,
         chromosome,
         feature)

ensemlbe_annotation_data = gene_annotation %>% 
  # filter(feature == 'gene') %>% 
  pull(gene_id) %>% 
  as_tibble() %>% 
  separate(value, 
           into = c('ensemble_id', 
                    'gene_name',
                    'relationship'), 
           sep = ';')

ensemble_annotation_genes = ensemlbe_annotation_data %>% 
  separate(ensemble_id, 
           into = c('trash', 
                    'ensemble_id'), 
           sep = '=') %>%
  separate(gene_name, 
           into = c('trash', 
                    'gene_name'), 
           sep = '=') %>% 
  dplyr::select(-trash) %>% 
  separate(ensemble_id, 
           into = c('ensemble_name', 
                    'trash'), 
           sep = '.CDS') %>% 
  dplyr::select(-trash) 

annotation_data = bind_cols(gene_metadata, 
                            ensemble_annotation_genes)


# Eco div 12 anno ---------------------------------------------------------
read_csv('Brain_eco_div_12.csv')%>% 
  filter(adj.P.Val <= 0.05) %>% 
  dplyr::rename(ensemble_name  = GeneID) %>% 
  inner_join(., 
           annotation_data, 
           by = 'ensemble_name') %>% 
  dplyr::select(gene_name, 
                logFC) %>% 
  write_tsv('Brain_eco_div_12_genes.txt')



# Eco div 18 anno ---------------------------------------------------------



read_csv('Brain_eco_div_18.csv')%>% 
  filter(adj.P.Val <= 0.05) %>% 
  dplyr::rename(ensemble_name  = GeneID) %>% 
  inner_join(., 
             annotation_data, 
             by = 'ensemble_name') %>% 
  dplyr::select(gene_name, 
                logFC) %>% 
  write_tsv('Brain_eco_div_18_genes.txt')



# Div amb hyb 12 ----------------------------------------------------------
read_csv('Brain_amb_hyb_12_div.csv')%>% 
  filter(adj.P.Val <= 0.05) %>% 
  dplyr::rename(ensemble_name  = GeneID) %>% 
  inner_join(., 
             annotation_data, 
             by = 'ensemble_name') %>% 
  dplyr::select(gene_name, 
                logFC) %>% 
  write_tsv('Brain_amb_hyb_12_div_genes.txt')

# Div amb hyb 18 ----------------------------------------------------------
read_csv('Brain_amb_hyb_18_div.csv')%>% 
  filter(adj.P.Val <= 0.05) %>% 
  dplyr::rename(ensemble_name  = GeneID) %>% 
  inner_join(., 
             annotation_data, 
             by = 'ensemble_name') %>% 
  dplyr::select(gene_name, 
                logFC) %>% 
  write_tsv('Brain_amb_hyb_18_div_genes.txt')


# Div geo hyb 12 ----------------------------------------------------------

read_csv('Brain_geo_hyb_12_div.csv')%>% 
  filter(adj.P.Val <= 0.05) %>% 
  dplyr::rename(ensemble_name  = GeneID) %>% 
  inner_join(., 
             annotation_data, 
             by = 'ensemble_name') %>% 
  dplyr::select(gene_name, 
                logFC) %>% 
  write_tsv('Brain_geo_hyb_12_div_genes.txt')


# Div geo hyb 18 ----------------------------------------------------------

read_csv('Brain_geo_hyb_18_div.csv')%>% 
  filter(adj.P.Val <= 0.05) %>% 
  dplyr::rename(ensemble_name  = GeneID) %>% 
  inner_join(., 
             annotation_data, 
             by = 'ensemble_name') %>% 
  dplyr::select(gene_name, 
                logFC) %>% 
  write_tsv('Brain_geo_hyb_18_div_genes.txt')



# trans div amb hyb 12 ----------------------------------------------------
read_csv('')%>% 
  filter(adj.P.Val <= 0.05) %>% 
  dplyr::rename(ensemble_name  = GeneID) %>% 
  inner_join(., 
             annotation_data, 
             by = 'ensemble_name') %>% 
  dplyr::select(gene_name, 
                logFC) %>% 
  write_tsv('Brain_geo_hyb_18_div_genes.txt')


# trans div amb hyb 18 ----------------------------------------------------


# trans div geo hyb 12 ----------------------------------------------------


# trans div geo hyb 18 ----------------------------------------------------



