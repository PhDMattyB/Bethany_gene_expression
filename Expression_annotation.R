##############################
## Genes of interest
##
## Matt Brachmann (PhDMattyB)
##
## 29.01.2025
##
##############################

library(tidyverse)
library(data.table)

setwd('~/Parsons_Postdoc/Stickleback_Genomic/Stickleback_Annotation_features/')

gene_annotation = read_tsv('stickleback_v5_ensembl_genes.gff3.gz', 
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
  na.omit()



gene_annotation %>% 
  arrange(chromosome) %>%
  filter(feature == 'gene')
  group_by(feature) %>% 
  distinct(feature)


goi = read_csv('~/Parsons_Postdoc/Bethany_gene_expression/GLMER_gene_expression_ecotemp_pval0.01.csv')

goi_names = goi %>% 
  select(gene_name, 
         mean_expression_relative, 
         ecow_temp18_pval) %>% 
  rename(ensemble_name = gene_name)

gene_metadata = gene_annotation %>% 
  filter(feature == 'gene') %>% 
  select(position,
         chromosome,
         feature)

ensemlbe_annotation_data = gene_annotation %>% 
  filter(feature == 'gene') %>% 
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
  select(-trash) %>% 
  separate(ensemble_id, 
           into = c('ensemble_name', 
                    'trash'), 
           sep = '.CDS') %>% 
  select(-trash) 

annotation_data = bind_cols(gene_metadata, 
          ensemble_annotation_genes)


cleaned_data = inner_join(goi_names, 
           annotation_data, 
           by = 'ensemble_name') %>% 
  select(ensemble_name, 
         gene_name, 
         chromosome, 
         position, 
         mean_expression_relative, 
         ecow_temp18_pval,
         feature, 
         relationship) %>% 
  separate(ensemble_name, 
           into = c('ensemble_name', 
                    'trash'), 
           sep = '_') %>% 
  select(-trash)


cleaned_data %>% 
  select(gene_name) %>% 
  # write_csv('Expression_gene_names_ecotemp_pval0.01.csv')
  write_tsv('Expression_gene_names_ecotemp_pval0.01.txt')
