## reanalysing Bethanys hybrid data from her second data chapter of her PhD
## the data from the last data chapter needs to be fixed apparently. 
## using a GLMM to model gene expression differences


# Big package -------------------------------------------------------------

library(tidyverse)
library(lme4)
library(optimx)
library(sjPlot)



# Read the data -----------------------------------------------------------

setwd("/media/shaun/2d49aafa-914e-40e6-a1ea-142b85301ef9/rna_aligned/F1 hybrids")

brain_meta = read_csv('hyb_data.csv')
brain_data = read_tsv('brain_gene_read_counts_table_all_final.tsv')

names(brain_data)

liver_meta = read_csv('hyb_data.csv')
liver_data = read_tsv('liver_gene_read_counts_table_all_final.tsv')


# Brain variant filter ----------------------------------------------------

brain_data_rot = brain_data %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID') %>%
  t() %>%
  as_tibble()

brain_full = bind_cols(brain_meta, 
                      brain_data_rot)

brain_var = brain_full %>%
  select(pop, 
         type,
         temp,
         starts_with('ENSGAC')) %>% 
  group_by(pop, 
           type, 
           temp) %>%
  summarise(across(everything(), var, na.rm = F))%>%
  ungroup() %>% 
  select(starts_with('ENSGAC')) %>% 
  colSums(abs(.)>0) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  rename(variance = 2) %>% 
  filter(variance > 0.000000e+00)

brain_orig_data = brain_data %>%
  rownames_to_column() %>% 
  as_tibble()

brain_cleaned = left_join(brain_var, 
                        brain_orig_data, 
                        by = 'rowname') %>% 
  select(-variance)

brain_cleaned = brain_cleaned %>% 
  as.data.frame() %>% 
  remove_rownames() %>% 
  column_to_rownames(var = 'rowname')


# Brain GLMER -------------------------------------------------------------


# liver variant filter ----------------------------------------------------

liver_data_rot = liver_data %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID') %>%
  t() %>%
  as_tibble()

liver_full = bind_cols(liver_meta, 
                       liver_data_rot)

liver_var = liver_full %>%
  select(pop, 
         type,
         temp,
         starts_with('ENSGAC')) %>% 
  group_by(pop, 
           type, 
           temp) %>%
  summarise(across(everything(), var, na.rm = F))%>%
  ungroup() %>% 
  select(starts_with('ENSGAC')) %>% 
  colSums(abs(.)>0) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  rename(variance = 2) %>% 
  filter(variance > 0.000000e+00)

liver_orig_data = liver_data %>%
  rownames_to_column() %>% 
  as_tibble()

liver_cleaned = left_join(liver_var, 
                          liver_orig_data, 
                          by = 'rowname') %>% 
  select(-variance)

liver_cleaned = liver_cleaned %>% 
  as.data.frame() %>% 
  remove_rownames() %>% 
  column_to_rownames(var = 'rowname')

# liver GLMER -------------------------------------------------------------


