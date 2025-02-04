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

brain_meta = read_csv('hyb_data.csv') %>% 
  separate(col = fam, 
           into = c('trash', 
                    'fam'), 
           sep = '_') %>% 
  select(-trash)
brain_data = read_tsv('brain_gene_read_counts_table_all_final.tsv')


liver_meta = read_csv('hyb_data.csv') %>% 
  separate(col = fam, 
           into = c('trash', 
                    'fam'), 
           sep = '_') %>% 
  select(-trash)
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
  as_tibble() %>% 
  select(-rowname) %>% 
  rename(rowname = GeneID)

liver_cleaned = left_join(liver_var, 
                          liver_orig_data, 
                          by = 'rowname') %>% 
  select(-variance)

liver_cleaned = liver_cleaned %>% 
  # select(-rowname)
  as.data.frame() %>% 
  remove_rownames() %>% 
  column_to_rownames(var = 'rowname')

# liver GLMER -------------------------------------------------------------
liver_cleaned2 <- data.frame(t(liver_cleaned))
liver_meta = as.data.frame(liver_meta)

liver_lme4_results <- as.data.frame( matrix( nrow =14797,  ncol = 10) )

for(i in 1:ncol(liver_cleaned2)){
  focal_gene = liver_cleaned2[,i]
  all_others = rowSums(liver_cleaned2[,-i])
  Y = cbind(focal_gene, 
            all_others)
  liver_Model = glmer(Y ~ type * temp * sex + (1|fam),
                family = "binomial",
                glmerControl(optimizer="bobyqa", 
                             optCtrl = list(maxfun = 100000)),
                data = liver_meta)
  liver_glmer_results = summary(liver_Model)
  # mod_coef = coef(summary(Model))
  ## dont need the second one. It does the same as summary
  # gene_results2 = anova(Model,
  #                       test = 'LRT')
  liver_lme4_results[i, 1] = colnames(liver_cleaned2)[i]
  liver_lme4_results[i, 2] = sum(liver_cleaned2[,i]/sum(liver_cleaned2))
  liver_lme4_results[i, 3] = liver_glmer_results$coefficients[1,4]
  liver_lme4_results[i, 4] = liver_glmer_results$coefficients[2,4]
  liver_lme4_results[i, 5] = liver_glmer_results$coefficients[3,4]
  liver_lme4_results[i, 6] = liver_glmer_results$coefficients[4,4]
  liver_lme4_results[i, 7] = liver_glmer_results$coefficients[5,4]
  liver_lme4_results[i, 8] = liver_glmer_results$coefficients[6,4]
  # liver_lme4_results[i, 9] = liver_glmer_results$coefficients[7,4]
  # liver_lme4_results[i, 10] = liver_glmer_results$coefficients[8,4]
  # liver_lme4_results[i, 11] = liver_glmer_results$coefficients[9,4]
  # liver_lme4_results[i, 12] = liver_glmer_results$coefficients[10,4]
  # lme4_results_table[i, 13] = gene_results2$coefficients[11,4]
}

liver_glmer_results

# head(results.table)
# tail(results.table)                   
names(liver_lme4_results) <- c('gene_name', 
                               'mean_expression_relative',  
                               'popMYV_pval', 
                               'popSKR_pval',  
                               'ecoW_pval', 
                               'temp18_pval',  
                               'popMYV_ecoW_pval', 
                               'popSKR_ecoW_pval',  
                               'popMYV_temp18_pval', 
                               'popSKR_temp18_pval',  
                               'ecow_temp18_pval', 
                               # 'popMYV_ecow_temp18_pval',  
                               'popSKR_ecow_temp18_pval')


liver_lme4_results[,1:12] %>% 
  as_tibble() %>%
  filter(ecow_temp18_pval < 0.01) %>% 
  write_csv('GLMER_LIVER_gene_expression_ecotemp_pval0.01.csv')
# write_csv('GLMER_gene_expression_fam_rand_effect.csv')


## The glmm with nested effects is the same as this one
