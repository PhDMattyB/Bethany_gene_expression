
library(tidyverse)
library(lme4)
library(optimx)
library(sjPlot)

setwd("/media/shaun/2d49aafa-914e-40e6-a1ea-142b85301ef9/rna_aligned/f1_labfish/big_batch_2/tsv files/brain/data_analysis/f1_DATA_ANALYSIS")
annot_table_brain <- read.csv("brain_background_list_2.csv")
brainrawcounts  <- read.delim("F1_lab_brain_over10Mreads_alltsv3.csv")
brainsampledata <- read.csv("F1_labfish_sampledata_brain_over10M(1).csv")
braincounttable <- read.csv("F1_lab_brain_over10Mreads_alltsv3.csv", row.names=1)



# data dive ---------------------------------------------------------------


gene_data = df2 %>% 
  as_tibble()
meta_data = meta %>% 
  as_tibble()

full_data = bind_cols(meta_data, 
          gene_data)

gene_var = full_data %>%
  select(pop, 
         eco,
         temp,
         starts_with('ENSGAC')) %>% 
  group_by(pop, 
           eco, 
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

orig_count_data = braincounttable %>%
  rownames_to_column() %>% 
  as_tibble()

cleand_data = left_join(gene_var, 
          orig_count_data, 
          by = 'rowname') %>% 
  select(-variance)

df = cleand_data %>% 
  as.data.frame() %>% 
  remove_rownames() %>% 
  column_to_rownames(var = 'rowname')

# glm model ---------------------------------------------------------------

df2 <- data.frame(t(df))
meta = brainsampledata

for(i in 1:ncol(df2)){
  focal_gene = df2[,i]
  all_others = rowSums(df2[,-i])
  Y = cbind(focal_gene, 
            all_others)
  Model = glm(Y ~ pop*eco*temp,
              family = quasibinomial,
              data = meta)
  # Model = glmer(Y ~ pop * eco * temp + (1|fam), 
  #               family = binomial(),
  #               data = meta)
  gene_results = summary(Model)
  # mod_coef = coef(summary(Model))
## dont need the second one. It does the same as summary
    # gene_results2 = anova(Model,
  #                       test = 'LRT')
  results.table[i, 1] = colnames(df2)[i]
  results.table[i, 2] = sum(df2[,i]/sum(df2))
  results.table[i, 3] = gene_results$coef[1,4]
  results.table[i, 4] = gene_results$coef[2,4]
  results.table[i, 5] = gene_results$coef[3,4]
  results.table[i, 6] = gene_results$coef[4,4]
  results.table[i, 7] = gene_results$coef[5,4]
  results.table[i, 8] = gene_results$coef[6,4]
  results.table[i, 9] = gene_results$coef[7,4]
  results.table[i, 10] = gene_results$coef[8,4]
  results.table[i, 11] = gene_results$coef[9,4]
  results.table[i, 12] = gene_results$coef[10,4]
  results.table[i, 13] = gene_results$coef[11,4]
}

gene_results

# head(results.table)
# tail(results.table)                   
names(results.table) <- c('gene_name', 
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
                          'popMYV_ecow_temp18_pval',  
                          'popSKR_ecow_temp18_pval')


## in this model only ECOW, popMYV*ECOW, and popSKR*ECOW is significant
## No temperature related treatment effects are significant
## have not accounted for family effects
results.table[,1:13] %>% 
  as_tibble() %>% 
  write_csv('GLM_gene_expression_Dans_code_No_random_effects.csv')
  arrange() %>% 
  filter(pvalue < 0.01)


# lme4 random fam effects -------------------------------------------------


df2 <- data.frame(t(df))

meta = brainsampledata
lme4_results_table <- as.data.frame( matrix( nrow =15694,  ncol = 15) )

for(i in 1:ncol(df2)){
  focal_gene = df2[,i]
  all_others = rowSums(df2[,-i])
  Y = cbind(focal_gene, 
            all_others)
  Model = glmer(Y ~ pop * eco * temp + (1|fam),
                family = binomial(),
                # REML = T,
                control = glmerControl(
                  optimizer ='optimx', optCtrl=list(method='nlminb')),
                data = meta)
  # Model = glmer(Y ~ eco*temp + (1|fam),
  #               family = binomial(),
  #               # REML = T,
  #               control = glmerControl(
  #                 optimizer ='optimx', optCtrl=list(method='nlminb')),
  #               data = meta)
  gene_results2 = summary(Model)
  # mod_coef = coef(summary(Model))
  ## dont need the second one. It does the same as summary
  # gene_results2 = anova(Model,
  #                       test = 'LRT')
  lme4_results_table[i, 1] = colnames(df2)[i]
  lme4_results_table[i, 2] = sum(df2[,i]/sum(df2))
  lme4_results_table[i, 3] = gene_results2$coefficients[1,4]
  lme4_results_table[i, 4] = gene_results2$coefficients[2,4]
  lme4_results_table[i, 5] = gene_results2$coefficients[3,4]
  lme4_results_table[i, 6] = gene_results2$coefficients[4,4]
  lme4_results_table[i, 7] = gene_results2$coefficients[5,4]
  lme4_results_table[i, 8] = gene_results2$coefficients[6,4]
  lme4_results_table[i, 9] = gene_results2$coefficients[7,4]
  lme4_results_table[i, 10] = gene_results2$coefficients[8,4]
  lme4_results_table[i, 11] = gene_results2$coefficients[9,4]
  lme4_results_table[i, 12] = gene_results2$coefficients[10,4]
  # lme4_results_table[i, 13] = gene_results2$coefficients[11,4]
}

gene_results2

# head(results.table)
# tail(results.table)                   
names(lme4_results_table) <- c('gene_name', 
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


lme4_results_table[,1:12] %>% 
  as_tibble() %>%
  filter(ecow_temp18_pval < 0.01) %>% 
  write_csv('GLMER_gene_expression_ecotemp_pval0.01.csv')
  # write_csv('GLMER_gene_expression_fam_rand_effect.csv')


## The glmm with nested effects is the same as this one
# glmm nested fam effects -------------------------------------------------

df2 <- data.frame(t(df))

meta = brainsampledata
lme4_results_table2 <- as.data.frame( matrix( nrow =15694,  ncol = 15) )

for(i in 1:ncol(df2)){
  focal_gene = df2[,i]
  all_others = rowSums(df2[,-i])
  Y = cbind(focal_gene, 
            all_others)
  Model2 = glmer(Y ~ pop * eco * temp + (1|eco/fam),
                family = binomial(),
                # REML = T,
                control = glmerControl(
                  optimizer ='optimx', optCtrl=list(method='nlminb')),
                data = meta)
  # Model = glmer(Y ~ eco*temp + (1|fam),
  #               family = binomial(),
  #               # REML = T,
  #               control = glmerControl(
  #                 optimizer ='optimx', optCtrl=list(method='nlminb')),
  #               data = meta)
  gene_results3 = summary(Model2)
  # mod_coef = coef(summary(Model))
  ## dont need the second one. It does the same as summary
  # gene_results2 = anova(Model,
  #                       test = 'LRT')
  lme4_results_table2[i, 1] = colnames(df2)[i]
  lme4_results_table2[i, 2] = sum(df2[,i]/sum(df2))
  lme4_results_table2[i, 3] = gene_results3$coefficients[1,4]
  lme4_results_table2[i, 4] = gene_results3$coefficients[2,4]
  lme4_results_table2[i, 5] = gene_results3$coefficients[3,4]
  lme4_results_table2[i, 6] = gene_results3$coefficients[4,4]
  lme4_results_table2[i, 7] = gene_results3$coefficients[5,4]
  lme4_results_table2[i, 8] = gene_results3$coefficients[6,4]
  lme4_results_table2[i, 9] = gene_results3$coefficients[7,4]
  lme4_results_table2[i, 10] = gene_results3$coefficients[8,4]
  lme4_results_table2[i, 11] = gene_results3$coefficients[9,4]
  lme4_results_table2[i, 12] = gene_results3$coefficients[10,4]
  # lme4_results_table[i, 13] = gene_results2$coefficients[11,4]
}

gene_results3


# head(results.table)
# tail(results.table)                   
names(lme4_results_table2) <- c('gene_name', 
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

# sort results table to highlight significant cases at the top
# results.table <- results.table[order(results.table$pvalue, decreasing = F) , ]
# dim(results.table)

## in this model only ECOW, popMYV*ECOW, and popSKR*ECOW is significant
## No temperature related treatment effects are significant
## have not accounted for family effects

lme4_results_table[,1:12] %>% 
  as_tibble() %>%
  filter(ecow_temp18_pval < 0.01) 
# write_csv('GLMER_gene_expression_fam_rand_effect.csv')





# liver data --------------------------------------------------------------

setwd('/media/shaun/2d49aafa-914e-40e6-a1ea-142b85301ef9/rna_aligned/f1_labfish/big_batch_2/tsv files/liver/data_analysis/')
liver_sample_data = read_csv('F1_labfish_sampledata10M_liver.csv')
liver_count_data = read_tsv('F1_lab_liver_over10Mreads_alltsv2.tsv') %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID') %>%
  t() %>%
  as_tibble()

liver_orig = read_tsv('F1_lab_liver_over10Mreads_alltsv2.tsv')
# liver remove invariant site ---------------------------------------------

full_data = bind_cols(liver_sample_data, 
                     liver_count_data)

gene_var = full_data %>%
  select(pop, 
         ecotype,
         temp,
         starts_with('ENSGAC')) %>% 
  group_by(pop, 
           ecotype, 
           temp) %>%
  summarise(across(everything(), var, na.rm = F))%>%
  ungroup() %>% 
  select(starts_with('ENSGAC')) %>% 
  colSums(abs(.)>0) %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  rename(variance = 2) %>% 
  filter(variance > 0.000000e+00) %>% 
  rename(GeneID = rowname)

# liver_orig_count_data = liver_count_data %>%
#   rownames_to_column() %>%
#   as_tibble()

liver_cleand_data = left_join(gene_var, 
                        liver_orig, 
                        by = 'GeneID') %>% 
  select(-variance)

liver_df = liver_cleand_data %>% 
  as.data.frame() %>% 
  remove_rownames() %>% 
  column_to_rownames(var = 'GeneID')


# liver glmer  ------------------------------------------------------------
liver_df2 <- data.frame(t(liver_df))

liver_lme4_results <- as.data.frame( matrix( nrow =14566,  ncol = 15) )

for(i in 1:ncol(liver_df2)){
  focal_gene = liver_df2[,i]
  all_others = rowSums(liver_df2[,-i])
  Y = cbind(focal_gene, 
            all_others)
  Model = glmer(Y ~ pop * ecotype * temp + (1|fam),
                family = binomial(),
                # REML = T,
                control = glmerControl(
                  optimizer ='optimx', optCtrl=list(method='nlminb')),
                data = liver_sample_data)
  # Model = glmer(Y ~ eco*temp + (1|fam),
  #               family = binomial(),
  #               # REML = T,
  #               control = glmerControl(
  #                 optimizer ='optimx', optCtrl=list(method='nlminb')),
  #               data = meta)
  liver_glmer_results = summary(Model)
  # mod_coef = coef(summary(Model))
  ## dont need the second one. It does the same as summary
  # gene_results2 = anova(Model,
  #                       test = 'LRT')
  liver_lme4_results[i, 1] = colnames(liver_df2)[i]
  liver_lme4_results[i, 2] = sum(df2[,i]/sum(liver_df2))
  liver_lme4_results[i, 3] = liver_glmer_results $coefficients[1,4]
  liver_lme4_results[i, 4] = liver_glmer_results $coefficients[2,4]
  liver_lme4_results[i, 5] = liver_glmer_results $coefficients[3,4]
  liver_lme4_results[i, 6] = liver_glmer_results $coefficients[4,4]
  liver_lme4_results[i, 7] = liver_glmer_results $coefficients[5,4]
  liver_lme4_results[i, 8] = liver_glmer_results $coefficients[6,4]
  liver_lme4_results[i, 9] = liver_glmer_results $coefficients[7,4]
  liver_lme4_results[i, 10] = liver_glmer_results $coefficients[8,4]
  liver_lme4_results[i, 11] = liver_glmer_results $coefficients[9,4]
  liver_lme4_results[i, 12] = liver_glmer_results $coefficients[10,4]
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
  write_csv('GLMER_gene_expression_ecotemp_pval0.01.csv')
  # write_csv('GLMER_gene_expression_fam_rand_effect.csv')


## The glmm with nested effects is the same as this one