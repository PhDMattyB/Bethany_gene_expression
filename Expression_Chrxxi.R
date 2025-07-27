##############################
## Gene expression differences chr21 inversion
##
## Matt Brachmann (PhDMattyB)
##
## 26.07.2025
##
##############################


library(tidyverse)

setwd('~/Parsons_Postdoc/Stickleback_Genomic/Stickleback_Annotation_features/')

brain_f1_exp_INVT = read_csv('BRAIN_BH_FDR_GLMER_ecotemp_pval0.05.csv') %>% 
  filter(chromosome == 'chrXXI') %>% 
  arrange(position) %>% 
  slice(7:12) %>% 
  mutate(mean_expression_log = log10(mean_expression_relative), 
         status = 'INVERSION')

brain_f1_exp_OUTVT = read_csv('BRAIN_BH_FDR_GLMER_ecotemp_pval0.05.csv') %>% 
  filter(chromosome == 'chrXXI') %>% 
  arrange(position) %>% 
  slice(1:6, 
        12:13) %>% 
  mutate(mean_expression_log = log10(mean_expression_relative), 
         status = 'OUTSIDE INVERSION')

brain_exp_xxi = bind_rows(brain_f1_exp_INVT, 
                          brain_f1_exp_OUTVT)

ggplot(data = brain_exp_xxi, 
       aes(x = position, 
           y = mean_expression_log, 
           col = status))+
  geom_point()+
  scale_y_reverse()

liver_f1_exp_INVT = read_csv('LIVER_BH_FDR_GLMER_ecotemp_pval0.05.csv') %>% 
  filter(chromosome == 'chrXXI') %>% 
  arrange(position) %>% 
  filter(position > 10023795.0, 
         position < 12011713.0) %>% 
  mutate(mean_expression_log = log10(mean_expression_relative), 
         status = 'INVERSION')
liver_f1_exp_OUTVT = read_csv('LIVER_BH_FDR_GLMER_ecotemp_pval0.05.csv') %>% 
  filter(chromosome == 'chrXXI') %>% 
  arrange(position) %>% 
  filter(position < 10023795.0, 
         position < 12011713.0) %>% 
  mutate(mean_expression_log = log10(mean_expression_relative), 
         status = 'OUTSIDE INVERSION')

liver_exp_xxi = bind_rows(liver_f1_exp_INVT, 
                          liver_f1_exp_OUTVT)

ggplot(data = liver_exp_xxi, 
       aes(x = position, 
           y = mean_expression_log, 
           col = status))+
  geom_point()+
  scale_y_reverse()
