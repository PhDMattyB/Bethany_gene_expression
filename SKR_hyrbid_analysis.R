##############################
## Bethany Gene expression
## SKR HYBRIDS
##
## Matt Brachmann (PhDMattyB)
##
## 22.05.2025
##
##############################


setwd('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/')

library(tidyverse)
library(edgeR)

theme_set(theme_bw())

brain_exp = read_tsv('brain_gene_read_counts_table_all_final.tsv')
liver_exp = read_tsv('liver_gene_read_counts_table_all_final.tsv')

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


brain_metadata %>% 
  group_by(ecotype, 
           temp) %>% 
  summarize(n = n())



# Brain normalization -----------------------------------------------------


brain_dge_list = DGEList(brain_exp)
brain_norm = calcNormFactors(brain_dge_list)

mm = model.matrix(~0 + ecotemp, 
                  data = metadata)

brain_keep = filterByExpr(brain_norm, 
                          min.count = 10,
                          mm)
sum(brain_keep) # number of genes retai
brain_keep = brain_norm[brain_keep,]

## EdgeR model
brain_dispersion = estimateDisp(brain_keep, 
                                mm) 
# Brain divergent gene expression -----------------------------------------

contrast = makeContrasts(eco12 = ecotempSKRC_12 - ecotempSKRW_12, 
                             eco18 = ecotempSKRC_18 - ecotempSKRW_18,
                             plast_amb = ecotempSKRC_12 - ecotempSKRC_18, 
                             plast_geo = ecotempSKRW_12 - ecotempSKRW_18, 
                             plast_hyb = ecotempSKRHYB_12 - ecotempSKRHYB_18, 
                             am_hyb_12 = ecotempSKRC_12 - ecotempSKRHYB_12, 
                             am_hyb_18 = ecotempSKRC_18 - ecotempSKRHYB_18, 
                             geo_hyb_12 = ecotempSKRW_12 - ecotempSKRHYB_12, 
                             geo_hyb_18 = ecotempSKRW_18 - ecotempSKRHYB_18,
                                  levels = mm)

brain_glm_div = glmQLFit(brain_dispersion, 
                     # contrast = ecotype.div.brain,
                     design = mm)

brain_glm_test = glmQLFTest(brain_glm_div, 
                                contrast = contrast)

brain_edger_results = topTags(brain_glm_test, 
        n = 13452,
        adjust.method = 'bonferroni', 
        p.value = 0.05)

brain_edger_results$table %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  write_csv('Brain_EdgeR_GLMQLFTest_results.csv')


# ## limma model
brain_voom = voom(brain_keep, 
                  mm, 
                  plot = T)
# 
brain_fit_limma <- limma::lmFit(brain_voom, 
                                contrast = contrast,
                          design=mm)
brain_fit_limma_contrast = contrasts.fit(brain_fit_limma, 
              contrasts = contrast)


brain_fit_ebayes = eBayes(brain_fit_limma_contrast)

brain_limma_results = topTable(brain_fit_ebayes, 
         n = 13452, 
         adjust.method = 'bonferroni', 
         p.value = 0.05)

brain_limma_results %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  write_csv('Brain_LIMMA_model_results.csv')


# liver normalization -----------------------------------------------------
liver_dge_list = DGEList(liver_exp)
liver_norm = calcNormFactors(liver_dge_list)

mm = model.matrix(~0 + ecotemp, 
                  data = metadata)

liver_keep = filterByExpr(liver_norm, 
                          min.count = 10,
                          mm)
sum(liver_keep) # number of genes retai
liver_keep = liver_norm[liver_keep,]

## EdgeR model
liver_dispersion = estimateDisp(liver_keep, 
                                mm) 
# Liver divergent gene expression -----------------------------------------

contrast = makeContrasts(eco12 = ecotempSKRC_12 - ecotempSKRW_12, 
                         eco18 = ecotempSKRC_18 - ecotempSKRW_18,
                         plast_amb = ecotempSKRC_12 - ecotempSKRC_18, 
                         plast_geo = ecotempSKRW_12 - ecotempSKRW_18, 
                         plast_hyb = ecotempSKRHYB_12 - ecotempSKRHYB_18, 
                         am_hyb_12 = ecotempSKRC_12 - ecotempSKRHYB_12, 
                         am_hyb_18 = ecotempSKRC_18 - ecotempSKRHYB_18, 
                         geo_hyb_12 = ecotempSKRW_12 - ecotempSKRHYB_12, 
                         geo_hyb_18 = ecotempSKRW_18 - ecotempSKRHYB_18,
                         levels = mm)

liver_glm_div = glmQLFit(liver_dispersion, 
                         # contrast = ecotype.div.liver,
                         design = mm)

liver_glm_test = glmQLFTest(liver_glm_div, 
                            contrast = contrast)

liver_edger_results = topTags(liver_glm_test, 
                              n = 13452,
                              adjust.method = 'bonferroni', 
                              p.value = 0.05)

liver_edger_results$table %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  write_csv('liver_EdgeR_GLMQLFTest_results.csv')


# ## limma model
liver_voom = voom(liver_keep, 
                  mm, 
                  plot = T)
# 
liver_fit_limma <- limma::lmFit(liver_voom, 
                                contrast = contrast,
                                design=mm)
liver_fit_limma_contrast = contrasts.fit(liver_fit_limma, 
                                         contrasts = contrast)


liver_fit_ebayes = eBayes(liver_fit_limma_contrast)

liver_limma_results = topTable(liver_fit_ebayes, 
                               n = 13452, 
                               adjust.method = 'bonferroni', 
                               p.value = 0.05)

liver_limma_results %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  write_csv('liver_LIMMA_model_results.csv')



# Brain gene expression - downstream --------------------------------------

brain_edge = read_csv('Brain_EdgeR_GLMQLFTest_results.csv')
brain_limma = read_csv('Brain_LIMMA_model_results.csv')

brain_common_genes = inner_join(brain_edge, 
           brain_limma, 
           by = 'GeneID')


### 12 degree inheritance pattern
exp_pattern_12 = brain_common_genes %>% 
  select(GeneID, 
         am_hyb_12, 
         geo_hyb_12)
Nothing = exp_pattern_12 %>% 
  filter(am_hyb_12 <= 0.32 & am_hyb_12 >= -0.32) %>% 
  filter(geo_hyb_12 < 0.32 & geo_hyb_12 > -0.32) %>% 
  mutate(exp_pattern = 'Neutral')
## 987 genes
dom_geo1 = exp_pattern_12 %>%
  filter(am_hyb_12 <= 0.32 &
         am_hyb_12 >= -0.32) %>%
  filter(geo_hyb_12 > 0.32) %>% 
  mutate(exp_pattern = 'Dominant Geothermal')

dom_geo2 = exp_pattern_12 %>%
  filter(am_hyb_12 <= 0.32 &
           am_hyb_12 >= -0.32) %>%
  filter(geo_hyb_12 < -0.32) %>% 
  mutate(exp_pattern = 'Dominant Geothermal')


dom_amb1 = exp_pattern_12 %>% 
  filter(geo_hyb_12 <= 0.32 & geo_hyb_12 >= -0.32) %>%
  filter(am_hyb_12 > 0.32) %>% 
  mutate(exp_pattern = 'Dominant Ambient')
dom_amb2 = exp_pattern_12 %>% 
  filter(geo_hyb_12 <= 0.32 & geo_hyb_12 >= 0.32) %>%
  filter(am_hyb_12 < -0.32) %>% 
    mutate(exp_pattern = 'Dominant Ambient')

additive1 = exp_pattern_12 %>% 
  filter(am_hyb_12 > 0.32 & geo_hyb_12 < -0.32) %>% 
  mutate(exp_pattern = 'Additive')
additive2 = exp_pattern_12 %>% 
  filter(am_hyb_12 < -0.32 & geo_hyb_12 > 0.32) %>% 
  mutate(exp_pattern = 'Additive')

transgressive1 = exp_pattern_12 %>% 
  filter(am_hyb_12 > 0.32 & geo_hyb_12 > 0.32) %>% 
  mutate(exp_pattern = 'Transgressive')

transgressive2 = exp_pattern_12 %>% 
  filter(am_hyb_12 < -0.32 & geo_hyb_12 < -0.32) %>% 
  mutate(exp_pattern = 'Transgressive')

exp_pattern_12_graph = bind_rows(dom_amb1, 
                                 dom_amb2, 
                                 dom_geo1,
                           dom_geo2, 
                           additive1, 
                           additive2, 
                           transgressive1, 
                           transgressive2, 
                           Nothing)

inheritance_pal = c('#ff9e00',
                    '#669bbc',
                    '#c1121f', 
                    '#adb5bd', 
                    '#7400b8')


Inheritance_pattern_12 = ggplot(data = exp_pattern_12_graph, 
       aes(x = geo_hyb_12, 
           y = am_hyb_12, 
           col = exp_pattern))+
  geom_point()+
  geom_hline(yintercept = 0, 
             col = 'black')+
  geom_vline(xintercept = 0, 
             col = 'black')+
  geom_hline(yintercept = 0.32, 
             linetype = 'dashed')+
  geom_hline(yintercept = -0.32, 
             linetype = 'dashed')+
  geom_vline(xintercept = 0.32, 
             linetype = 'dashed')+
  geom_vline(xintercept = -0.32, 
             linetype = 'dashed')+
  labs(x = 'Geothermal 12 - Hybrid 12', 
       y = 'Ambient 12 - Hybrid 12',
       title = 'A)')+
  xlim(-3.5, 3.5)+
  ylim(-3.5, 3.5)+
  scale_color_manual(values = inheritance_pal)+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')


### 18 degree inheritance pattern
exp_pattern_18 = brain_common_genes %>% 
  select(GeneID, 
         am_hyb_18, 
         geo_hyb_18)
Nothing = exp_pattern_18 %>% 
  filter(am_hyb_18 <= 0.32 & am_hyb_18 >= -0.32) %>% 
  filter(geo_hyb_18 < 0.32 & geo_hyb_18 > -0.32) %>% 
  mutate(exp_pattern = 'Neutral')
## 987 genes
dom_geo1 = exp_pattern_18 %>%
  filter(am_hyb_18 <= 0.32 &
           am_hyb_18 >= -0.32) %>%
  filter(geo_hyb_18 > 0.32) %>% 
  mutate(exp_pattern = 'Dominant Geothermal')

dom_geo2 = exp_pattern_18 %>%
  filter(am_hyb_18 <= 0.32 &
           am_hyb_18 >= -0.32) %>%
  filter(geo_hyb_18 < -0.32) %>% 
  mutate(exp_pattern = 'Dominant Geothermal')


dom_amb1 = exp_pattern_18 %>% 
  filter(geo_hyb_18 <= 0.32 & geo_hyb_18 >= -0.32) %>%
  filter(am_hyb_18 > 0.32) %>% 
  mutate(exp_pattern = 'Dominant Ambient')
dom_amb2 = exp_pattern_18 %>% 
  filter(geo_hyb_18 <= 0.32 & geo_hyb_18 >= 0.32) %>%
  filter(am_hyb_18 < -0.32) %>% 
  mutate(exp_pattern = 'Dominant Ambient')

additive1 = exp_pattern_18 %>% 
  filter(am_hyb_18 > 0.32 & geo_hyb_18 < -0.32) %>% 
  mutate(exp_pattern = 'Additive')
additive2 = exp_pattern_18 %>% 
  filter(am_hyb_18 < -0.32 & geo_hyb_18 > 0.32) %>% 
  mutate(exp_pattern = 'Additive')

transgressive1 = exp_pattern_18 %>% 
  filter(am_hyb_18 > 0.32 & geo_hyb_18 > 0.32) %>% 
  mutate(exp_pattern = 'Transgressive')

transgressive2 = exp_pattern_18 %>% 
  filter(am_hyb_18 < -0.32 & geo_hyb_18 < -0.32) %>% 
  mutate(exp_pattern = 'Transgressive')

exp_pattern_18_graph = bind_rows(dom_amb1, 
                                 dom_amb2, 
                                 dom_geo1,
                                 dom_geo2, 
                                 additive1, 
                                 additive2, 
                                 transgressive1, 
                                 transgressive2, 
                                 Nothing)

inheritance_pal = c('#ff9e00',
                    '#669bbc',
                    '#c1121f', 
                    '#adb5bd', 
                    '#7400b8')


Inheritance_pattern_18 = ggplot(data = exp_pattern_18_graph, 
       aes(x = geo_hyb_18, 
           y = am_hyb_18, 
           col = exp_pattern))+
  geom_point()+
  geom_hline(yintercept = 0, 
             col = 'black')+
  geom_vline(xintercept = 0, 
             col = 'black')+
  geom_hline(yintercept = 0.32, 
             linetype = 'dashed')+
  geom_hline(yintercept = -0.32, 
             linetype = 'dashed')+
  geom_vline(xintercept = 0.32, 
             linetype = 'dashed')+
  geom_vline(xintercept = -0.32, 
             linetype = 'dashed')+
  labs(x = 'Geothermal 12 - Hybrid 12', 
       y = 'Ambient 12 - Hybrid 12', 
       title = 'B)')+
  xlim(-3.5, 3.5)+
  ylim(-3.5, 3.5)+
  scale_color_manual(values = inheritance_pal)+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')



# Liver inheritance pattern -----------------------------------------------


liver_edge = read_csv('Liver_EdgeR_GLMQLFTest_results.csv')
liver_limma = read_csv('Liver_LIMMA_model_results.csv')

liver_common_genes = inner_join(liver_edge, 
                                liver_limma, 
                                by = 'GeneID')


### 12 degree inheritance pattern
liver_pattern_12 = liver_common_genes %>% 
  select(GeneID, 
         am_hyb_12, 
         geo_hyb_12)
Nothing = liver_pattern_12 %>% 
  filter(am_hyb_12 <= 0.32 & am_hyb_12 >= -0.32) %>% 
  filter(geo_hyb_12 < 0.32 & geo_hyb_12 > -0.32) %>% 
  mutate(liver_pattern = 'Neutral')
## 987 genes
dom_geo1 = liver_pattern_12 %>%
  filter(am_hyb_12 <= 0.32 &
           am_hyb_12 >= -0.32) %>%
  filter(geo_hyb_12 > 0.32) %>% 
  mutate(liver_pattern = 'Dominant Geothermal')

dom_geo2 = liver_pattern_12 %>%
  filter(am_hyb_12 <= 0.32 &
           am_hyb_12 >= -0.32) %>%
  filter(geo_hyb_12 < -0.32) %>% 
  mutate(liver_pattern = 'Dominant Geothermal')


dom_amb1 = liver_pattern_12 %>% 
  filter(geo_hyb_12 <= 0.32 & geo_hyb_12 >= -0.32) %>%
  filter(am_hyb_12 > 0.32) %>% 
  mutate(liver_pattern = 'Dominant Ambient')
dom_amb2 = liver_pattern_12 %>% 
  filter(geo_hyb_12 <= 0.32 & geo_hyb_12 >= 0.32) %>%
  filter(am_hyb_12 < -0.32) %>% 
  mutate(liver_pattern = 'Dominant Ambient')

additive1 = liver_pattern_12 %>% 
  filter(am_hyb_12 > 0.32 & geo_hyb_12 < -0.32) %>% 
  mutate(liver_pattern = 'Additive')
additive2 = liver_pattern_12 %>% 
  filter(am_hyb_12 < -0.32 & geo_hyb_12 > 0.32) %>% 
  mutate(liver_pattern = 'Additive')

transgressive1 = liver_pattern_12 %>% 
  filter(am_hyb_12 > 0.32 & geo_hyb_12 > 0.32) %>% 
  mutate(liver_pattern = 'Transgressive')

transgressive2 = liver_pattern_12 %>% 
  filter(am_hyb_12 < -0.32 & geo_hyb_12 < -0.32) %>% 
  mutate(liver_pattern = 'Transgressive')

liver_pattern_12_graph = bind_rows(dom_amb1, 
                                 dom_amb2, 
                                 dom_geo1,
                                 dom_geo2, 
                                 additive1, 
                                 additive2, 
                                 transgressive1, 
                                 transgressive2, 
                                 Nothing)

inheritance_pal = c('#ff9e00',
                    '#669bbc',
                    '#c1121f', 
                    '#adb5bd', 
                    '#7400b8')


liver_Inheritance_pattern_12 = ggplot(data = liver_pattern_12_graph, 
                                aes(x = geo_hyb_12, 
                                    y = am_hyb_12, 
                                    col = liver_pattern))+
  geom_point()+
  geom_hline(yintercept = 0, 
             col = 'black')+
  geom_vline(xintercept = 0, 
             col = 'black')+
  geom_hline(yintercept = 0.32, 
             linetype = 'dashed')+
  geom_hline(yintercept = -0.32, 
             linetype = 'dashed')+
  geom_vline(xintercept = 0.32, 
             linetype = 'dashed')+
  geom_vline(xintercept = -0.32, 
             linetype = 'dashed')+
  labs(x = 'Geothermal 12 - Hybrid 12', 
       y = 'Ambient 12 - Hybrid 12',
       title = 'C)')+
  xlim(-3.5, 3.5)+
  ylim(-3.5, 3.5)+
  scale_color_manual(values = inheritance_pal)+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')


### 18 degree inheritance pattern
liver_pattern_18 = liver_common_genes %>% 
  select(GeneID, 
         am_hyb_18, 
         geo_hyb_18)
Nothing = liver_pattern_18 %>% 
  filter(am_hyb_18 <= 0.32 & am_hyb_18 >= -0.32) %>% 
  filter(geo_hyb_18 < 0.32 & geo_hyb_18 > -0.32) %>% 
  mutate(liver_pattern = 'Neutral')
## 987 genes
dom_geo1 = liver_pattern_18 %>%
  filter(am_hyb_18 <= 0.32 &
           am_hyb_18 >= -0.32) %>%
  filter(geo_hyb_18 > 0.32) %>% 
  mutate(liver_pattern = 'Dominant Geothermal')

dom_geo2 = liver_pattern_18 %>%
  filter(am_hyb_18 <= 0.32 &
           am_hyb_18 >= -0.32) %>%
  filter(geo_hyb_18 < -0.32) %>% 
  mutate(liver_pattern = 'Dominant Geothermal')


dom_amb1 = liver_pattern_18 %>% 
  filter(geo_hyb_18 <= 0.32 & geo_hyb_18 >= -0.32) %>%
  filter(am_hyb_18 > 0.32) %>% 
  mutate(liver_pattern = 'Dominant Ambient')
dom_amb2 = liver_pattern_18 %>% 
  filter(geo_hyb_18 <= 0.32 & geo_hyb_18 >= 0.32) %>%
  filter(am_hyb_18 < -0.32) %>% 
  mutate(liver_pattern = 'Dominant Ambient')

additive1 = liver_pattern_18 %>% 
  filter(am_hyb_18 > 0.32 & geo_hyb_18 < -0.32) %>% 
  mutate(liver_pattern = 'Additive')
additive2 = liver_pattern_18 %>% 
  filter(am_hyb_18 < -0.32 & geo_hyb_18 > 0.32) %>% 
  mutate(liver_pattern = 'Additive')

transgressive1 = liver_pattern_18 %>% 
  filter(am_hyb_18 > 0.32 & geo_hyb_18 > 0.32) %>% 
  mutate(liver_pattern = 'Transgressive')

transgressive2 = liver_pattern_18 %>% 
  filter(am_hyb_18 < -0.32 & geo_hyb_18 < -0.32) %>% 
  mutate(liver_pattern = 'Transgressive')

liver_pattern_18_graph = bind_rows(dom_amb1, 
                                 dom_amb2, 
                                 dom_geo1,
                                 dom_geo2, 
                                 additive1, 
                                 additive2, 
                                 transgressive1, 
                                 transgressive2, 
                                 Nothing)

liv_inheritance_pal = c('#669bbc',
                    '#c1121f', 
                    '#adb5bd', 
                    '#7400b8')


liver_Inheritance_pattern_18 = ggplot(data = liver_pattern_18_graph, 
                                aes(x = geo_hyb_18, 
                                    y = am_hyb_18, 
                                    col = liver_pattern))+
  geom_point()+
  geom_hline(yintercept = 0, 
             col = 'black')+
  geom_vline(xintercept = 0, 
             col = 'black')+
  geom_hline(yintercept = 0.32, 
             linetype = 'dashed')+
  geom_hline(yintercept = -0.32, 
             linetype = 'dashed')+
  geom_vline(xintercept = 0.32, 
             linetype = 'dashed')+
  geom_vline(xintercept = -0.32, 
             linetype = 'dashed')+
  labs(x = 'Geothermal 12 - Hybrid 12', 
       y = 'Ambient 12 - Hybrid 12', 
       title = 'D)')+
  xlim(-3.5, 3.5)+
  ylim(-3.5, 3.5)+
  scale_color_manual(values = liv_inheritance_pal)+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')



# Combine the inheritance graphs ------------------------------------------

Inherit_deez_nuts = (Inheritance_pattern_12 + Inheritance_pattern_18)/(liver_Inheritance_pattern_12+liver_Inheritance_pattern_18)

ggsave('Hybrid_inheritance_pattern.tiff', 
       plot = Inherit_deez_nuts, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 15)

