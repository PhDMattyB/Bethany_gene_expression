##############################
## Bethany Gene expression
## SKR HYBRIDS
## SEX EFFECTS - MOLECOL Review 1
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

sex_metadata = read_csv('Sex_metadata.csv')%>%
  select(-...7, 
         -...8,
         -...9)
sex_metadata$temp = as.character(sex_metadata$temp)
sex_metadata$sample = as.character(sex_metadata$sample)
sex_data = sex_metadata %>%
  select(sex)

metadata = names(brain_exp) %>% 
  as_tibble() %>% 
  dplyr::slice(-1) %>% 
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
  tidyr::unite(col = ecotemp, 
               c('ecotype',
                 'temp'),
               sep = '_',
               remove = F)%>%
  bind_cols(.,
            sex_data)


metadata %>% 
  group_by(ecotype, 
           temp, 
           sex) %>% 
  summarize(n = n())



# Brain - No sex effect models ----------------------------------------------------
brain_dge_list = DGEList(brain_exp)
brain_norm = calcNormFactors(brain_dge_list)

mm2 = model.matrix(~0 + ecotemp, 
                   data = metadata)

brain_keep = filterByExpr(brain_norm, 
                          min.count = 10,
                          mm2)
sum(brain_keep) # number of genes retai
brain_keep = brain_norm[brain_keep,]

## EdgeR model
brain2_dispersion = estimateDisp(brain_keep, 
                                 mm2)
contrast2 = makeContrasts(eco12 = ecotempSKRC_12 - ecotempSKRW_12, 
                          eco18 = ecotempSKRC_18 - ecotempSKRW_18,
                          plast_amb = ecotempSKRC_12 - ecotempSKRC_18, 
                          plast_geo = ecotempSKRW_12 - ecotempSKRW_18, 
                          plast_hyb = ecotempSKRHYB_12 - ecotempSKRHYB_18, 
                          am_hyb_12 = ecotempSKRC_12 - ecotempSKRHYB_12, 
                          am_hyb_18 = ecotempSKRC_18 - ecotempSKRHYB_18, 
                          geo_hyb_12 = ecotempSKRW_12 - ecotempSKRHYB_12, 
                          geo_hyb_18 = ecotempSKRW_18 - ecotempSKRHYB_18,
                          levels = mm2)

brain2_glm_div = glmQLFit(brain2_dispersion, 
                          # contrast = ecotype.div.brain,
                          design = mm2)

brain2_glm_test = glmQLFTest(brain2_glm_div, 
                             contrast = contrast2)

brain2_edger_results = topTags(brain2_glm_test, 
                               n = 13452,
                               adjust.method = 'bonferroni', 
                               p.value = 0.05)

brain_edger_results$table %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  write_csv('Sex_Brain_EdgeR_GLMQLFTest_results.csv')


# ## limma model
brain2_voom = voom(brain_keep, 
                   mm2, 
                   plot = T)
# 
brain2_fit_limma <- limma::lmFit(brain2_voom, 
                                 contrast = contrast2,
                                 design=mm2)
brain2_fit_limma_contrast = contrasts.fit(brain2_fit_limma, 
                                          contrasts = contrast2)


brain2_fit_ebayes = eBayes(brain2_fit_limma_contrast)

topTable(brain_fit_ebayes, coef = "sexM")

brain_nosex_limma_results = topTable(brain2_fit_ebayes, 
                                     n = 13452, 
                                     adjust.method = 'bonferroni', 
                                     p.value = 0.05)


brain_nosex_limma_results_all = topTable(brain2_fit_ebayes, 
                                         n = 13452, 
                                         adjust.method = 'bonferroni')




# sex included models -----------------------------------------------------

# Brain normalization -----------------------------------------------------


brain_dge_list = DGEList(brain_exp)
brain_norm = calcNormFactors(brain_dge_list)

mm = model.matrix(~0 + ecotemp + sex, 
                  data = metadata)

# mm2 = model.matrix(~0 + ecotemp, 
#                   data = metadata)

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
                              n = 13511,
                              adjust.method = 'bonferroni', 
                              p.value = 0.05)

# brain_edger_results$table %>% 
#   as.data.frame() %>% 
#   as_tibble() %>% 
#   write_csv('Sex_Brain_EdgeR_GLMQLFTest_results.csv')


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

topTable(brain_fit_ebayes, coef = "sexM")


brain_limma_results = topTable(brain_fit_ebayes, 
                               n = 13511, 
                               adjust.method = 'bonferroni', 
                               p.value = 0.05)


brain_limma_results_all = topTable(brain_fit_ebayes, 
                                   n = 13511, 
                                   adjust.method = 'bonferroni')

eco12 = topTable(fit = brain_fit_ebayes,
                 coef = which(colnames(brain_fit_ebayes$coefficients) == 'eco12'),
                 adjust.method = 'bonferroni',
                 number = 13511)

eco12 %>%
  as_tibble() %>%
  write_csv("SEX_Brain_eco_div_12.csv")

eco12 %>%
  as_tibble() %>%
  filter(adj.P.Val <= 0.05) %>%
  write_csv("SEX_Brain_eco_div_12_significant.csv")

eco18 = topTable(fit = brain_fit_ebayes,
                 coef = which(colnames(brain_fit_ebayes$coefficients) == 'eco18'),
                 adjust.method = 'bonferroni',
                 number = 13511)

eco18 %>%
  as_tibble() %>%
  write_csv("SEX_Brain_eco_div_18.csv")

eco18 %>%
  as_tibble() %>%
  filter(adj.P.Val <= 0.05) %>%
  write_csv("SEX_Brain_eco_div_18_significant.csv")


plast_amb = topTable(fit = brain_fit_ebayes,
                     coef = which(colnames(brain_fit_ebayes$coefficients) == 'plast_amb'),
                     adjust.method = 'bonferroni',
                     number = 13511)

plast_amb %>%
  as_tibble() %>%
  write_csv("Sex_Brain_ambient_plastic.csv")

plast_amb %>%
  as_tibble() %>%
  filter(adj.P.Val <= 0.05) %>%
  write_csv("Sex_Brain_ambient_plastic_significant.csv")


plast_geo = topTable(fit = brain_fit_ebayes,
                     coef = which(colnames(brain_fit_ebayes$coefficients) == 'plast_geo'),
                     adjust.method = 'bonferroni',
                     number = 13511)

plast_geo %>%
  as_tibble() %>%
  write_csv('Sex_Brain_geothermal_plastic.csv')

plast_geo %>%
  as_tibble() %>%
  filter(adj.P.Val <= 0.05) %>%
  write_csv('Sex_Brain_geothermal_plastic_significant.csv')


plast_hyb = topTable(fit = brain_fit_ebayes,
                     coef = which(colnames(brain_fit_ebayes$coefficients) == 'plast_hyb'),
                     adjust.method = 'bonferroni',
                     number = 13511)

plast_hyb %>%
  as_tibble() %>%
  write_csv('Sex_Brain_hybrid_plastic.csv')

plast_hyb %>%
  as_tibble() %>%
  filter(adj.P.Val <= 0.05) %>%
  write_csv('Sex_Brain_hybrid_plastic_significant.csv')


amb_hyb_12 = topTable(fit = brain_fit_ebayes,
                      coef = which(colnames(brain_fit_ebayes$coefficients) == 'am_hyb_12'),
                      adjust.method = 'bonferroni',
                      number = 13511)

amb_hyb_12 %>%
  as_tibble() %>%
  write_csv('Sex_Brain_amb_hyb_12_div.csv')

amb_hyb_12 %>%
  as_tibble() %>%
  filter(adj.P.Val <= 0.05) %>%
  write_csv('Sex_Brain_amb_hyb_12_div_significant.csv')

amb_hyb_18 = topTable(fit = brain_fit_ebayes,
                      coef = which(colnames(brain_fit_ebayes$coefficients) == 'am_hyb_18'),
                      adjust.method = 'bonferroni',
                      number = 13511)

amb_hyb_18 %>%
  as_tibble() %>%
  write_csv('Sex_Brain_amb_hyb_18_div.csv')

amb_hyb_18 %>%
  as_tibble() %>%
  filter(adj.P.Val <= 0.05) %>%
  write_csv('Sex_Brain_amb_hyb_18_div_significant.csv')


geo_hyb_12 = topTable(fit = brain_fit_ebayes,
                      coef = which(colnames(brain_fit_ebayes$coefficients) == 'geo_hyb_12'),
                      adjust.method = 'bonferroni',
                      number = 13511)

geo_hyb_12 %>%
  as_tibble() %>%
  write_csv('Sex_Brain_geo_hyb_12_div.csv')

geo_hyb_12 %>%
  as_tibble() %>%
  filter(adj.P.Val <= 0.05) %>%
  write_csv('Sex_Brain_geo_hyb_12_div_significant.csv')


geo_hyb_18 = topTable(fit = brain_fit_ebayes,
                      coef = which(colnames(brain_fit_ebayes$coefficients) == 'geo_hyb_18'),
                      adjust.method = 'bonferroni',
                      number = 13511)

geo_hyb_18 %>%
  as_tibble() %>%
  write_csv('Sex_Brain_geo_hyb_18_div.csv')

geo_hyb_18 %>%
  as_tibble() %>%
  filter(adj.P.Val <= 0.05) %>%
  write_csv('Sex_Brain_geo_hyb_18_div_significant.csv')


brain_limma_results %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  write_csv('Sex_Brain_LIMMA_model_results.csv')

brain_limma_results_all %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  write_csv('Sex_Brain_LIMMA_model_results_all.csv')

# liver normalization -----------------------------------------------------
liver_dge_list = DGEList(liver_exp)
liver_norm = calcNormFactors(liver_dge_list)

mm = model.matrix(~0 + ecotemp + sex, 
                  data = metadata)

liver_keep = filterByExpr(liver_norm, 
                          min.count = 10,
                          mm)
sum(liver_keep) # number of genes retain
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
                              n = 10590,
                              adjust.method = 'bonferroni', 
                              p.value = 0.05)

liver_edger_results$table %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  write_csv('Sex_liver_EdgeR_GLMQLFTest_results.csv')


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

Liver_limma_results_all = topTable(liver_fit_ebayes, 
                                   n = 10590, 
                                   adjust.method = 'bonferroni')

liver_eco12 = topTable(fit = liver_fit_ebayes,
                       coef = which(colnames(liver_fit_ebayes$coefficients) == 'eco12'),
                       adjust.method = 'bonferroni',
                       number = 10590)

liver_eco12 %>%
  as_tibble() %>%
  write_csv("Sex_Liver_eco_div_12.csv")

liver_eco12 %>%
  as_tibble() %>%
  filter(adj.P.Val <= 0.05) %>%
  write_csv("Sex_Liver_eco_div_12_significant.csv")



liver_eco18 = topTable(fit = liver_fit_ebayes,
                       coef = which(colnames(liver_fit_ebayes$coefficients) == 'eco18'),
                       adjust.method = 'bonferroni',
                       number = 10590)

liver_eco18 %>%
  as_tibble() %>%
  write_csv("Sex_Liver_eco_div_18.csv")
liver_eco18 %>%
  as_tibble() %>%
  filter(adj.P.Val <= 0.05) %>%
  write_csv("Sex_Liver_eco_div_18_significant.csv")


liver_plast_amb = topTable(fit = liver_fit_ebayes,
                           coef = which(colnames(liver_fit_ebayes$coefficients) == 'plast_amb'),
                           adjust.method = 'bonferroni',
                           number = 10590)

liver_plast_amb %>%
  as_tibble() %>%
  write_csv("Sex_Liver_ambient_plastic.csv")

liver_plast_amb %>%
  as_tibble() %>%
  filter(adj.P.Val <= 0.05) %>%
  write_csv("Sex_Liver_ambient_plastic_significant.csv")


liver_plast_geo = topTable(fit = liver_fit_ebayes,
                           coef = which(colnames(liver_fit_ebayes$coefficients) == 'plast_geo'),
                           adjust.method = 'bonferroni',
                           number = 10590)

liver_plast_geo %>%
  as_tibble() %>%
  write_csv('Sex_Liver_geothermal_plastic.csv')
liver_plast_geo %>%
  as_tibble() %>%
  filter(adj.P.Val <= 0.05) %>%
  write_csv('Sex_Liver_geothermal_plastic_significant.csv')

liver_plast_hyb = topTable(fit = liver_fit_ebayes,
                           coef = which(colnames(liver_fit_ebayes$coefficients) == 'plast_hyb'),
                           adjust.method = 'bonferroni',
                           number = 10590)

liver_plast_hyb %>%
  as_tibble() %>%
  write_csv('Sex_Liver_hybrid_plastic.csv')

liver_plast_hyb %>%
  as_tibble() %>%
  filter(adj.P.Val <= 0.05) %>%
  write_csv('Sex_Liver_hybrid_plastic_significant.csv')


liver_amb_hyb_12 = topTable(fit = liver_fit_ebayes,
                            coef = which(colnames(liver_fit_ebayes$coefficients) == 'am_hyb_12'),
                            adjust.method = 'bonferroni',
                            number = 10590)

liver_amb_hyb_12 %>%
  as_tibble() %>%
  write_csv('Sex_liver_amb_hyb_12_div.csv')

liver_amb_hyb_12 %>%
  as_tibble() %>%
  filter(adj.P.Val <= 0.05) %>%
  write_csv('Sex_liver_amb_hyb_12_div_significant.csv')


liver_amb_hyb_18 = topTable(fit = liver_fit_ebayes,
                            coef = which(colnames(liver_fit_ebayes$coefficients) == 'am_hyb_18'),
                            adjust.method = 'bonferroni',
                            number = 10590)

liver_amb_hyb_18 %>%
  as_tibble() %>%
  write_csv('Sex_liver_amb_hyb_18_div.csv')

liver_amb_hyb_18 %>%
  as_tibble() %>%
  filter(adj.P.Val <= 0.05) %>%
  write_csv('Sex_liver_amb_hyb_18_div_significant.csv')


liver_geo_hyb_12 = topTable(fit = liver_fit_ebayes,
                            coef = which(colnames(liver_fit_ebayes$coefficients) == 'geo_hyb_12'),
                            adjust.method = 'bonferroni',
                            number = 10590)

liver_geo_hyb_12 %>%
  as_tibble() %>%
  write_csv('Sex_liver_geo_hyb_12_div.csv')

liver_geo_hyb_12 %>%
  as_tibble() %>%
  filter(adj.P.Val <= 0.05) %>%
  write_csv('Sex_liver_geo_hyb_12_div_significant.csv')


liver_geo_hyb_18 = topTable(fit = liver_fit_ebayes,
                            coef = which(colnames(liver_fit_ebayes$coefficients) == 'geo_hyb_18'),
                            adjust.method = 'bonferroni',
                            number = 10590)

liver_geo_hyb_18 %>%
  as_tibble() %>%
  write_csv('Sex_liver_geo_hyb_18_div.csv')

liver_geo_hyb_18 %>%
  as_tibble() %>%
  filter(adj.P.Val <= 0.05) %>%
  write_csv('Sex_liver_geo_hyb_18_div_significant.csv')




Liver_limma_results_all %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  write_csv('Sex_Liver_LIMMA_model_results_all.csv')




liver_limma_results = topTable(liver_fit_ebayes, 
                               n = 10590, 
                               adjust.method = 'bonferroni', 
                               p.value = 0.05)

liver_limma_results %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  write_csv('Sex_liver_LIMMA_model_results.csv')
