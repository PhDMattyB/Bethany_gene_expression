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
