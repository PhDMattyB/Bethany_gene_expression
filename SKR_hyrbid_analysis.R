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

# brain_glm = glmQLFit(brain_dispersion, 
#          design = mm, 
#          dispersion = brain_dispersion$common.dispersion)
# 
# brain_glm_test = glmQLFTest(brain_glm)
# 
# topTags(brain_glm_test)
# 
# ## limma model
# brain_voom = voom(brain_keep, mm, plot = T)
# 
# fit_limma <- limma::lmFit(brain_voom, design=mm)


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

topTags(eco_div_brain_test)

eco_div_brain_test$coefficients

# liver normalization -----------------------------------------------------


liver_dge_list = DGEList(liver_exp)
liver_norm = calcNormFactors(liver_dge_list)

mm = model.matrix(~0 + ecotype*temp, 
                  data = brain_metadata)

liver_keep = filterByExpr(liver_norm, 
                          min.count = 10,
                          mm)
sum(liver_keep) # number of genes retai
liver_keep = liver_norm[liver_keep,]

liver_voom = voom(liver_keep, mm, plot = T)

liver_dispersion = estimateDisp(liver_keep)
