##############################
## Bethany Gene expression
## SKR HYBRIDS
## SEX EFFECTS - MOLECOL Review 1
##
## Matt Brachmann (PhDMattyB)
##
## 29.07.2026
##
##############################

### Comparing model outputs for including and accounting for sex effects


# Brain results -----------------------------------------------------------


brain_sex_effects_eco12 = read_csv('Sex_Brain_eco_div_12_significant.csv') %>% 
  mutate(status = 'Outlier') %>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))

brain_nosex_eco12 = read_csv('Brain_eco_div_12_significant.csv') %>% 
  mutate(status = 'Outlier') %>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))

brain_sex_effects_eco12$GeneID == brain_nosex_eco12$GeneID

brain_sex_effects_eco18 = read_csv('Sex_Brain_eco_div_18_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))

brain_nosex_eco18 = read_csv('Brain_eco_div_18_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))


brain_sex_effects_plast_amb = read_csv('Sex_Brain_ambient_plastic_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
brain_nosex_plast_amb = read_csv('Brain_ambient_plastic_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))


brain_sex_effects_plast_geo = read_csv('Sex_Brain_geothermal_plastic_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
brain_nosex_plast_geo = read_csv('Brain_geothermal_plastic_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))

brain_sex_effects_plast_hyb = read_csv('Sex_Brain_hybrid_plastic_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
brain_nosex_plast_hyb = read_csv('Brain_hybrid_plastic_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))

brain_sex_effects_amb_hyb_div_12 = read_csv('Sex_Brain_amb_hyb_12_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
brain_nosex_amb_hyb_div_12 = read_csv('Brain_amb_hyb_12_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))

brain_sex_effects_amb_hyb_div_18 = read_csv('Sex_Brain_amb_hyb_18_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
brain_nosex_amb_hyb_div_18 = read_csv('Brain_amb_hyb_18_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))


brain_sex_effects_geo_hyb_div_12 = read_csv('Sex_Brain_geo_hyb_12_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
brain_nosex_geo_hyb_div_12 = read_csv('Brain_geo_hyb_12_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))

brain_sex_effects_geo_hyb_div_18 = read_csv('Sex_Brain_geo_hyb_18_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))

brain_nosex_geo_hyb_div_18 = read_csv('Brain_geo_hyb_18_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))


# Liver results -----------------------------------------------------------


liver_sex_effects_eco12 = read_csv('Sex_liver_eco_div_12_significant.csv') %>% 
  mutate(status = 'Outlier') %>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))

liver_nosex_eco12 = read_csv('liver_eco_div_12_significant.csv') %>% 
  mutate(status = 'Outlier') %>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))

liver_sex_effects_eco12$GeneID == liver_nosex_eco12$GeneID

liver_sex_effects_eco18 = read_csv('Sex_liver_eco_div_18_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))

liver_nosex_eco18 = read_csv('liver_eco_div_18_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))


liver_sex_effects_plast_amb = read_csv('Sex_liver_ambient_plastic_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
liver_nosex_plast_amb = read_csv('liver_ambient_plastic_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))


liver_sex_effects_plast_geo = read_csv('Sex_liver_geothermal_plastic_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
liver_nosex_plast_geo = read_csv('liver_geothermal_plastic_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))

liver_sex_effects_plast_hyb = read_csv('Sex_liver_hybrid_plastic_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
liver_nosex_plast_hyb = read_csv('liver_hybrid_plastic_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))

liver_sex_effects_amb_hyb_div_12 = read_csv('Sex_liver_amb_hyb_12_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
liver_nosex_amb_hyb_div_12 = read_csv('liver_amb_hyb_12_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))

liver_sex_effects_amb_hyb_div_18 = read_csv('Sex_liver_amb_hyb_18_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
liver_nosex_amb_hyb_div_18 = read_csv('liver_amb_hyb_18_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))


liver_sex_effects_geo_hyb_div_12 = read_csv('Sex_liver_geo_hyb_12_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
liver_nosex_geo_hyb_div_12 = read_csv('liver_geo_hyb_12_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))

liver_sex_effects_geo_hyb_div_18 = read_csv('Sex_liver_geo_hyb_18_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))

liver_nosex_geo_hyb_div_18 = read_csv('liver_geo_hyb_18_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))


# Variance Partitioning ---------------------------------------------------

library(variancePartition)

brain_dge_list = DGEList(brain_exp)
brain_norm = calcNormFactors(brain_dge_list)

mm2 = model.matrix(~0 + ecotemp + sex, 
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
                               n = 13511,
                               adjust.method = 'bonferroni', 
                               p.value = 0.05)

# brain_edger_results$table %>% 
#   as.data.frame() %>% 
#   as_tibble() %>% 
#   write_csv('Sex_Brain_EdgeR_GLMQLFTest_results.csv')


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

brain_nosex_limma_results = topTable(brain2_fit_ebayes, 
                                     n = 13511, 
                                     adjust.method = 'bonferroni', 
                                     p.value = 0.05)


brain_nosex_limma_results_all = topTable(brain2_fit_ebayes, 
                                         n = 13511, 
                                         adjust.method = 'bonferroni')



# form <- ~ Age + (1 | Individual) + (1 | Tissue) + (1 | Batch)

metadata

form2 = model.matrix(~ecotype * temp + sex, 
                   data = metadata)


brain_keep = filterByExpr(brain_norm, 
                          min.count = 10,
                          form2)
sum(brain_keep) # number of genes retai
brain_keep = brain_norm[brain_keep,]

varpart_voom = voom(brain_keep, 
                   form2, 
                   plot = T)

form = ~ ecotype * temp + sex 

varPart <- fitExtractVarPartModel(varpart_voom, form, metadata)

vp <- sortCols(varPart)

plotPercentBars(vp[1:10, ])

plotVarPart(vp)

results <- fitVarPartModel(varpart_voom, form, metadata)
varPart_res <- extractVarPart(results)
pSummaries <- fitVarPartModel(varpart_voom, form, metadata, summary)

pSummaries[[1]]





brain_dge_list = DGEList(brain_exp)
brain_norm = calcNormFactors(brain_dge_list)
# Specify variables to be included in the voom() estimates of
# uncertainty.
# Recommend including variables with a small number of categories
# that explain a substantial amount of variation
form2 = model.matrix(~ecotype * temp + sex, 
                     data = metadata)

# Estimate precision weights for each gene and sample
# This models uncertainty in expression measurements
vobjGenes <- voom(brain_norm, form2)

form = ~ ecotype * temp + sex 

varPart <- fitExtractVarPartModel(vobjGenes, form, metadata)

vp <- sortCols(varPart)

plotPercentBars(vp[1:10, ])

varpart_cols = c('#005f73', 
                 '#e9d8a6', 
                 '#bb3e03', 
                 '#fb6f92', 
                 '#9f9aa4')

varpart_plot = plotVarPart(vp)+
  scale_fill_manual(values = varpart_cols)

ggsave('Variance_partition_plot.svg', 
       plot = varpart_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20, 
       height = 15)
