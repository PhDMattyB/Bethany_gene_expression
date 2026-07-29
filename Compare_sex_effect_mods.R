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

