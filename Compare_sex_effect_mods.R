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

