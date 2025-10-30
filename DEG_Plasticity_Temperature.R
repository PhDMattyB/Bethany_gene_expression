##############################
## DEG plasticity between temperatures
##
## Matt Brachmann (PhDMattyB)
##
## 27.05.2025
##
##############################

setwd('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/')

library(tidyverse)
library(patchwork)
library(edgeR)

theme_set(theme_bw())

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

## Limma and edger analyses were done in another script
## This is a follow up script from
## the SKR_hybrid_analysis script to look at temperature
## and ecotype effects more effectively


liver_common_genes = read_csv('Liver_DEG_Common_edger_limma.csv')


# Brain temp plast --------------------------------------------------------
brain_limma_all = read_csv('Brain_LIMMA_model_results_all.csv')

brain_limma_all %>% 
  dplyr::select(GeneID) %>% 
  write_tsv('brain_limma_gene_list.txt')

## All genes 

brain_eco12_neutral = read_csv('Brain_eco_div_12.csv') %>% 
  filter(adj.P.Val > 0.05) %>% 
  mutate(status = 'Neutral') %>% 
  mutate(Regulated = 'Neutral')
brain_eco18_neutral = read_csv('Brain_eco_div_18.csv') %>% 
  filter(adj.P.Val > 0.05)%>% 
  mutate(status = 'Neutral')%>% 
  mutate(Regulated = 'Neutral')
brain_plast_amb_neutral = read_csv('Brain_ambient_plastic.csv') %>% 
  filter(adj.P.Val > 0.05)%>% 
  mutate(status = 'Neutral')%>% 
  mutate(Regulated = 'Neutral')
brain_plast_geo_neutral = read_csv('Brain_geothermal_plastic.csv') %>% 
  filter(adj.P.Val > 0.05)%>% 
  mutate(status = 'Neutral')%>% 
  mutate(Regulated = 'Neutral')
brain_plast_hyb_neutral = read_csv('Brain_hybrid_plastic.csv') %>% 
  filter(adj.P.Val > 0.05)%>% 
  mutate(status = 'Neutral')%>% 
  mutate(Regulated = 'Neutral')

brain_amb_hyb_div_12_neutral = read_csv('Brain_amb_hyb_12_div.csv') %>% 
  filter(adj.P.Val > 0.05)%>% 
  mutate(status = 'Neutral')%>% 
  mutate(Regulated = 'Neutral')
brain_amb_hyb_div_18_neutral = read_csv('Brain_amb_hyb_18_div.csv') %>% 
  filter(adj.P.Val > 0.05)%>% 
  mutate(status = 'Neutral')%>% 
  mutate(Regulated = 'Neutral')


brain_geo_hyb_div_12_neutral = read_csv('Brain_geo_hyb_12_div.csv') %>% 
  filter(adj.P.Val > 0.05)%>% 
  mutate(status = 'Neutral')%>% 
  mutate(Regulated = 'Neutral')
brain_geo_hyb_div_18_neutral = read_csv('Brain_geo_hyb_18_div.csv') %>% 
  filter(adj.P.Val > 0.05)%>% 
  mutate(status = 'Neutral')%>% 
  mutate(Regulated = 'Neutral')


##  brain significantly differentially expressed genes
brain_eco12 = read_csv('Brain_eco_div_12_significant.csv') %>% 
  mutate(status = 'Outlier') %>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
brain_eco18 = read_csv('Brain_eco_div_18_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
brain_plast_amb = read_csv('Brain_ambient_plastic_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
brain_plast_geo = read_csv('Brain_geothermal_plastic_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
brain_plast_hyb = read_csv('Brain_hybrid_plastic_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))

brain_amb_hyb_div_12 = read_csv('Brain_amb_hyb_12_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
brain_amb_hyb_div_18 = read_csv('Brain_amb_hyb_18_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
brain_geo_hyb_div_12 = read_csv('Brain_geo_hyb_12_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
brain_geo_hyb_div_18 = read_csv('Brain_geo_hyb_18_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))

brain_eco12_clean = bind_rows(brain_eco12, 
                              brain_eco12_neutral)

brain_eco18_clean = bind_rows(brain_eco18, 
                              brain_eco18_neutral)

brain_plast_amb_clean = bind_rows(brain_plast_amb, 
                              brain_plast_amb_neutral)

brain_plast_geo_clean = bind_rows(brain_plast_geo, 
                                  brain_plast_geo_neutral)

brain_plast_hyb_clean = bind_rows(brain_plast_hyb, 
                                  brain_plast_hyb_neutral)

brain_amb_hyb_clean12 = bind_rows(brain_amb_hyb_div_12, 
                                brain_amb_hyb_div_12_neutral)
brain_amb_hyb_clean18 = bind_rows(brain_amb_hyb_div_18, 
                                  brain_amb_hyb_div_18_neutral)

brain_geo_hyb_clean12 = bind_rows(brain_geo_hyb_div_12, 
                                  brain_geo_hyb_div_12_neutral)

brain_geo_hyb_clean18 = bind_rows(brain_geo_hyb_div_18, 
                                  brain_geo_hyb_div_18_neutral)



# brain overlap differential exppression ----------------------------------------



inner_join(brain_eco12, 
           brain_eco18, 
           by = 'GeneID') 
## only 2 genes are differentially divergent between ecotypes
## across the two temperatures. 

inner_join(brain_plast_amb, 
           brain_plast_geo, 
           by = 'GeneID') %>% 
  select(GeneID) %>% 
  write_tsv('Brain_plasticity_amb_vs_geo.txt', 
            col_names = F)
## 47 genes are plastic in both the ambient and geothermal ecotypes

inner_join(brain_plast_amb, 
           brain_plast_hyb, 
           by = 'GeneID')%>% 
  select(GeneID) %>% 
  write_tsv('Brain_plasticity_amb_vs_hyb.txt', 
            col_names = F)
## 43 genes are plastic in both the amient and hybridl ecotypes

inner_join(brain_plast_geo, 
           brain_plast_hyb, 
           by = 'GeneID')%>% 
  select(GeneID) %>% 
  write_tsv('Brain_plasticity_geo_vs_hyb.txt', 
            col_names = F)
## 85 genes are plastic in both the geothermal and hybrid ecotypes

inner_join(brain_plast_amb, 
           brain_plast_geo, 
           by = 'GeneID') %>% 
  inner_join(.,
             brain_plast_hyb, 
             by = 'GeneID')%>% 
  select(GeneID) %>% 
  write_tsv('Brain_plasticity_amb_vs_geo_vs_hyb.txt', 
            col_names = F)

## 38 that are plastic between all three ecotypes


# brain quick volcano plots -----------------------------------------------------

volcano_cols = c('#127475',
                 '#dedbd2', 
                 '#f2542d')


# eco divergence vol plots ------------------------------------------------


Brain_eco_div12_volplot = ggplot(data = brain_eco12_clean, 
       aes(x = logFC, 
           y = adj.P.Val))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             aes(col = Regulated))+
  scale_y_reverse()+
  scale_color_manual(values = volcano_cols)+
  labs(x = 'log Fold Change', 
       y = 'Adjusted p-value', 
       title = 'A) Geothermal - ambient divergence 12˚C')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  geom_hline(yintercept = 0.05, 
             col = '#f72585', 
             size = 1, 
             linetype = 'dashed')+
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))
  

Brain_eco_div18_volplot = ggplot(data = brain_eco18_clean, 
       aes(x = logFC, 
           y = adj.P.Val))+
  geom_point(size = 2.5,
             col = 'black')+
  geom_point(size = 2, 
             aes(col = Regulated))+
  scale_y_reverse()+
  scale_color_manual(values = volcano_cols)+
  labs(x = 'log Fold Change', 
       y = 'Adjusted p-value', 
       title = 'B) Geothermal - ambient divergence 18˚C')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  geom_hline(yintercept = 0.05, 
             col = '#f72585', 
             size = 1, 
             linetype = 'dashed')+
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        axis.title.y = element_blank())

Brain_eco_div_vol_plot = Brain_eco_div12_volplot + Brain_eco_div18_volplot

ggsave('EcoDiv_Volcanoe_Plot.svg', 
       plot = Brain_eco_div_vol_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 10)



# overlap plastic gene expression -----------------------------------------

brain_plast_amb = read_csv('Brain_ambient_plastic_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
brain_plast_geo = read_csv('Brain_geothermal_plastic_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
brain_plast_hyb = read_csv('Brain_hybrid_plastic_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))

inner_join(brain_plast_amb, 
           brain_plast_geo, 
           by = 'GeneID') %>% 
  select(GeneID,
         logFC.x, 
         Regulated.x,
         logFC.y, 
         Regulated.y) %>% 
  filter(Regulated.x == 'Up-regulated')

inner_join(brain_plast_amb, 
           brain_plast_geo, 
           by = 'GeneID') %>% 
  inner_join(., 
             brain_plast_hyb, 
             by = 'GeneID') %>%
  select(GeneID, 
         Regulated.x, 
         Regulated.y, 
         Regulated) %>% 
  filter(Regulated == 'Up-regulated')

inner_join(brain_plast_amb, 
           brain_plast_hyb, 
           by = 'GeneID') %>% 
  select(GeneID, 
         Regulated.x, 
         Regulated.y) %>% 
  filter(Regulated.x == 'Up-regulated')

inner_join(brain_plast_geo, 
           brain_plast_hyb, 
           by = 'GeneID') %>% 
  select(GeneID, 
         Regulated.x, 
         Regulated.y) %>% 
  filter(Regulated.x == 'Up-regulated')

# plasticity vol plots ----------------------------------------------------


Brain_amb_plast_volplot = ggplot(data = brain_plast_amb_clean, 
       aes(x = logFC, 
           y = adj.P.Val))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             aes(col = Regulated))+
  scale_y_reverse()+
  scale_color_manual(values = volcano_cols)+
  labs(x = 'log Fold Change', 
       y = 'Adjusted p-value', 
       title = 'A) Ambient plasticity')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  geom_hline(yintercept = 0.05, 
             col = '#f72585', 
             size = 1, 
             linetype = 'dashed')+
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

Brain_geo_plast_volplot = ggplot(data = brain_plast_geo_clean, 
       aes(x = logFC, 
           y = adj.P.Val))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             aes(col = Regulated))+
  scale_y_reverse()+
  scale_color_manual(values = volcano_cols)+
  labs(x = 'log Fold Change', 
       y = 'Adjusted p-value', 
       title = 'B) Geothermal plasticity')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  geom_hline(yintercept = 0.05, 
             col = '#f72585', 
             size = 1, 
             linetype = 'dashed')+
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        axis.title.y = element_blank())


Brain_hyb_plast_volplot = ggplot(data = brain_plast_hyb_clean, 
       aes(x = logFC, 
           y = adj.P.Val))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             aes(col = Regulated))+
  scale_y_reverse()+
  scale_color_manual(values = volcano_cols)+
  labs(x = 'log Fold Change', 
       y = 'Adjusted p-value', 
       title = 'C) Hybrid plasticity')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  geom_hline(yintercept = 0.05, 
             col = '#f72585', 
             size = 1, 
             linetype = 'dashed')+
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        axis.title.y = element_blank())
  
plasticity_volplot = Brain_amb_plast_volplot+Brain_geo_plast_volplot+Brain_hyb_plast_volplot

ggsave('Plasticity_Volcanoe_Plot.svg', 
       plot = plasticity_volplot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 40, 
       height = 10)

# hybrid divergence vol plots ---------------------------------------------

Brain_amb_hyb_12_volplot = ggplot(data = brain_amb_hyb_clean12, 
                                 aes(x = logFC, 
                                     y = adj.P.Val))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             aes(col = Regulated))+
  scale_y_reverse()+
  scale_color_manual(values = volcano_cols)+
  labs(x = 'log Fold Change', 
       y = 'Adjusted p-value', 
       title = 'A) Ambient - Hybrid divergence 12˚C')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  geom_hline(yintercept = 0.05, 
             col = '#f72585', 
             size = 1, 
             linetype = 'dashed')+
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

Brain_amb_hyb_18_volplot = ggplot(data = brain_amb_hyb_clean18, 
                                  aes(x = logFC, 
                                      y = adj.P.Val))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             aes(col = Regulated))+
  scale_y_reverse()+
  scale_color_manual(values = volcano_cols)+
  labs(x = 'log Fold Change', 
       y = 'Adjusted p-value', 
       title = 'B) Ambient - Hybrid divergence 18˚C')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  geom_hline(yintercept = 0.05, 
             col = '#f72585', 
             size = 1, 
             linetype = 'dashed')+
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        axis.title.y = element_blank())


Brain_geo_hyb_12_volplot = ggplot(data = brain_geo_hyb_clean12, 
                                  aes(x = logFC, 
                                      y = adj.P.Val))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             aes(col = Regulated))+
  scale_y_reverse()+
  scale_color_manual(values = volcano_cols)+
  labs(x = 'log Fold Change', 
       y = 'Adjusted p-value', 
       title = 'C) Geothermal - Hybrid divergence 12˚C')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  geom_hline(yintercept = 0.05, 
             col = '#f72585', 
             size = 1, 
             linetype = 'dashed')+
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

Brain_geo_hyb_18_volplot = ggplot(data = brain_geo_hyb_clean18, 
                                  aes(x = logFC, 
                                      y = adj.P.Val))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             aes(col = Regulated))+
  scale_y_reverse()+
  scale_color_manual(values = volcano_cols)+
  labs(x = 'log Fold Change', 
       y = 'Adjusted p-value', 
       title = 'D) Geothermal - Hybrid divergence 18˚C')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  geom_hline(yintercept = 0.05, 
             col = '#f72585', 
             size = 1, 
             linetype = 'dashed')+
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        axis.title.y = element_blank())

hybrid_div_volplots = Brain_amb_hyb_12_volplot+Brain_amb_hyb_18_volplot+Brain_geo_hyb_12_volplot+Brain_geo_hyb_18_volplot

ggsave('Hybrid_divergence_Volcanoe_Plot.svg', 
       plot = hybrid_div_volplots, 
       dpi = 'retina', 
       units = 'cm', 
       width = 40, 
       height = 20)



# brain hybrid vs pure divergence and plasticity ----------------------------------------------

neutral_amb_hyb_12 = read_csv('Brain_amb_hyb_12_div.csv') %>% 
  filter(adj.P.Val > 0.05) %>% 
  mutate(status = 'Neutral')
significant_amb_hyb_12 = read_csv('Brain_amb_hyb_12_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')

neutral_amb_hyb_18 = read_csv('Brain_amb_hyb_18_div.csv') %>% 
  filter(adj.P.Val > 0.05) %>% 
  mutate(status = 'Neutral')
significant_amb_hyb_18 = read_csv('Brain_amb_hyb_18_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')

neutral_geo_hyb_12 = read_csv('Brain_geo_hyb_12_div.csv') %>% 
  filter(adj.P.Val > 0.05) %>% 
  mutate(status = 'Neutral')
significant_geo_hyb_12 = read_csv('Brain_geo_hyb_12_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')

neutral_geo_hyb_18 = read_csv('Brain_geo_hyb_18_div.csv') %>% 
  filter(adj.P.Val > 0.05) %>% 
  mutate(status = 'Neutral')
significant_geo_hyb_18 = read_csv('Brain_geo_hyb_18_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')

brain_amb_hyb_12_clean = bind_rows(significant_amb_hyb_12, 
                                   neutral_amb_hyb_12)

brain_amb_hyb_18_clean = bind_rows(significant_amb_hyb_18, 
                                   neutral_amb_hyb_18)

brain_geo_hyb_12_clean = bind_rows(significant_geo_hyb_12, 
                                   neutral_geo_hyb_12)
brain_geo_hyb_18_clean = bind_rows(significant_geo_hyb_18, 
                                   neutral_geo_hyb_18)


ggplot(data = brain_amb_hyb_12_clean, 
       aes(x = logFC, 
           y = adj.P.Val, 
           col = status))+
  geom_point()+
  scale_y_reverse()

ggplot(data = brain_amb_hyb_18_clean, 
       aes(x = logFC, 
           y = adj.P.Val, 
           col = status))+
  geom_point()+
  scale_y_reverse()

ggplot(data = brain_geo_hyb_12_clean, 
       aes(x = logFC, 
           y = adj.P.Val, 
           col = status))+
  geom_point()+
  scale_y_reverse()

ggplot(data = brain_geo_hyb_18_clean, 
       aes(x = logFC, 
           y = adj.P.Val, 
           col = status))+
  geom_point()+
  scale_y_reverse()


ggplot(data = brain_plast_hyb_clean, 
       aes(x = logFC, 
           y = adj.P.Val, 
           col = status))+
  geom_point()+
  scale_y_reverse()




## quantify overlap in divergence at same temperature
div_pure_hyb_12 = inner_join(significant_amb_hyb_12, 
           significant_geo_hyb_12, 
           by = 'GeneID')

ggplot(data = div_pure_hyb_12, 
       aes(x = logFC.y, 
           y = logFC.x))+
  geom_point()


div_pure_hyb_18 = inner_join(significant_amb_hyb_18, 
           significant_geo_hyb_18, 
           by = 'GeneID')
ggplot(data = div_pure_hyb_18, 
       aes(x = logFC.y, 
           y = logFC.x))+
  geom_point()


## plasticity in divergence
amb_hyb_div_common = inner_join(significant_amb_hyb_12, 
           significant_amb_hyb_18, 
           by = 'GeneID') 

# amb_hyb_12_only = anti_join(significant_amb_hyb_12, 
#           amb_hyb_div_common, 
#           by = 'GeneID')
# 
# amb_hyb_18_only = anti_join(significant_amb_hyb_18, 
#                             amb_hyb_div_common, 
#                             by = 'GeneID')


geo_hyb_div_common = inner_join(significant_geo_hyb_12, 
           significant_geo_hyb_18, 
           by = 'GeneID')

# geo_hyb_12_only = anti_join(significant_geo_hyb_12, 
#                             geo_hyb_div_common, 
#                             by = 'GeneID')
# 
# geo_hyb_18_only = anti_join(significant_geo_hyb_18, 
#                             geo_hyb_div_common, 
#                             by = 'GeneID')



Div_eco_temp = inner_join(amb_hyb_div_common, 
           geo_hyb_div_common, 
           by = 'GeneID')

amb_hyb_div_common_temp = anti_join(amb_hyb_div_common,
          Div_eco_temp,
          by = 'GeneID')

geo_hyb_div_common_temp = anti_join(geo_hyb_div_common,
                                    Div_eco_temp,
                                    by = 'GeneID')

amb_hyb_12_only = anti_join(significant_amb_hyb_12, 
                            Div_eco_temp, 
                            by = 'GeneID')
amb_hyb_18_only = anti_join(significant_amb_hyb_18, 
                            Div_eco_temp, 
                            by = 'GeneID')


geo_hyb_12_only = anti_join(significant_geo_hyb_12, 
                            Div_eco_temp, 
                            by = 'GeneID')
geo_hyb_18_only = anti_join(significant_geo_hyb_18, 
                            Div_eco_temp, 
                            by = 'GeneID')
# Liver temp plast --------------------------------------------------------
Liver_limma_all = read_csv('Liver_LIMMA_model_results_all.csv')

## All genes 

Liver_eco12_neutral = read_csv('Liver_eco_div_12.csv') %>% 
  filter(adj.P.Val > 0.05) %>% 
  mutate(status = 'Neutral') %>% 
  mutate(Regulated = 'Neutral')
Liver_eco18_neutral = read_csv('Liver_eco_div_18.csv') %>% 
  filter(adj.P.Val > 0.05)%>% 
  mutate(status = 'Neutral') %>% 
  mutate(Regulated = 'Neutral')
Liver_plast_amb_neutral = read_csv('Liver_ambient_plastic.csv') %>% 
  filter(adj.P.Val > 0.05)%>% 
  mutate(status = 'Neutral') %>% 
  mutate(Regulated = 'Neutral')
Liver_plast_geo_neutral = read_csv('Liver_geothermal_plastic.csv') %>% 
  filter(adj.P.Val > 0.05)%>% 
  mutate(status = 'Neutral') %>% 
  mutate(Regulated = 'Neutral')
Liver_plast_hyb_neutral = read_csv('Liver_hybrid_plastic.csv') %>% 
  filter(adj.P.Val > 0.05)%>% 
  mutate(status = 'Neutral') %>% 
  mutate(Regulated = 'Neutral')

##  Liver significantly differentially expressed genes
Liver_eco12 = read_csv('Liver_eco_div_12.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
Liver_eco18 = read_csv('Liver_eco_div_18.csv')%>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
Liver_plast_amb = read_csv('Liver_ambient_plastic.csv')%>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
Liver_plast_geo = read_csv('Liver_geothermal_plastic.csv')%>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))
Liver_plast_hyb = read_csv('Liver_hybrid_plastic.csv')%>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))


Liver_eco12_clean = bind_rows(Liver_eco12, 
                              Liver_eco12_neutral)

Liver_eco18_clean = bind_rows(Liver_eco18, 
                              Liver_eco18_neutral)

Liver_plast_amb_clean = bind_rows(Liver_plast_amb, 
                                  Liver_plast_amb_neutral)

Liver_plast_geo_clean = bind_rows(Liver_plast_geo, 
                                  Liver_plast_geo_neutral)

Liver_plast_hyb_clean = bind_rows(Liver_plast_hyb, 
                                  Liver_plast_hyb_neutral)


# Liver overlap differential exppression ----------------------------------------

inner_join(Liver_eco12, 
           Liver_eco18, 
           by = 'GeneID') 

inner_join(Liver_plast_amb, 
           Liver_plast_geo, 
           by = 'GeneID') %>% 
  select(GeneID) %>% 
  write_tsv('Liver_plasticity_amb_vs_geo.txt', 
            col_names = F)
## 47 genes are plastic in both the ambient and geothermal ecotypes

inner_join(Liver_plast_amb, 
           Liver_plast_hyb, 
           by = 'GeneID')%>% 
  select(GeneID) %>% 
  write_tsv('Liver_plasticity_amb_vs_hyb.txt', 
            col_names = F)
## 43 genes are plastic in both the amient and hybridl ecotypes

inner_join(Liver_plast_geo, 
           Liver_plast_hyb, 
           by = 'GeneID')%>% 
  select(GeneID) %>% 
  write_tsv('Liver_plasticity_geo_vs_hyb.txt', 
            col_names = F)
## 85 genes are plastic in both the geothermal and hybrid ecotypes

inner_join(Liver_plast_amb, 
           Liver_plast_geo, 
           by = 'GeneID') %>% 
  inner_join(.,
             Liver_plast_hyb, 
             by = 'GeneID')%>% 
  select(GeneID) %>% 
  write_tsv('Liver_plasticity_amb_vs_geo_vs_hyb.txt', 
            col_names = F)

## 38 that are plastic between all three ecotypes

# Liver quick volcano plots -----------------------------------------------------
volcano_cols = c('#127475',
                 '#dedbd2', 
                 '#f2542d')
vol_plots_no_DEG = c('#dedbd2')

# liver eco div vol plots -------------------------------------------------


liver_eco_div_12_volplot = ggplot(data = Liver_eco12_clean, 
       aes(x = logFC, 
           y = adj.P.Val))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             aes(col = Regulated))+
  scale_y_reverse()+
  scale_color_manual(values = vol_plots_no_DEG)+
  labs(x = 'log Fold Change', 
       y = 'Adjusted p-value', 
       title = 'A) Geothermal - ambient divergence 12˚C')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  geom_hline(yintercept = 0.05, 
             col = '#f72585', 
             size = 1, 
             linetype = 'dashed')+
  theme(
    legend.position = 'none',
        panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

liver_eco_div_18_volplot = ggplot(data = Liver_eco18_clean, 
       aes(x = logFC, 
           y = adj.P.Val))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             aes(col = Regulated))+
  scale_y_reverse()+
  scale_color_manual(values = vol_plots_no_DEG)+
  labs(x = 'log Fold Change', 
       y = 'Adjusted p-value', 
       title = 'B) Geothermal - ambient divergence 18˚C')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  geom_hline(yintercept = 0.05, 
             col = '#f72585', 
             size = 1, 
             linetype = 'dashed')+
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        axis.title.y = element_blank())


liver_eco_div_volplot = liver_eco_div_12_volplot+liver_eco_div_18_volplot

ggsave('Liver_EcoDiv_Volcanoe_Plot.svg', 
       plot = liver_eco_div_volplot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 10)

# liver plasticity vol plots ----------------------------------------------


liver_amb_plast_volplot = ggplot(data = Liver_plast_amb_clean, 
       aes(x = logFC, 
           y = adj.P.Val))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             aes(col = Regulated))+
  scale_y_reverse()+
  scale_color_manual(values = volcano_cols)+
  labs(x = 'log Fold Change', 
       y = 'Adjusted p-value', 
       title = 'A) Ambient plasticity')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  geom_hline(yintercept = 0.05, 
             col = '#f72585', 
             size = 1, 
             linetype = 'dashed')+
  theme(
    legend.position = 'none',
    panel.grid = element_blank(), 
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12))


liver_geo_plast_volplot = ggplot(data = Liver_plast_geo_clean, 
       aes(x = logFC, 
           y = adj.P.Val))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             aes(col = Regulated))+
  scale_y_reverse()+
  scale_color_manual(values = volcano_cols)+
  labs(x = 'log Fold Change', 
       y = 'Adjusted p-value', 
       title = 'B) Geothermal plasticity')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  geom_hline(yintercept = 0.05, 
             col = '#f72585', 
             size = 1, 
             linetype = 'dashed')+
  theme(
    legend.position = 'none',
    panel.grid = element_blank(), 
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12), 
    axis.title.y = element_blank())


liver_hyb_plast_volplot = ggplot(data = Liver_plast_hyb_clean, 
       aes(x = logFC, 
           y = adj.P.Val))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             aes(col = Regulated))+
  scale_y_reverse()+
  scale_color_manual(values = volcano_cols)+
  labs(x = 'log Fold Change', 
       y = 'Adjusted p-value', 
       title = 'C) Hybrid plasticity')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  geom_hline(yintercept = 0.05, 
             col = '#f72585', 
             size = 1, 
             linetype = 'dashed')+
  theme(
    legend.position = 'none',
    panel.grid = element_blank(), 
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12), 
    axis.title.y = element_blank())

liver_plasticity_volplots = liver_amb_plast_volplot+liver_geo_plast_volplot+liver_hyb_plast_volplot

ggsave('Liver_Plasticity_Volcanoe_Plot.svg', 
       plot = liver_plasticity_volplots, 
       dpi = 'retina', 
       units = 'cm', 
       width = 40, 
       height = 10)

#liver hybrid vs pure divergence and plasticity ----------------------------------------------
  
liv_neutral_amb_hyb_12 = read_csv('liver_amb_hyb_12_div.csv') %>% 
  filter(adj.P.Val > 0.05) %>% 
  mutate(status = 'Neutral')%>% 
  mutate(Regulated = 'Neutral')

liv_significant_amb_hyb_12 = read_csv('liver_amb_hyb_12_div.csv') %>% 
  filter(adj.P.Val <= 0.05)  %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))

liv_neutral_amb_hyb_18 = read_csv('liver_amb_hyb_18_div.csv') %>% 
  filter(adj.P.Val > 0.05) %>% 
  mutate(status = 'Neutral')%>% 
  mutate(Regulated = 'Neutral')
liv_significant_amb_hyb_18 = read_csv('liver_amb_hyb_18_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))

liv_neutral_geo_hyb_12 = read_csv('liver_geo_hyb_12_div.csv') %>% 
  filter(adj.P.Val > 0.05) %>% 
  mutate(status = 'Neutral')%>% 
  mutate(Regulated = 'Neutral')
liv_significant_geo_hyb_12 = read_csv('liver_geo_hyb_12_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))

liv_neutral_geo_hyb_18 = read_csv('liver_geo_hyb_18_div.csv') %>% 
  filter(adj.P.Val > 0.05) %>% 
  mutate(status = 'Neutral')%>% 
  mutate(Regulated = 'Neutral')
liv_significant_geo_hyb_18 = read_csv('liver_geo_hyb_18_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier')%>% 
  mutate(Regulated = case_when(
    logFC >=0 ~ "Up-regulated",
    logFC <= 0 ~ "Down-regulated"
  ))

liver_amb_hyb_12_clean = bind_rows(liv_significant_amb_hyb_12, 
                                   liv_neutral_amb_hyb_12)

liver_amb_hyb_18_clean = bind_rows(liv_significant_amb_hyb_18,
                                   liv_neutral_amb_hyb_18)

liver_geo_hyb_12_clean = bind_rows(liv_significant_geo_hyb_12, 
                                   liv_neutral_geo_hyb_12)
liver_geo_hyb_18_clean = bind_rows(liv_significant_geo_hyb_18,
                                   liv_neutral_geo_hyb_18)


liver_amb_hyb_div_12_volplot = ggplot(data = liver_amb_hyb_12_clean, 
       aes(x = logFC, 
           y = adj.P.Val))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             aes(col = Regulated))+
  scale_y_reverse()+
  scale_color_manual(values = volcano_cols)+
  labs(x = 'log Fold Change', 
       y = 'Adjusted p-value', 
       title = 'A) Ambient - Hybrid divergence 12˚C')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  geom_hline(yintercept = 0.05, 
             col = '#f72585', 
             size = 1, 
             linetype = 'dashed')+
  theme(
    legend.position = 'none',
    panel.grid = element_blank(), 
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12), 
    # axis.title.y = element_blank()
    )

liver_amb_hyb_div_18_volplot = ggplot(data = liver_amb_hyb_18_clean,
       aes(x = logFC,
           y = adj.P.Val))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             aes(col = Regulated))+
  scale_y_reverse()+
  scale_color_manual(values = vol_plots_no_DEG)+
  labs(x = 'log Fold Change', 
       y = 'Adjusted p-value', 
       title = 'B) Ambient - Hybrid divergence 18˚C')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  geom_hline(yintercept = 0.05, 
             col = '#f72585', 
             size = 1, 
             linetype = 'dashed')+
  theme(
    legend.position = 'none',
    panel.grid = element_blank(), 
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12), 
    axis.title.y = element_blank()
  )

liver_geo_hyb_div_12_volplot = ggplot(data = liver_geo_hyb_12_clean, 
       aes(x = logFC, 
           y = adj.P.Val))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             aes(col = Regulated))+
  scale_y_reverse()+
  scale_color_manual(values = volcano_cols)+
  labs(x = 'log Fold Change', 
       y = 'Adjusted p-value', 
       title = 'C) Geothermal - Hybrid divergence 12˚C')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  geom_hline(yintercept = 0.05, 
             col = '#f72585', 
             size = 1, 
             linetype = 'dashed')+
  theme(
    legend.position = 'none',
    panel.grid = element_blank(), 
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12), 
    # axis.title.y = element_blank()
  )


liver_geo_hyb_div_18_volplot = ggplot(data = liver_geo_hyb_18_clean,
       aes(x = logFC,
           y = adj.P.Val))+
  geom_point(size = 2.5, 
             col = 'black')+
  geom_point(size = 2, 
             aes(col = Regulated))+
  scale_y_reverse()+
  scale_color_manual(values = vol_plots_no_DEG)+
  labs(x = 'log Fold Change', 
       y = 'Adjusted p-value', 
       title = 'D) Geothermal - Hybrid divergence 18˚C')+
  geom_vline(xintercept = 0, 
             col = 'black', 
             size = 1, 
             linetype = 'dashed')+
  geom_hline(yintercept = 0.05, 
             col = '#f72585', 
             size = 1, 
             linetype = 'dashed')+
  theme(
    legend.position = 'none',
    panel.grid = element_blank(), 
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12), 
    axis.title.y = element_blank()
  )

liver_hyb_div_volplot = (liver_amb_hyb_div_12_volplot+liver_amb_hyb_div_18_volplot)/(liver_geo_hyb_div_12_volplot+liver_geo_hyb_div_18_volplot)


ggsave('Liver_Hybrid_divergence_Volcanoe_Plot.tiff', 
       plot = liver_hyb_div_volplot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 40, 
       height = 20)


## quantify overlap in divergence at same temperature
liv_div_pure_hyb_12 = inner_join(liv_significant_amb_hyb_12, 
                             liv_significant_geo_hyb_12, 
                             by = 'GeneID')




# amb_hyb_12_only = anti_join(significant_amb_hyb_12, 
#           amb_hyb_div_common, 
#           by = 'GeneID')
# 
# amb_hyb_18_only = anti_join(significant_amb_hyb_18, 
#                             amb_hyb_div_common, 
#                             by = 'GeneID')


# geo_hyb_div_common = inner_join(significant_geo_hyb_12, 
#                                 significant_geo_hyb_18, 
#                                 by = 'GeneID')

# geo_hyb_12_only = anti_join(significant_geo_hyb_12, 
#                             geo_hyb_div_common, 
#                             by = 'GeneID')
# 
# geo_hyb_18_only = anti_join(significant_geo_hyb_18, 
#                             geo_hyb_div_common, 
#                             by = 'GeneID')



# Div_eco_temp = inner_join(amb_hyb_div_common, 
#                           geo_hyb_div_common, 
#                           by = 'GeneID')
# 
# amb_hyb_div_common_temp = anti_join(amb_hyb_div_common,
#                                     Div_eco_temp,
#                                     by = 'GeneID')
# 
# geo_hyb_div_common_temp = anti_join(geo_hyb_div_common,
#                                     Div_eco_temp,
#                                     by = 'GeneID')
# 
# amb_hyb_12_only = anti_join(significant_amb_hyb_12, 
#                             Div_eco_temp, 
#                             by = 'GeneID')
# amb_hyb_18_only = anti_join(significant_amb_hyb_18, 
#                             Div_eco_temp, 
#                             by = 'GeneID')
# 
# 
# geo_hyb_12_only = anti_join(significant_geo_hyb_12, 
#                             Div_eco_temp, 
#                             by = 'GeneID')
# geo_hyb_18_only = anti_join(significant_geo_hyb_18, 
#                             Div_eco_temp, 
#                             by = 'GeneID')

# stickleback v5 annotation -----------------------------------------------

gene_annotation = read_tsv('~/Parsons_Postdoc/Stickleback_Genomic/Stickleback_Annotation_features/stickleback_v5_ensembl_genes.gff3.gz', 
                           col_names = F, 
                           skip = 1) %>% 
  # filter(X3 %in% c('gene', 
  #                  'exon', 
  #                  'CDS')) %>% 
  group_by(X1) %>% 
  arrange(X4, 
          X5) %>% 
  ## arrange each gene by its start and end points on each chromosome
  mutate(mid = X4 + (X5-X4)/2) %>% 
  dplyr::select(X1, 
                X3:X5, 
                X9:mid) %>% 
  rename(chromosome = X1, 
         feature = X3, 
         start = X4, 
         end = X5, 
         gene_id = X9, 
         position = mid) %>% 
  na.omit() %>% 
  separate(col = gene_id, 
           into = c('ensemble_id', 
                    'gene_name', 
                    'parent_code', 
                    'gene_name2'), 
           sep = ';') %>%
  separate(col = gene_name, 
           into = c('Garbage', 
                    'gene_name'), 
           sep = '=') %>% 
  dplyr::select(chromosome, 
                position, 
                start, ,
                end, 
                feature,
                gene_name) %>% 
  na.omit()

gene_annotation %>% 
  select(chromosome, 
         gene_id) %>% 
  pull(gene_id)


gene_annotation %>% 
  ungroup() %>% 
  select(gene_name) %>% 
  write_tsv('exonID.txt')



