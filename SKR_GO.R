################################################
## SKR gene ontology graphs
##
## Matthew Brachmann (MKB) @phdmattyb
##
## 17.11.2025
###############################################

library(tidyverse)
library(patchwork)

theme_set(theme_bw())

setwd('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/')

db_col_pal = c(
  # '#031d44', 
  '#04395e', 
  '#70a288',
  '#dab785', 
  '#d5896f')
# ambient plastic go ------------------------------------------------------

amb_plast_go = read_csv('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/Ambient_plastic_unique_GO/Enrichment_GO/GO_AllLists.csv') %>% 
  dplyr::select(GO, 
                `Z-score`, 
                LogP, 
                Enrichment, 
                `Log(q-value)`, 
                Category, 
                Description) %>% 
  unite(GO_cat, 
        c('GO', 
          'Description'), 
        sep = ' - ', 
        remove = F)

amb_plast_go_plot = amb_plast_go %>% 
  # mutate_if(is.character, toupper) %>% 
  arrange(`Z-score`) %>% 
  # mutate(Category = factor(Category,levels = Category)) %>%
  mutate(GO_cat = factor(GO_cat, unique(GO_cat))) %>% 
  # rename(Database = database) %>% 
  ggplot(aes(y = GO_cat,
             x = `Z-score`, 
             fill = Category))+
  # ggplot(aes(y = Term, 
  #            x = log_adj_pval, 
  #            col = Combined.Score, 
  #            fill = Combined.Score))+
  geom_col()+
  # scale_x_reverse()+
  # scale_y_reverse()+
  scale_fill_manual(values = db_col_pal)+
  labs(x = 'Log q-value', 
       title = 'D)')+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_blank(), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')

ggsave('Metascape_GO_Ambient_plastic_plots.svg',
       plot = amb_plast_go_plot,
       dpi = 'retina',
       units = 'cm',
       width = 15,
       height = 5)


# geothermal plastic go ---------------------------------------------------

geo_plast_go = read_csv('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/Geothermal_plastic_unique_GO/Enrichment_GO/GO_AllLists.csv') %>% 
  dplyr::select(GO, 
                `Z-score`, 
                LogP, 
                Enrichment, 
                `Log(q-value)`, 
                Category, 
                Description) %>% 
  unite(GO_cat, 
        c('GO', 
          'Description'), 
        sep = ' - ', 
        remove = F)

geo_plast_go_plot = geo_plast_go %>% 
  # mutate_if(is.character, toupper) %>% 
  arrange(`Z-score`) %>% 
  # mutate(Category = factor(Category,levels = Category)) %>%
  mutate(GO_cat = factor(GO_cat, unique(GO_cat))) %>% 
  # rename(Database = database) %>% 
  ggplot(aes(y = GO_cat,
             x = `Z-score`, 
             fill = Category))+
  # ggplot(aes(y = Term, 
  #            x = log_adj_pval, 
  #            col = Combined.Score, 
  #            fill = Combined.Score))+
  geom_col()+
  # scale_x_reverse()+
  # scale_y_reverse()+
  scale_fill_manual(values = db_col_pal)+
  labs(x = 'Log q-value', 
       title = 'E)')+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_blank(), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')

ggsave('Metascape_GO_Geothermal_plastic_plots.svg',
       plot = geo_plast_go_plot,
       dpi = 'retina',
       units = 'cm',
       width = 30,
       height = 10)

# hybrid plastic go -------------------------------------------------------

hyb_plast_go = read_csv('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/Hybrid_plastic_unique_GO/Enrichment_GO/GO_AllLists.csv') %>% 
  dplyr::select(GO, 
                `Z-score`, 
                LogP, 
                Enrichment, 
                `Log(q-value)`, 
                Category, 
                Description) %>% 
  unite(GO_cat, 
        c('GO', 
          'Description'), 
        sep = ' - ', 
        remove = F)

hyb_plast_go_plot = hyb_plast_go %>% 
  # mutate_if(is.character, toupper) %>% 
  arrange(`Z-score`) %>% 
  # mutate(Category = factor(Category,levels = Category)) %>%
  mutate(GO_cat = factor(GO_cat, unique(GO_cat))) %>% 
  # rename(Database = database) %>% 
  ggplot(aes(y = GO_cat,
             x = `Z-score`, 
             fill = Category))+
  # ggplot(aes(y = Term, 
  #            x = log_adj_pval, 
  #            col = Combined.Score, 
  #            fill = Combined.Score))+
  geom_col()+
  # scale_x_reverse()+
  # scale_y_reverse()+
  scale_fill_manual(values = db_col_pal)+
  labs(x = 'Log q-value', 
       title = 'F)')+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_blank(), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')

ggsave('Metascape_GO_hybrid_plastic_plots.svg',
       plot = hyb_plast_go_plot,
       dpi = 'retina',
       units = 'cm',
       width = 40,
       height = 15)

# overlap plastic go ------------------------------------------------------
overlap_plast_go = read_csv('~/Parsons_Postdoc/SKR_hybrid_Gene_expression/overlap_plastic_responses_GO/Enrichment_GO/GO_AllLists.csv') %>% 
  dplyr::select(GO, 
                `Z-score`, 
                LogP, 
                Enrichment, 
                `Log(q-value)`, 
                Category, 
                Description) %>% 
  unite(GO_cat, 
        c('GO', 
          'Description'), 
        sep = ' - ', 
        remove = F)

overlap_plast_go_plot = overlap_plast_go %>% 
  # mutate_if(is.character, toupper) %>% 
  arrange(`Z-score`) %>% 
  # mutate(Category = factor(Category,levels = Category)) %>%
  mutate(GO_cat = factor(GO_cat, unique(GO_cat))) %>% 
  # rename(Database = database) %>% 
  ggplot(aes(y = GO_cat,
             x = `Z-score`, 
             fill = Category))+
  # ggplot(aes(y = Term, 
  #            x = log_adj_pval, 
  #            col = Combined.Score, 
  #            fill = Combined.Score))+
  geom_col()+
  # scale_x_reverse()+
  # scale_y_reverse()+
  scale_fill_manual(values = db_col_pal)+
  labs(x = 'Log q-value', 
       title = 'G)')+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_blank(), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')

ggsave('Metascape_GO_Overlap_plastic_plots.svg',
       plot = overlap_plast_go_plot,
       dpi = 'retina',
       units = 'cm',
       width = 30,
       height = 10)


# combine plastic plots ---------------------------------------------------

plastic_combo_plots = (amb_plast_go_plot+geo_plast_go_plot+hyb_plast_go_plot)/overlap_plast_go_plot

ggsave('TEST_Metascape_GO_plastic_plots.svg',
       plot = plastic_combo_plots,
       dpi = 'retina',
       units = 'cm',
       width = 100,
       height = 30)

