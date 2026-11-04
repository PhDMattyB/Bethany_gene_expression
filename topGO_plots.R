### topGO SKR gene expression graphs

setwd('C:/Users/phdma/OneDrive/OneDrive - University of Glasgow/Documents/Parsons_Postdoc/SKR_Hybrid_Gene_expression')

library(tidyverse)

go_pal = c('#84a59d', 
           '#f28482', 
           '#f6bd60')

go_pal2 = c('#184e77', 
            '#52b69a', 
            '#b5e48c')

go_pal3 = c('#5f0f40',
            '#184e77', 
            '#52b69a')

go_pal4 = c('#d1495b',
            '#00798c', 
            '#edae49')

# amb geo div 12 ----------------------------------------------------------
div_12_GO_BP = read_csv('Div_amb_geo_12_topGO_BP.csv') %>% 
  mutate(Type = 'Biological process')
  # write_csv('Div_amb_geo_12_topGO_BP.csv')
div_12_GO_MF = read_csv("Div_amb_geo_12_topGO_MF.csv") %>% 
  mutate(Type = 'Molecular Function')
  # write_csv("Div_amb_geo_12_topGO_MF.csv")
div_12_GO_CC = read_csv("Div_amb_geo_12_topGO_CC.csv") %>% 
  mutate(Type = 'Cellular Component')
# %>% 
#   write_csv("Div_amb_geo_12_topGO_CC.csv")

div_12_topGO = bind_rows(div_12_GO_BP, 
                         div_12_GO_MF, 
                         div_12_GO_CC) %>% 
  unite(col = 'GO_Term', 
        c('GO.ID', 
          'Term'), 
        sep = ': ') %>% 
  mutate(neg_log_pval = -log10(weight01Fisher))

div_12_topGO_plot = ggplot(data = div_12_topGO, 
       aes(y = GO_Term, 
           x = neg_log_pval))+
  geom_col(aes(fill = Type), 
           col = 'black')+
  labs(x = '-Log(p-value)')+
  scale_fill_manual(values = go_pal4)+
  theme(axis.title.y = element_blank(), 
        panel.grid.major.y = element_blank(),
        legend.position = 'none')

ggsave('Div_12_amb_geo_topGO_plot.svg', 
       plot = div_12_topGO_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20, 
       height = 10)


# amb geo div 18 ----------------------------------------------------------

div_18_GO_BP %>% write_csv('Div_amb_geo_18_topGO_BP.csv')
div_18_GO_MF %>% write_csv("Div_amb_geo_18_topGO_MF.csv")
div_18_GO_CC %>% write_csv("Div_amb_geo_18_topGO_CC.csv")

div_18_GO_BP = read_csv('Div_amb_geo_18_topGO_BP.csv') %>% 
  mutate(Type = 'Biological process')
div_18_GO_MF = read_csv("Div_amb_geo_18_topGO_MF.csv") %>% 
  mutate(Type = 'Molecular Function')
# div_18_GO_CC = read_csv("Div_amb_geo_18_topGO_CC.csv") %>% 
#   mutate(Type = 'Cellular Component')

div_18_topGO = bind_rows(div_18_GO_BP, 
                         div_18_GO_MF) %>% 
  unite(col = 'GO_Term', 
        c('GO.ID', 
          'Term'), 
        sep = ': ') %>% 
  mutate(neg_log_pval = -log10(weight01Fisher))

div_18_pal = c('#d1495b',
            '#edae49')

div_18_topGO_plot = ggplot(data = div_18_topGO, 
                           aes(y = GO_Term, 
                               x = neg_log_pval))+
  geom_col(aes(fill = Type), 
           col = 'black')+
  labs(x = '-Log(p-value)')+
  scale_fill_manual(values = div_18_pal)+
  theme(axis.title.y = element_blank(), 
        panel.grid.major.y = element_blank(),
        legend.position = 'none')


ggsave('Div_18_amb_geo_topGO_plot.svg', 
       plot = div_18_topGO_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20, 
       height = 10)



# Ambient plasticity ------------------------------------------------------
amb_plast_GO_BP %>% write_csv('Amb_plast_topGO_BP.csv')
amb_plast_GO_MF %>% write_csv('Amb_plast_topGO_MF.csv')
amb_plast_GO_CC %>% write_csv('Amb_plast_topGO_CC.csv')


amb_plast_GO_BP = read_csv('Amb_plast_topGO_BP.csv') %>% 
  mutate(Type = 'Biological Process')
amb_plast_GO_MF = read_csv("Amb_plast_topGO_MF.csv") %>% 
  mutate(Type = 'Molecular Function')
amb_plast_GO_CC = read_csv("Amb_plast_topGO_CC.csv") %>%
  mutate(Type = 'Cellular Component')

amb_plast_topGO = bind_rows(amb_plast_GO_BP, 
                         amb_plast_GO_MF, 
                         amb_plast_GO_CC) %>% 
  unite(col = 'GO_Term', 
        c('GO.ID', 
          'Term'), 
        sep = ': ') %>% 
  mutate(neg_log_pval = -log10(weight01Fisher), 
         Type = factor(Type,
           levels = c("Biological Process", 
                      "Molecular Function", 
                      "Cellular Component")))%>% 
  arrange(desc(Type), neg_log_pval) %>% 
  mutate(GO_Term = factor(GO_Term,
                          levels = unique(GO_Term))) 



amb_plast_topGO_plot = ggplot(data = amb_plast_topGO, 
                           aes(y = GO_Term, 
                               x = neg_log_pval, 
                               group = Type))+
  geom_col(aes(fill = Type),
           col = 'black')+
  labs(x = '-Log(p-value)')+
  scale_fill_manual(values = go_pal4)+
  theme(axis.title.y = element_blank(), 
        panel.grid.major.y = element_blank(), 
        legend.position = 'none')


ggsave('Amb_plast_topGO_plot.svg', 
       plot = amb_plast_topGO_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20, 
       height = 10)



# geothermal plasticity ---------------------------------------------------


# hyrbid plasticity -------------------------------------------------------


# overlap plasticity ------------------------------------------------------


