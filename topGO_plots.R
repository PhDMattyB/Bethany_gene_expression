### topGO SKR gene expression graphs

setwd('C:/Users/phdma/OneDrive/OneDrive - University of Glasgow/Documents/Parsons_Postdoc/SKR_Hybrid_Gene_expression')

library(tidyverse)

go_pal = c('#84a59d', 
           '#f28482', 
           '#f6bd60')

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

ggplot(data = div_12_topGO, 
       aes(y = GO_Term, 
           x = neg_log_pval))+
  geom_col(aes(fill = Type), 
           col = 'black')+
  labs(x = '-Log(p-value)')+
  scale_fill_manual(values = go_pal)+
  theme(axis.title.y = element_blank())




# amb geo div 18 ----------------------------------------------------------


