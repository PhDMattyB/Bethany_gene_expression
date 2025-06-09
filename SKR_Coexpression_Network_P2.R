##############################
## Gene coexpression analysis - part 2
##
## Matt Brachmann (PhDMattyB)
##
## 09.06.2025
##
##############################


setwd('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/')

library(tidyverse)
library(igraph)
library(ggraph)

library(readxl)
library(patchwork)
library(RColorBrewer)
library(viridis)

theme_set(theme_bw())

set.seed(1738)

brain_exp = read_tsv('Brain_Normalized_expression_counts.txt')
brain_genes = read_tsv('Brain_Normalized_expression_gene_list.txt')

brain_exp = bind_cols(brain_genes, 
                      brain_exp)


brain_limma = read_tsv('brain_limma_gene_list.txt')

brain_count_limma = inner_join(brain_exp, 
                               brain_limma, 
                               by = 'GeneID')



# liver_exp = read_tsv('Liver_Normalized_expression.txt')
# liver_exp_gene = read_tsv("Liver_Normalized_expression_gene_list.txt")
# 
# liver_exp = bind_cols(liver_exp_gene, 
#                       liver_exp)
# 
# 
# liver_limma = read_tsv("Liver_limma_gene_list.txt")
# 
# liver_count_limma = inner_join(liver_exp, 
#                                liver_limma, 
#                                by = 'GeneID')

metadata = names(brain_exp) %>% 
  as_tibble() %>% 
  slice(-1) %>% 
  separate(col = value, 
           into = c('ecotype', 
                    'temp', 
                    'family', 
                    'sample', 
                    'tissue'), 
           sep = '_', 
           remove = F) %>% 
  separate(col = ecotype, 
           into = c('sample_num', 
                    'ecotype'), 
           sep = '-') %>% 
  unite(col = ecotemp, 
        c('ecotype',
          'temp'),
        sep = '_',
        remove = F)


# brain_long = brain_count_limma %>% 
#   # rename(gene_ID = `...1`) %>% 
#   pivot_longer(cols = !GeneID, 
#                names_to = "library", 
#                values_to = "tpm") %>% 
#   mutate(logTPM = log10(tpm + 1)) 


# brain_long_wide = brain_long %>% 
#   select(GeneID, 
#          library, 
#          logTPM) %>% 
#   pivot_wider(names_from = library, 
#               values_from = logTPM)

brain_pca_data = brain_long_wide %>% 
  select(-GeneID)

brain_pca = prcomp(t(brain_count_limma[,-1]))
pc_importance = as.data.frame(t(summary(brain_pca)$importance))
head(pc_importance, 20)


brain_PCA_coord = brain_pca$x[, 1:10] %>% 
  as.data.frame() %>% 
  mutate(value = row.names(.)) %>% 
  full_join(metadata %>% 
              select(value,
                     sample_num, 
                     ecotemp, 
                     ecotype, 
                     temp, 
                     family), by = "value")

PCA_cols = c('#559fe4', 
             '#f94144',
             '#086375', 
             '#fee440',
             '#1565c0', 
             '#a0001c')


brain_PCA_coord %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = ecotemp), 
             color = "grey20", 
             shape = 21, 
             size = 3, 
             alpha = 0.8) +
  # scale_color_manual(values = PCA_cols)+
  scale_fill_manual(values = PCA_cols)+
  labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = NULL) +  
  theme_bw() +
  geom_hline(yintercept = 0, 
             col = 'black')+
  geom_vline(xintercept = 0, 
             col = 'black')+
  theme(panel.grid = element_blank(),
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )
