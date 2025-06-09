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
    axis.text = element_text(color = "black"), 
    legend.position = 'none')

ggsave('Brain_PCA_All_DEG_Axes1_2.tiff', 
       plot = last_plot(),
       dpi = 'retina', 
       width = 4, 
       height = 3)


brain_PCA_coord %>% 
  ggplot(aes(x = PC2, y = PC3)) +
  geom_point(aes(fill = ecotemp), 
             color = "grey20", 
             shape = 21, 
             size = 3, 
             alpha = 0.8) +
  # scale_color_manual(values = PCA_cols)+
  scale_fill_manual(values = PCA_cols)+
  labs(x = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC3 (", pc_importance[3, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = NULL) +  
  theme_bw() +
  geom_hline(yintercept = 0, 
             col = 'black')+
  geom_vline(xintercept = 0, 
             col = 'black')+
  theme(panel.grid = element_blank(),
        text = element_text(size= 14),
        axis.text = element_text(color = "black"))

ggsave('Brain_PCA_All_DEG_Axes2_3.tiff', 
       plot = last_plot(),
       dpi = 'retina', 
       width = 4, 
       height = 3)





## variance filtering

## grabbing the most variable genes
brain_var_genes = brain_count_limma %>% 
  pivot_longer(!GeneID) %>% 
  group_by(GeneID) %>% 
  summarize(var = var(value)) %>% 
  ungroup() %>% 
  filter(var > quantile(var, 0.5))

brain_var_genes_data = brain_count_limma %>% 
  filter(GeneID %in% brain_var_genes$GeneID) %>% 
  as.data.frame()

# brain_var_long = brain_var_genes_data %>% 
#   pivot_longer(!GeneID) %>% 
#   as.data.frame()

row.names(brain_var_genes_data) = brain_var_genes_data$GeneID


brain_cor_mat = cor(t(brain_var_genes_data[,-1]))

number_comparisons = ncol(brain_var_genes_data) - 1

brain_cor_mat_upper <- brain_cor_mat
brain_cor_mat_upper[lower.tri(brain_cor_mat_upper)] <- NA

brain_edge_table = brain_cor_mat_upper %>% 
  as.data.frame() %>% 
  mutate(from = row.names(brain_cor_mat)) %>% 
  pivot_longer(cols = !from, 
               names_to = 'to', 
               values_to = 'r') %>% 
  filter(is.na(r) == F) %>% 
  filter(from != to) %>% 
  mutate(t = r*sqrt((number_comparisons-2)/(1-r^2))) %>% 
  mutate(p.value = case_when(
    t > 0 ~ pt(t, df = number_comparisons-2, lower.tail = F),
    t <= 0 ~ pt(t, df = number_comparisons-2, lower.tail = F)
  )) %>% 
  mutate(FDR = p.adjust(p.value, method = 'fdr'))


brain_edge_table %>% 
  filter(r > 0) %>% 
  filter(FDR < 0.05) %>% 
  slice_min(order_by = abs(r), n = 10)


brain_edge_table %>% 
  # slice_sample(n = 20000) %>% 
  ggplot(aes(x = r)) +
  geom_histogram(color = "white", bins = 100) +
  geom_vline(xintercept = 0.7, color = "tomato1", size = 1.2) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )


brain_edge_table_select = brain_edge_table %>% 
  filter(r > 0.7)


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
  na.omit()


annotation_ensemble_genes = gene_annotation %>% 
  filter(feature == 'gene') %>% 
  pull(gene_id) %>% 
  as_tibble() %>% 
  separate(col = value, 
           into = c('ID', 
                    'Gene_name'), 
           sep = ";") %>% 
  # select(-parent) %>% 
  separate(col = ID, 
           into = c('trash', 
                    'ID'), 
           sep = '=') %>% 
  select(-trash) %>% 
  separate(col = Gene_name, 
           into = c('trash', 
                    'gene_name'), 
           sep = '=') %>% 
  select(-trash) %>% 
  rename(GeneID = ID)
# %>% 
#   select(ID) %>% 
#   rename(GeneID = ID)


brain_node_tab = data.frame(
  GeneID = c(brain_edge_table_select$from, 
             brain_edge_table_select$to) %>% 
    unique()
) %>% 
  left_join(annotation_ensemble_genes, 
            by = c('GeneID')) %>% 
  rename(functional_annotation = gene_name)


brain_network = graph_from_data_frame(
  brain_edge_table_select,
  vertices = brain_node_tab,
  directed = F
)

brain_modules <- cluster_leiden(brain_network, 
                                resolution = 2, 
                          objective_function = "modularity")



optimize_resolution <- function(network, resolution){
  modules = network %>% 
    cluster_leiden(resolution_parameter = resolution,
                   objective_function = "modularity")
  
  parsed_modules = data.frame(
    gene_ID = names(membership(modules)),
    module = as.vector(membership(modules)) 
  )
  
  num_module_5 = parsed_modules %>% 
    group_by(module) %>% 
    count() %>% 
    arrange(-n) %>% 
    filter(n >= 5) %>% 
    nrow() %>% 
    as.numeric()
  
  num_genes_contained = parsed_modules %>% 
    group_by(module) %>% 
    count() %>% 
    arrange(-n) %>% 
    filter(n >= 5) %>% 
    ungroup() %>% 
    summarise(sum = sum(n)) %>% 
    as.numeric()
  
  c(num_module_5, num_genes_contained)
  
}

optimization_results <- purrr::map_dfc(
  .x = seq(from = 0.25, to = 5, by = 0.25),
  .f = optimize_resolution, 
  network = brain_network
) %>% 
  t() %>% 
  cbind(
    resolution = seq(from = 0.25, to = 5, by = 0.25)
  ) %>% 
  as.data.frame() %>% 
  rename(num_module = V1,
         num_contained_gene = V2)


Optimize_num_module <- optimization_results %>% 
  ggplot(aes(x = resolution, y = num_module)) +
  geom_line(size = 1.1, alpha = 0.8, color = "dodgerblue2") +
  geom_point(size = 3, alpha = 0.7) +
  geom_vline(xintercept = 2, size = 1, linetype = 4) +
  labs(x = "resolution parameter",
       y = "num. modules\nw/ >=5 genes") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

Optimize_num_gene <- optimization_results %>% 
  ggplot(aes(x = resolution, y = num_contained_gene)) +
  geom_line(size = 1.1, alpha = 0.8, color = "violetred2") +
  geom_point(size = 3, alpha = 0.7) +
  geom_vline(xintercept = 2, size = 1, linetype = 4) +
  labs(x = "resolution parameter",
       y = "num. genes in\nmodules w/ >=5 genes") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

wrap_plots(Optimize_num_module, Optimize_num_gene, nrow = 2)


brain_network_modules <- data.frame(
  GeneID = names(membership(brain_modules)),
  module = as.vector(membership(brain_modules)) 
) %>% 
  inner_join(brain_node_tab, by = "GeneID")

brain_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5)

brain_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5) %>% 
  ungroup() %>% 
  summarise(sum = sum(n))


brain_modules_greater_5 <- brain_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5)

brain_network_modules <- brain_network_modules %>% 
  filter(module %in% brain_modules_greater_5$module)

brain_var_genes_data_long = brain_var_genes_data %>% 
  pivot_longer(!GeneID) %>% 
  as_tibble() %>% 
  # slice(-1) %>% 
  separate(col = name, 
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


brain_high_var_modules = brain_var_genes_data_long %>% 
  inner_join(brain_network_modules,
             by = 'GeneID')


brain_modules_mean_exp = brain_high_var_modules %>% 
  group_by(module, ecotype, temp) %>% 
  summarise(mean_exp = mean(value)) %>% 
  ungroup()


brain_module_peak_exp = brain_modules_mean_exp %>% 
  group_by(module) %>%
  slice_max(order_by = mean_exp, n = 1)

