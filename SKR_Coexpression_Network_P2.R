##############################
## Gene coexpression analysis - part 2
##
## Matt Brachmann (PhDMattyB)
##
## 09.06.2025
##
##############################

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


setwd('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/')

library(tidyverse)
library(igraph)
library(ggraph)

library(readxl)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(RCy3)

theme_set(theme_bw())

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



# PCA ---------------------------------------------------------------------


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




# variance filtering ------------------------------------------------------



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


# module optimization -----------------------------------------------------

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



# module expression -------------------------------------------------------



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
  group_by(module, ecotemp) %>% 
  summarise(mean_exp = mean(value)) %>% 
  ungroup()


brain_module_peak_exp = brain_modules_mean_exp %>% 
  group_by(module) %>%
  slice_max(order_by = mean_exp, n = 1) 

brain_high_var_modules %>% 
  filter(module == 5 | module == 6) %>%
ggplot(aes(x = ecotemp, y = value)) +
  geom_line(aes(group = GeneID), alpha = 0.3, color = "grey70") +
  geom_line(data = brain_modules_mean_exp %>%  
            filter(module == 5 | module == 6), 
            aes(x = ecotemp, 
                y = mean_exp, 
                group = module), 
            size = 2)+
  facet_grid(~module) 



# Heatmap -----------------------------------------------------------------


brain_modules_mean_exp$mean_exp %>% summary()

quantile(brain_modules_mean_exp$mean_exp, 0.95)

brain_modules_mean_exp = brain_modules_mean_exp %>% 
  mutate(mean_exp_clipped = case_when(
    mean_exp > 5.15 ~ 5.15, 
    mean_exp < -5.15 ~ -5.15, 
    T ~ mean_exp
  ))

brain_modules_mean_exp_reordered = brain_modules_mean_exp %>% 
  full_join(brain_module_peak_exp %>% 
              select(module,
                     mean_exp), 
            by = 'module') 

brain_modules_mean_exp_reordered %>% 
  group_by(module) %>% 
  summarize(n = n()) %>% View()

brain_heatmap = brain_modules_mean_exp_reordered %>% 
  separate(col = ecotemp, 
           into = c('ecotype', 
                    'temp'), 
           sep = '_', 
           remove = F) %>% 
  ggplot(aes(x = temp, 
             y = as.factor(module)))+
  facet_grid(.~ ecotype, 
             scales = "free", 
             space = "free") +
  geom_tile(aes(fill = mean_exp_clipped), 
            color = 'grey80')+
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")),
                       limits = c(-5.15, 5.15),
                       breaks = c(-5.15, 0, 5.15),
                       labels = c("< -5.15", "0", "> 5.15"))+
  labs(x = NULL,
       y = "Module",
       fill = "Normalized expression")+
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    # axis.text.x = element_blank(),
    axis.text.x = element_text(size = 14),
    strip.text = element_blank(),
    legend.position = "bottom",
    panel.spacing = unit(0.5, "lines") 
  )

heat_strip1 = expand.grid(
  ecotype = unique(metadata$ecotype),
  temp = unique(metadata$temp),
  stringsAsFactors = F
) %>% 
  mutate(ecotype = factor(ecotype, levels = c(
    "SKRC",
    "SKRHYB",
    "SKRW"))) %>% 
  mutate(temp = factor(temp, levels = c(
    "12",
    "18"))) %>% 
  ggplot(aes(x = ecotype, 
             y = 1)) +
  facet_grid(.~ ecotype, 
             scales = "free", 
             space = "free") +
  geom_tile(aes(fill = ecotype)) +
  scale_fill_manual(values = brewer.pal(8, "Set2")) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    panel.spacing = unit(0.5, "lines"),
    legend.key.height = unit(0.75, "lines")
  )

wrap_plots(brain_heatmap, 
           heat_strip1, 
           nrow = 2, 
           heights = c(1, 0.08, 0.08), 
           guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )



# network -----------------------------------------------------------------

subnetwork_edges = brain_edge_table_select %>% 
  # filter(from %in% names(neighbors_of_bait) &
  #          to %in% names(neighbors_of_bait)) %>% 
  group_by(from) %>% 
  slice_max(order_by = r, n = 5) %>% 
  ungroup() %>% 
  group_by(to) %>% 
  slice_max(order_by = r, n = 5) %>% 
  ungroup()

subnetwork_genes = c(subnetwork_edges$from, 
                     subnetwork_edges$to) %>% 
  unique()

length(subnetwork_genes)
dim(subnetwork_edges)


subnetwork_nodes <- brain_node_tab %>% 
  filter(GeneID %in% subnetwork_genes) %>% 
  left_join(brain_network_modules, by = "GeneID") %>% 
  left_join(brain_module_peak_exp, by = "module") %>% 
  mutate(module_annotation = case_when(
    str_detect(module, '1|2|3|4|5|6') ~ 'Mod 1-6', 
    str_detect(module, '7|8|10|13|14|15') ~ 'Mod 7-12', 
    str_detect(module, '16|18|19|28|29|30') ~ 'Mod 13-18', 
    str_detect(module, '33|36|43|45|47|52') ~ 'Mod 19-24', 
    str_detect(module, '53|60|64|68|72|73') ~ 'Mod 25-30', 
    str_detect(module, '76|79|82|89|96|97') ~ 'Mod 31-36', 
    str_detect(module, '101|102|126|130|134') ~ 'Mod 37-41', 
    str_detect(module, '189') ~ 'Mod 42'
  ))

dim(subnetwork_nodes)


brain_subnetwork <- graph_from_data_frame(subnetwork_edges,
                                       vertices = subnetwork_nodes,
                                       directed = F)



brain_subnetwork %>% 
  ggraph(layout = "kk",
         circular = F) +
  geom_edge_diagonal(color = "grey70", 
                     width = 0.5, 
                     alpha = 0.5) +
  geom_node_point(alpha = 0.8, 
                  color = "white", 
                  shape = 21, 
                  size = 2,
                  aes(fill = module_annotation)) + 
  scale_fill_manual(values = c(brewer.pal(8, "Accent")[c(1,3,6)], 
                               "grey30")) +
  labs(fill = "Modules") +
  guides(size = "none",
         fill = guide_legend(override.aes = list(size = 4), 
                             title.position = "top", nrow = 2)) +
  theme_void()+
  theme(
    text = element_text(size = 14), 
    legend.position = "bottom",
    legend.justification = 1,
    title = element_text(size = 12)
  )



# plasticity in expression networks ---------------------------------------

# Ambient plasticity - gene network ---------------------------------------
brain_plast_amb = read_csv('Brain_ambient_plastic_significant.csv')%>% 
  mutate(status = 'Outlier') %>% 
  left_join(., 
            brain_count_limma, 
            by = 'GeneID') %>% 
  select(GeneID, 
         9:56) %>% 
  as.data.frame()



# Ambient plasticity - Gene correlations ----------------------------------


row.names(brain_plast_amb) = brain_plast_amb$GeneID

amb_plast_cor_mat = cor(t(brain_plast_amb[,-1]))

number_comparisons = ncol(brain_plast_amb) - 1

plast_amb_cor_mat_upper = amb_plast_cor_mat
plast_amb_cor_mat_upper[lower.tri(plast_amb_cor_mat_upper)] <- NA

plast_amb_edge_table = plast_amb_cor_mat_upper %>% 
  as.data.frame() %>% 
  mutate(from = row.names(amb_plast_cor_mat)) %>% 
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


plast_amb_edge_table %>% 
  filter(r > 0) %>% 
  filter(FDR < 0.05) %>% 
  slice_min(order_by = abs(r), n = 10)


plast_amb_edge_table %>% 
  # slice_sample(n = 20000) %>% 
  ggplot(aes(x = r)) +
  geom_histogram(color = "white", bins = 100) +
  geom_vline(xintercept = 0.7, color = "tomato1", size = 1.2) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )


plast_amb_edge_table_select = plast_amb_edge_table %>% 
  filter(r > 0.7 | r < -0.7)







# Ambient plasticity - module optimization --------------------------------

plast_amb_node_tab = data.frame(
  GeneID = c(plast_amb_edge_table_select$from, 
             plast_amb_edge_table_select$to) %>% 
    unique()
) %>% 
  left_join(annotation_ensemble_genes, 
            by = c('GeneID')) %>% 
  rename(functional_annotation = gene_name)


plast_amb_network = graph_from_data_frame(
  plast_amb_edge_table_select,
  vertices = plast_amb_node_tab,
  directed = F
)

plast_amb_modules = cluster_leiden(plast_amb_network, 
                                resolution = 1.5, 
                                objective_function = "modularity")




plast_amb_optimization = purrr::map_dfc(
  .x = seq(from = 0.25, to = 5, by = 0.25),
  .f = optimize_resolution, 
  network = plast_amb_network) %>% 
  t() %>% 
  cbind(
    resolution = seq(from = 0.25, to = 5, by = 0.25)
  ) %>% 
  as.data.frame() %>% 
  rename(num_module = V1,
         num_contained_gene = V2)


plast_amb_optimize_num_module <- plast_amb_optimization %>% 
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

plast_amb_optimize_num_gene = plast_amb_optimization %>% 
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

wrap_plots(plast_amb_optimize_num_module, 
           plast_amb_optimize_num_gene, nrow = 2)





# Ambient plasticity - expression network ---------------------------------

plast_amb_network_modules <- data.frame(
  GeneID = names(membership(plast_amb_modules)),
  module = as.vector(membership(plast_amb_modules)) 
) %>% 
  inner_join(plast_amb_node_tab, by = "GeneID")

plast_amb_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5)

plast_amb_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5) %>% 
  ungroup() %>% 
  summarise(sum = sum(n))

# 
plast_amb_modules_greater_3 <- plast_amb_network_modules %>%
  group_by(module) %>%
  count() %>%
  arrange(-n) %>%
  filter(n >= 3)

plast_amb_network_modules <- plast_amb_network_modules %>%
  filter(module %in% plast_amb_modules_greater_3$module)

plast_amb_long = brain_plast_amb %>% 
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


plast_amb_high_var_modules = plast_amb_long %>% 
  inner_join(plast_amb_network_modules,
             by = 'GeneID')


plast_amb_modules_mean_exp = plast_amb_high_var_modules %>% 
  group_by(module, ecotemp) %>% 
  summarise(mean_exp = mean(value)) %>% 
  ungroup()


plast_amb_module_peak_exp = plast_amb_modules_mean_exp %>% 
  group_by(module) %>%
  slice_max(order_by = mean_exp, n = 1) 

plast_amb_high_var_modules %>% 
  # filter(module == 5 | module == 6) %>%
  ggplot(aes(x = ecotemp, y = value)) +
  geom_line(aes(group = GeneID), alpha = 0.3, color = "grey70") +
  geom_line(data = plast_amb_modules_mean_exp,  
              # filter(module == 5 | module == 6), 
            aes(x = ecotemp, 
                y = mean_exp, 
                group = module), 
            size = 2)+
  facet_grid(~module) 



# Ambient plasticity - heatmap --------------------------------------------

plast_amb_modules_mean_exp$mean_exp %>% summary()

quantile(plast_amb_modules_mean_exp$mean_exp, 0.95)

plast_amb_modules_mean_exp = plast_amb_modules_mean_exp %>% 
  mutate(mean_exp_clipped = case_when(
    mean_exp > 6.94 ~ 6.94, 
    mean_exp < -6.94 ~ -6.94, 
    T ~ mean_exp
  ))

plast_amb_modules_mean_exp_reordered = plast_amb_modules_mean_exp %>% 
  full_join(plast_amb_module_peak_exp %>% 
              select(module,
                     mean_exp), 
            by = 'module') 
# %>% 
#   mutate(module_rename = case_when(
#     module == '1' ~ '1', 
#     module == '2' ~ '2', 
#     module == '4' ~ '3', 
#     module == '5' ~ '4', 
#     module == '6' ~ '5', 
#     module == '8' ~ '6'
#   ))

# plast_amb_modules_mean_exp_reordered %>% 
#   group_by(module) %>% 
#   summarize(n = n()) %>% View()

plast_amb_heatmap = plast_amb_modules_mean_exp_reordered %>% 
  separate(col = ecotemp, 
           into = c('ecotype', 
                    'temp'), 
           sep = '_', 
           remove = F) %>% 
  ggplot(aes(x = temp, 
             y = as.factor(module)))+
  facet_grid(.~ ecotype, 
             scales = "free", 
             space = "free") +
  geom_tile(aes(fill = mean_exp_clipped), 
            color = 'grey80')+
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")),
                       limits = c(-6.94, 6.94),
                       breaks = c(-6.94, 0, 6.94),
                       labels = c("< -6.94", "0", "> 6.94"))+
  labs(x = NULL,
       y = "Module",
       fill = "Normalized expression")+
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    # axis.text.x = element_blank(),
    axis.text.x = element_text(size = 14),
    strip.text = element_blank(),
    legend.position = "bottom",
    panel.spacing = unit(0.5, "lines") 
  )

heat_strip1 = expand.grid(
  ecotype = unique(metadata$ecotype),
  temp = unique(metadata$temp),
  stringsAsFactors = F
) %>% 
  mutate(ecotype = factor(ecotype, levels = c(
    "SKRC",
    "SKRHYB",
    "SKRW"))) %>% 
  mutate(temp = factor(temp, levels = c(
    "12",
    "18"))) %>% 
  ggplot(aes(x = ecotype, 
             y = 1)) +
  facet_grid(.~ ecotype, 
             scales = "free", 
             space = "free") +
  geom_tile(aes(fill = ecotype)) +
  scale_fill_manual(values = brewer.pal(8, "Set2")) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    panel.spacing = unit(0.5, "lines"),
    legend.key.height = unit(0.75, "lines")
  )

wrap_plots(plast_amb_heatmap, 
           heat_strip1, 
           nrow = 2, 
           heights = c(1, 0.08, 0.08), 
           guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )


# Ambient plasticity - network --------------------------------------------

plast_amb_subnetwork_edges = plast_amb_edge_table_select %>% 
  # filter(from %in% names(neighbors_of_bait) &
  #          to %in% names(neighbors_of_bait)) %>% 
  group_by(from) %>% 
  slice_max(order_by = r, n = 3) %>% 
  ungroup() %>% 
  group_by(to) %>% 
  slice_max(order_by = r, n = 3) %>% 
  ungroup()

plast_amb_subnetwork_genes = c(plast_amb_subnetwork_edges$from, 
                     plast_amb_subnetwork_edges$to) %>% 
  unique()

# length(subnetwork_genes)
# dim(subnetwork_edges)


plast_amb_subnetwork_nodes <- plast_amb_node_tab %>% 
  filter(GeneID %in% plast_amb_subnetwork_genes) %>% 
  left_join(plast_amb_network_modules, by = "GeneID") %>% 
  left_join(plast_amb_module_peak_exp, by = "module") 

plast_amb_subnetwork_nodes$module = as.character(plast_amb_subnetwork_nodes$module)

# %>% 
#   mutate(module_rename = case_when(
#     module == '1' ~ '1', 
#     module == '2' ~ '2', 
#     module == '4' ~ '3', 
#     module == '5' ~ '4', 
#     module == '6' ~ '5', 
#     module == '8' ~ '6'
#   ))

dim(plast_amb_subnetwork_nodes)


plast_amb_subnetwork <- graph_from_data_frame(plast_amb_subnetwork_edges,
                                          vertices = plast_amb_subnetwork_nodes,
                                          directed = T)



ambient_plast_gene_network = plast_amb_subnetwork %>% 
  # ggraph()+
  ggraph(layout = "linear",
         circular = T) +
  # geom_edge_link(aes(color = factor(module_rename))) + 
  geom_edge_diagonal(color = "grey70", 
                     width = 0.5, 
                     alpha = 0.5) +
  # geom_node_point(alpha = 0.8, 
  #                 color = "white", 
  #                 shape = 21, 
  #                 size = 4,
  #                 aes(fill = module)) + 
  geom_node_point(alpha = 0.8, 
                  color = "black", 
                  shape = 21, 
                  size = 4,
                  fill = '#023e8a') + 
  scale_fill_manual(values = c(brewer.pal(8, "Accent"), 
                               "grey10")) +
  labs(fill = "Modules") +
  guides(size = "none",
         fill = guide_legend(override.aes = list(size = 4), 
                             title.position = "top", nrow = 2)) +
  theme_void()+
  theme(
    text = element_text(size = 14), 
    legend.position = "bottom",
    legend.justification = 1,
    title = element_text(size = 12)
  )


# Geo plasticity - gene network analysis ----------------------------------


brain_plast_geo = read_csv('Brain_geothermal_plastic_significant.csv')%>% 
  mutate(status = 'Outlier')%>% 
  left_join(., 
            brain_count_limma, 
            by = 'GeneID') %>% 
  select(GeneID, 
         9:56) %>% 
  as.data.frame()

# Geo plasticity - Gene correlations ----------------------------------


row.names(brain_plast_geo) = brain_plast_geo$GeneID

geo_plast_cor_mat = cor(t(brain_plast_geo[,-1]))

number_comparisons = ncol(brain_plast_geo) - 1

plast_geo_cor_mat_upper = geo_plast_cor_mat
plast_geo_cor_mat_upper[lower.tri(plast_geo_cor_mat_upper)] <- NA

plast_geo_edge_table = plast_geo_cor_mat_upper %>% 
  as.data.frame() %>% 
  mutate(from = row.names(geo_plast_cor_mat)) %>% 
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


plast_geo_edge_table %>% 
  filter(r > 0) %>% 
  filter(FDR < 0.05) %>% 
  slice_min(order_by = abs(r), n = 10)


plast_geo_edge_table %>% 
  # slice_sample(n = 20000) %>% 
  ggplot(aes(x = r)) +
  geom_histogram(color = "white", bins = 100) +
  geom_vline(xintercept = 0.7, color = "tomato1", size = 1.2) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )


plast_geo_edge_table_select = plast_geo_edge_table %>% 
  filter(r > 0.7 | r < -0.7)


# geo plasticity - module optimization --------------------------------

plast_geo_node_tab = data.frame(
  GeneID = c(plast_geo_edge_table_select$from, 
             plast_geo_edge_table_select$to) %>% 
    unique()
) %>% 
  left_join(annotation_ensemble_genes, 
            by = c('GeneID')) %>% 
  rename(functional_annotation = gene_name)


plast_geo_network = graph_from_data_frame(
  plast_geo_edge_table_select,
  vertices = plast_geo_node_tab,
  directed = F
)

plast_geo_modules = cluster_leiden(plast_geo_network, 
                                   resolution = 1.5, 
                                   objective_function = "modularity")




plast_geo_optimization = purrr::map_dfc(
  .x = seq(from = 0.25, to = 5, by = 0.25),
  .f = optimize_resolution, 
  network = plast_geo_network) %>% 
  t() %>% 
  cbind(
    resolution = seq(from = 0.25, to = 5, by = 0.25)
  ) %>% 
  as.data.frame() %>% 
  rename(num_module = V1,
         num_contained_gene = V2)


plast_geo_optimize_num_module <- plast_geo_optimization %>% 
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

plast_geo_optimize_num_gene = plast_geo_optimization %>% 
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

wrap_plots(plast_geo_optimize_num_module, 
           plast_geo_optimize_num_gene, nrow = 2)





# geo plasticity - expression network ---------------------------------

plast_geo_network_modules <- data.frame(
  GeneID = names(membership(plast_geo_modules)),
  module = as.vector(membership(plast_geo_modules)) 
) %>% 
  inner_join(plast_geo_node_tab, by = "GeneID")

plast_geo_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5)

plast_geo_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5) %>% 
  ungroup() %>% 
  summarise(sum = sum(n))

# 
plast_geo_modules_greater_3 <- plast_geo_network_modules %>%
  group_by(module) %>%
  count() %>%
  arrange(-n) %>%
  filter(n >= 3)

plast_geo_network_modules <- plast_geo_network_modules %>%
  filter(module %in% plast_geo_modules_greater_3$module)

plast_geo_long = brain_plast_geo %>% 
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


plast_geo_high_var_modules = plast_geo_long %>% 
  inner_join(plast_geo_network_modules,
             by = 'GeneID')


plast_geo_modules_mean_exp = plast_geo_high_var_modules %>% 
  group_by(module, ecotemp) %>% 
  summarise(mean_exp = mean(value)) %>% 
  ungroup()


plast_geo_module_peak_exp = plast_geo_modules_mean_exp %>% 
  group_by(module) %>%
  slice_max(order_by = mean_exp, n = 1) 

plast_geo_high_var_modules %>% 
  # filter(module == 5 | module == 6) %>%
  ggplot(aes(x = ecotemp, y = value)) +
  geom_line(aes(group = GeneID), alpha = 0.3, color = "grey70") +
  geom_line(data = plast_geo_modules_mean_exp,  
            # filter(module == 5 | module == 6), 
            aes(x = ecotemp, 
                y = mean_exp, 
                group = module), 
            size = 2)+
  facet_grid(~module) 



# geo plasticity - heatmap --------------------------------------------

plast_geo_modules_mean_exp$mean_exp %>% summary()

quantile(plast_geo_modules_mean_exp$mean_exp, 0.95)

plast_geo_modules_mean_exp = plast_geo_modules_mean_exp %>% 
  mutate(mean_exp_clipped = case_when(
    mean_exp > 6.94 ~ 6.94, 
    mean_exp < -6.94 ~ -6.94, 
    T ~ mean_exp
  ))

plast_geo_modules_mean_exp_reordered = plast_geo_modules_mean_exp %>% 
  full_join(plast_geo_module_peak_exp %>% 
              select(module,
                     mean_exp), 
            by = 'module') 
# %>% 
#   mutate(module_rename = case_when(
#     module == '1' ~ '1', 
#     module == '2' ~ '2', 
#     module == '4' ~ '3', 
#     module == '5' ~ '4', 
#     module == '6' ~ '5', 
#     module == '8' ~ '6'
#   ))

# plast_geo_modules_mean_exp_reordered %>% 
#   group_by(module) %>% 
#   summarize(n = n()) %>% View()

plast_geo_heatmap = plast_geo_modules_mean_exp_reordered %>% 
  separate(col = ecotemp, 
           into = c('ecotype', 
                    'temp'), 
           sep = '_', 
           remove = F) %>% 
  ggplot(aes(x = temp, 
             y = as.factor(module)))+
  facet_grid(.~ ecotype, 
             scales = "free", 
             space = "free") +
  geom_tile(aes(fill = mean_exp_clipped), 
            color = 'grey80')+
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")),
                       limits = c(-6.94, 6.94),
                       breaks = c(-6.94, 0, 6.94),
                       labels = c("< -6.94", "0", "> 6.94"))+
  labs(x = NULL,
       y = "Module",
       fill = "Normalized expression")+
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    # axis.text.x = element_blank(),
    axis.text.x = element_text(size = 14),
    strip.text = element_blank(),
    legend.position = "bottom",
    panel.spacing = unit(0.5, "lines") 
  )

heat_strip1 = expand.grid(
  ecotype = unique(metadata$ecotype),
  temp = unique(metadata$temp),
  stringsAsFactors = F
) %>% 
  mutate(ecotype = factor(ecotype, levels = c(
    "SKRC",
    "SKRHYB",
    "SKRW"))) %>% 
  mutate(temp = factor(temp, levels = c(
    "12",
    "18"))) %>% 
  ggplot(aes(x = ecotype, 
             y = 1)) +
  facet_grid(.~ ecotype, 
             scales = "free", 
             space = "free") +
  geom_tile(aes(fill = ecotype)) +
  scale_fill_manual(values = brewer.pal(8, "Set2")) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    panel.spacing = unit(0.5, "lines"),
    legend.key.height = unit(0.75, "lines")
  )

wrap_plots(plast_geo_heatmap, 
           heat_strip1, 
           nrow = 2, 
           heights = c(1, 0.08, 0.08), 
           guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )


# geo plasticity - network --------------------------------------------

plast_geo_subnetwork_edges = plast_geo_edge_table_select %>% 
  # filter(from %in% names(neighbors_of_bait) &
  #          to %in% names(neighbors_of_bait)) %>% 
  group_by(from) %>% 
  slice_max(order_by = r, n = 3) %>% 
  ungroup() %>% 
  group_by(to) %>% 
  slice_max(order_by = r, n = 3) %>% 
  ungroup()

plast_geo_subnetwork_genes = c(plast_geo_subnetwork_edges$from, 
                               plast_geo_subnetwork_edges$to) %>% 
  unique()

# length(subnetwork_genes)
# dim(subnetwork_edges)


plast_geo_subnetwork_nodes <- plast_geo_node_tab %>% 
  filter(GeneID %in% plast_geo_subnetwork_genes) %>% 
  left_join(plast_geo_network_modules, by = "GeneID") %>% 
  left_join(plast_geo_module_peak_exp, by = "module") 

plast_geo_subnetwork_nodes$module = as.character(plast_geo_subnetwork_nodes$module)

# %>% 
#   mutate(module_rename = case_when(
#     module == '1' ~ '1', 
#     module == '2' ~ '2', 
#     module == '4' ~ '3', 
#     module == '5' ~ '4', 
#     module == '6' ~ '5', 
#     module == '8' ~ '6'
#   ))

dim(plast_geo_subnetwork_nodes)


plast_geo_subnetwork <- graph_from_data_frame(plast_geo_subnetwork_edges,
                                              vertices = plast_geo_subnetwork_nodes,
                                              directed = T)



plast_geo_gene_network = plast_geo_subnetwork %>% 
  # ggraph()+
  ggraph(layout = "linear",
         circular = T) +
  # geom_edge_link(aes(color = factor(module_rename))) + 
  geom_edge_diagonal(color = "grey70", 
                     width = 0.5, 
                     alpha = 0.5) +
  # geom_node_point(alpha = 0.8, 
  #                 color = "white", 
  #                 shape = 21, 
  #                 size = 4,
  #                 aes(fill = module)) + 
  geom_node_point(alpha = 0.8, 
                  color = "black", 
                  shape = 21, 
                  size = 4,
                  fill = '#e63946') + 
  scale_fill_manual(values = c(brewer.pal(8, "Accent"), 
                               "grey10")) +
  labs(fill = "Modules") +
  guides(size = "none",
         fill = guide_legend(override.aes = list(size = 4), 
                             title.position = "top", nrow = 2)) +
  theme_void()+
  theme(
    text = element_text(size = 14), 
    legend.position = "bottom",
    legend.justification = 1,
    title = element_text(size = 12)
  )


# hybrid plasticity - gene network analysis -------------------------------


brain_plast_hyb = read_csv('Brain_hybrid_plastic_significant.csv') %>% 
  mutate(status = 'Outlier')%>% 
  left_join(., 
            brain_count_limma, 
            by = 'GeneID') %>% 
  select(GeneID, 
         9:56) %>% 
  as.data.frame()

# hyb plasticity - Gene correlations ----------------------------------


row.names(brain_plast_hyb) = brain_plast_hyb$GeneID

hyb_plast_cor_mat = cor(t(brain_plast_hyb[,-1]))

number_comparisons = ncol(brain_plast_hyb) - 1

plast_hyb_cor_mat_upper = hyb_plast_cor_mat
plast_hyb_cor_mat_upper[lower.tri(plast_hyb_cor_mat_upper)] <- NA

plast_hyb_edge_table = plast_hyb_cor_mat_upper %>% 
  as.data.frame() %>% 
  mutate(from = row.names(hyb_plast_cor_mat)) %>% 
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


plast_hyb_edge_table %>% 
  filter(r > 0) %>% 
  filter(FDR < 0.05) %>% 
  slice_min(order_by = abs(r), n = 10)


plast_hyb_edge_table %>% 
  # slice_sample(n = 20000) %>% 
  ggplot(aes(x = r)) +
  geom_histogram(color = "white", bins = 100) +
  geom_vline(xintercept = 0.7, color = "tomato1", size = 1.2) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )


plast_hyb_edge_table_select = plast_hyb_edge_table %>% 
  filter(r > 0.7 | r < -0.7)


# hyb plasticity - module optimization --------------------------------

plast_hyb_node_tab = data.frame(
  GeneID = c(plast_hyb_edge_table_select$from, 
             plast_hyb_edge_table_select$to) %>% 
    unique()
) %>% 
  left_join(annotation_ensemble_genes, 
            by = c('GeneID')) %>% 
  rename(functional_annotation = gene_name)


plast_hyb_network = graph_from_data_frame(
  plast_hyb_edge_table_select,
  vertices = plast_hyb_node_tab,
  directed = F
)

plast_hyb_modules = cluster_leiden(plast_hyb_network, 
                                   resolution = 1.5, 
                                   objective_function = "modularity")




plast_hyb_optimization = purrr::map_dfc(
  .x = seq(from = 0.25, to = 5, by = 0.25),
  .f = optimize_resolution, 
  network = plast_hyb_network) %>% 
  t() %>% 
  cbind(
    resolution = seq(from = 0.25, to = 5, by = 0.25)
  ) %>% 
  as.data.frame() %>% 
  rename(num_module = V1,
         num_contained_gene = V2)


plast_hyb_optimize_num_module <- plast_hyb_optimization %>% 
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

plast_hyb_optimize_num_gene = plast_hyb_optimization %>% 
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

wrap_plots(plast_hyb_optimize_num_module, 
           plast_hyb_optimize_num_gene, nrow = 2)





# hyb plasticity - expression network ---------------------------------

plast_hyb_network_modules <- data.frame(
  GeneID = names(membership(plast_hyb_modules)),
  module = as.vector(membership(plast_hyb_modules)) 
) %>% 
  inner_join(plast_hyb_node_tab, by = "GeneID")

plast_hyb_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5)

plast_hyb_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5) %>% 
  ungroup() %>% 
  summarise(sum = sum(n))

# 
plast_hyb_modules_greater_3 <- plast_hyb_network_modules %>%
  group_by(module) %>%
  count() %>%
  arrange(-n) %>%
  filter(n >= 3)

plast_hyb_network_modules <- plast_hyb_network_modules %>%
  filter(module %in% plast_hyb_modules_greater_3$module)

plast_hyb_long = brain_plast_hyb %>% 
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


plast_hyb_high_var_modules = plast_hyb_long %>% 
  inner_join(plast_hyb_network_modules,
             by = 'GeneID')


plast_hyb_modules_mean_exp = plast_hyb_high_var_modules %>% 
  group_by(module, ecotemp) %>% 
  summarise(mean_exp = mean(value)) %>% 
  ungroup()


plast_hyb_module_peak_exp = plast_hyb_modules_mean_exp %>% 
  group_by(module) %>%
  slice_max(order_by = mean_exp, n = 1) 

plast_hyb_high_var_modules %>% 
  # filter(module == 5 | module == 6) %>%
  ggplot(aes(x = ecotemp, y = value)) +
  geom_line(aes(group = GeneID), alpha = 0.3, color = "grey70") +
  geom_line(data = plast_hyb_modules_mean_exp,  
            # filter(module == 5 | module == 6), 
            aes(x = ecotemp, 
                y = mean_exp, 
                group = module), 
            size = 2)+
  facet_grid(~module) 



# hyb plasticity - heatmap --------------------------------------------

plast_hyb_modules_mean_exp$mean_exp %>% summary()

quantile(plast_hyb_modules_mean_exp$mean_exp, 0.95)

plast_hyb_modules_mean_exp = plast_hyb_modules_mean_exp %>% 
  mutate(mean_exp_clipped = case_when(
    mean_exp > 6.94 ~ 6.94, 
    mean_exp < -6.94 ~ -6.94, 
    T ~ mean_exp
  ))

plast_hyb_modules_mean_exp_reordered = plast_hyb_modules_mean_exp %>% 
  full_join(plast_hyb_module_peak_exp %>% 
              select(module,
                     mean_exp), 
            by = 'module') 
# %>% 
#   mutate(module_rename = case_when(
#     module == '1' ~ '1', 
#     module == '2' ~ '2', 
#     module == '4' ~ '3', 
#     module == '5' ~ '4', 
#     module == '6' ~ '5', 
#     module == '8' ~ '6'
#   ))

# plast_hyb_modules_mean_exp_reordered %>% 
#   group_by(module) %>% 
#   summarize(n = n()) %>% View()

plast_hyb_heatmap = plast_hyb_modules_mean_exp_reordered %>% 
  separate(col = ecotemp, 
           into = c('ecotype', 
                    'temp'), 
           sep = '_', 
           remove = F) %>% 
  ggplot(aes(x = temp, 
             y = as.factor(module)))+
  facet_grid(.~ ecotype, 
             scales = "free", 
             space = "free") +
  geom_tile(aes(fill = mean_exp_clipped), 
            color = 'grey80')+
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")),
                       limits = c(-6.94, 6.94),
                       breaks = c(-6.94, 0, 6.94),
                       labels = c("< -6.94", "0", "> 6.94"))+
  labs(x = NULL,
       y = "Module",
       fill = "Normalized expression")+
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    # axis.text.x = element_blank(),
    axis.text.x = element_text(size = 14),
    strip.text = element_blank(),
    legend.position = "bottom",
    panel.spacing = unit(0.5, "lines") 
  )

heat_strip1 = expand.grid(
  ecotype = unique(metadata$ecotype),
  temp = unique(metadata$temp),
  stringsAsFactors = F
) %>% 
  mutate(ecotype = factor(ecotype, levels = c(
    "SKRC",
    "SKRHYB",
    "SKRW"))) %>% 
  mutate(temp = factor(temp, levels = c(
    "12",
    "18"))) %>% 
  ggplot(aes(x = ecotype, 
             y = 1)) +
  facet_grid(.~ ecotype, 
             scales = "free", 
             space = "free") +
  geomm_tile(aes(fill = ecotype)) +
  scale_fill_manual(values = brewer.pal(8, "Set2")) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    panel.spacing = unit(0.5, "lines"),
    legend.key.height = unit(0.75, "lines")
  )

wrap_plots(plast_hyb_heatmap, 
           heat_strip1, 
           nrow = 2, 
           heights = c(1, 0.08, 0.08), 
           guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )


# hyb plasticity - network --------------------------------------------

plast_hyb_subnetwork_edges = plast_hyb_edge_table_select %>% 
  # filter(from %in% names(neighbors_of_bait) &
  #          to %in% names(neighbors_of_bait)) %>% 
  group_by(from) %>% 
  slice_max(order_by = r, n = 3) %>% 
  ungroup() %>% 
  group_by(to) %>% 
  slice_max(order_by = r, n = 3) %>% 
  ungroup()

plast_hyb_subnetwork_genes = c(plast_hyb_subnetwork_edges$from, 
                               plast_hyb_subnetwork_edges$to) %>% 
  unique()

# length(subnetwork_genes)
# dim(subnetwork_edges)


plast_hyb_subnetwork_nodes <- plast_hyb_node_tab %>% 
  filter(GeneID %in% plast_hyb_subnetwork_genes) %>% 
  left_join(plast_hyb_network_modules, by = "GeneID") %>% 
  left_join(plast_hyb_module_peak_exp, by = "module") 

plast_hyb_subnetwork_nodes$module = as.character(plast_hyb_subnetwork_nodes$module)

# %>% 
#   mutate(module_rename = case_when(
#     module == '1' ~ '1', 
#     module == '2' ~ '2', 
#     module == '4' ~ '3', 
#     module == '5' ~ '4', 
#     module == '6' ~ '5', 
#     module == '8' ~ '6'
#   ))

dim(plast_hyb_subnetwork_nodes)


plast_hyb_subnetwork <- graph_from_data_frame(plast_hyb_subnetwork_edges,
                                              vertices = plast_hyb_subnetwork_nodes,
                                              directed = T)



plast_hyb_gene_network = plast_hyb_subnetwork %>% 
  # ggraph()+
  ggraph(layout = "linear",
         circular = T) +
  # hybm_edge_link(aes(color = factor(module_rename))) + 
  geom_edge_diagonal(color = "grey70", 
                     width = 0.5, 
                     alpha = 0.5) +
  # hybm_node_point(alpha = 0.8, 
  #                 color = "white", 
  #                 shape = 21, 
  #                 size = 4,
  #                 aes(fill = module)) + 
  geom_node_point(alpha = 0.8, 
                  color = "black", 
                  shape = 21, 
                  size = 4,
                  fill = '#90be6d') + 
  scale_fill_manual(values = c(brewer.pal(8, "Accent"), 
                               "grey10")) +
  labs(fill = "Modules") +
  guides(size = "none",
         fill = guide_legend(override.aes = list(size = 4), 
                             title.position = "top", nrow = 2)) +
  theme_void()+
  theme(
    text = element_text(size = 14), 
    legend.position = "bottom",
    legend.justification = 1,
    title = element_text(size = 12)
  )


(ambient_plast_gene_network|plast_geo_gene_network)/(plast_hyb_gene_network)

wrap_plots(ambient_plast_gene_network, 
           plast_geo_gene_network, 
           plast_hyb_gene_network,
           nrow = 3, 
           heights = c(1, 0.08, 0.08), 
           guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )



# Combine the three graphs ------------------------------------------------


## Need to output the Igraph files as Networks. 
## These can then be put into Cytoscape to be merged
## and to better visualize the network
createNetworkFromIgraph(plast_amb_subnetwork)

createNetworkFromIgraph(plast_geo_subnetwork)

createNetworkFromIgraph(plast_hyb_subnetwork)





# ecological divergence networks ------------------------------------------

# amb vs geo Divergence at 12 degrees -------------------------------------------------

brain_eco12 = read_csv('Brain_eco_div_12.csv')%>% 
  filter(adj.P.Val <= 0.05) %>% 
  mutate(status = 'Outlier') %>% 
  left_join(., 
            brain_count_limma, 
            by = 'GeneID') %>% 
  select(GeneID, 
         9:56) %>% 
  as.data.frame()







row.names(brain_eco12) = brain_eco12$GeneID

amb_plast_cor_mat = cor(t(brain_eco12[,-1]))

number_comparisons = ncol(brain_eco12) - 1

eco12_cor_mat_upper = amb_plast_cor_mat
eco12_cor_mat_upper[lower.tri(eco12_cor_mat_upper)] <- NA

eco12_edge_table = eco12_cor_mat_upper %>% 
  as.data.frame() %>% 
  mutate(from = row.names(amb_plast_cor_mat)) %>% 
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


eco12_edge_table %>% 
  filter(r > 0) %>% 
  filter(FDR < 0.05) %>% 
  slice_min(order_by = abs(r), n = 10)


eco12_edge_table %>% 
  # slice_sample(n = 20000) %>% 
  ggplot(aes(x = r)) +
  geom_histogram(color = "white", bins = 100) +
  geom_vline(xintercept = 0.7, color = "tomato1", size = 1.2) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )


eco12_edge_table_select = eco12_edge_table %>% 
  filter(r > 0.7 | r < -0.7)










eco12_node_tab = data.frame(
  GeneID = c(eco12_edge_table_select$from, 
             eco12_edge_table_select$to) %>% 
    unique()
) %>% 
  left_join(annotation_ensemble_genes, 
            by = c('GeneID')) %>% 
  rename(functional_annotation = gene_name)


eco12_network = graph_from_data_frame(
  eco12_edge_table_select,
  vertices = eco12_node_tab,
  directed = F
)

eco12_modules = cluster_leiden(eco12_network, 
                                   resolution = 1, 
                                   objective_function = "modularity")

eco12_optimization = purrr::map_dfc(
  .x = seq(from = 0.25, to = 5, by = 0.25),
  .f = optimize_resolution, 
  network = eco12_network) %>% 
  t() %>% 
  cbind(
    resolution = seq(from = 0.25, to = 5, by = 0.25)
  ) %>% 
  as.data.frame() %>% 
  rename(num_module = V1,
         num_contained_gene = V2)


eco12_optimize_num_module <- eco12_optimization %>% 
  ggplot(aes(x = resolution, y = num_module)) +
  geom_line(size = 1.1, alpha = 0.8, color = "dodgerblue2") +
  geom_point(size = 3, alpha = 0.7) +
  geom_vline(xintercept = 1, size = 1, linetype = 4) +
  labs(x = "resolution parameter",
       y = "num. modules\nw/ >=5 genes") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

eco12_optimize_num_gene = eco12_optimization %>% 
  ggplot(aes(x = resolution, y = num_contained_gene)) +
  geom_line(size = 1.1, alpha = 0.8, color = "violetred2") +
  geom_point(size = 3, alpha = 0.7) +
  geom_vline(xintercept = 1, size = 1, linetype = 4) +
  labs(x = "resolution parameter",
       y = "num. genes in\nmodules w/ >=5 genes") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

wrap_plots(eco12_optimize_num_module, 
           eco12_optimize_num_gene, nrow = 2)








eco12_network_modules <- data.frame(
  GeneID = names(membership(eco12_modules)),
  module = as.vector(membership(eco12_modules)) 
) %>% 
  inner_join(eco12_node_tab, by = "GeneID")

eco12_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5)

eco12_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5) %>% 
  ungroup() %>% 
  summarise(sum = sum(n))
eco12_modules_greater_3 <- eco12_network_modules %>%
  group_by(module) %>%
  count() %>%
  arrange(-n) %>%
  filter(n >= 3)
eco12_network_modules <- eco12_network_modules %>%
  filter(module %in% eco12_modules_greater_3$module)
eco12_long = brain_eco12 %>% 
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
eco12_high_var_modules = eco12_long %>% 
  inner_join(eco12_network_modules,
             by = 'GeneID')
eco12_modules_mean_exp = eco12_high_var_modules %>% 
  group_by(module, ecotemp) %>% 
  summarise(mean_exp = mean(value)) %>% 
  ungroup()

eco12_module_peak_exp = eco12_modules_mean_exp %>% 
  group_by(module) %>%
  slice_max(order_by = mean_exp, n = 1) 
# 
# eco12_high_var_modules %>% 
#   # filter(module == 5 | module == 6) %>%
#   ggplot(aes(x = ecotemp, y = value)) +
#   geom_line(aes(group = GeneID), alpha = 0.3, color = "grey70") +
#   geom_line(data = eco12_modules_mean_exp,  
#             # filter(module == 5 | module == 6), 
#             aes(x = ecotemp, 
#                 y = mean_exp, 
#                 group = module), 
#             size = 2)+
#   facet_grid(~module) 






eco12_subnetwork_edges = eco12_edge_table_select %>% 
  # filter(from %in% names(neighbors_of_bait) &
  #          to %in% names(neighbors_of_bait)) %>% 
  group_by(from) %>% 
  slice_max(order_by = r, n = 3) %>% 
  ungroup() %>% 
  group_by(to) %>% 
  slice_max(order_by = r, n = 3) %>% 
  ungroup()

eco12_subnetwork_genes = c(eco12_subnetwork_edges$from, 
                               eco12_subnetwork_edges$to) %>% 
  unique()

# length(eco12_subnetwork_genes)
# dim(eco12_subnetwork_edges)


eco12_subnetwork_nodes <- eco12_node_tab %>% 
  filter(GeneID %in% eco12_subnetwork_genes) %>% 
  left_join(eco12_network_modules, by = "GeneID") %>% 
  left_join(eco12_module_peak_exp, by = "module") 

eco12_subnetwork_nodes$module = as.character(eco12_subnetwork_nodes$module)

# %>% 
#   mutate(module_rename = case_when(
#     module == '1' ~ '1', 
#     module == '2' ~ '2', 
#     module == '4' ~ '3', 
#     module == '5' ~ '4', 
#     module == '6' ~ '5', 
#     module == '8' ~ '6'
#   ))

dim(eco12_subnetwork_nodes)


eco12_subnetwork <- graph_from_data_frame(eco12_subnetwork_edges,
                                              vertices = eco12_subnetwork_nodes,
                                              directed = T)



eco12_divergence_gene_network = eco12_subnetwork %>% 
  # ggraph()+
  ggraph(layout = "linear",
         circular = T) +
  # geom_edge_link(aes(color = factor(module_rename))) + 
  geom_edge_diagonal(color = "grey70", 
                     width = 0.5, 
                     alpha = 0.5) +
  # geom_node_point(alpha = 0.8, 
  #                 color = "white", 
  #                 shape = 21, 
  #                 size = 4,
  #                 aes(fill = module)) + 
  geom_node_point(alpha = 0.8, 
                  color = "black", 
                  shape = 21, 
                  size = 4,
                  fill = '#023e8a') + 
  scale_fill_manual(values = c(brewer.pal(8, "Accent"), 
                               "grey10")) +
  labs(fill = "Modules") +
  guides(size = "none",
         fill = guide_legend(override.aes = list(size = 4), 
                             title.position = "top", nrow = 2)) +
  theme_void()+
  theme(
    text = element_text(size = 14), 
    legend.position = "bottom",
    legend.justification = 1,
    title = element_text(size = 12)
  )



# Divergence amb vs hyb @12 degrees  --------------------------------------


