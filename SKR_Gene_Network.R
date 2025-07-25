##############################
## gene network analysis
##
## Matt Brachmann (PhDMattyB)
##
## 09.06.2025
##
##############################


setwd('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/')

library(tidyverse)
library(WGCNA)

brain_exp = read_tsv('Brain_Normalized_expression_counts.txt')
brain_genes = read_tsv('Brain_Normalized_expression_gene_list.txt')

brain_exp = bind_cols(brain_genes, 
                      brain_exp)


brain_limma = read_tsv('brain_limma_gene_list.txt')

brain_count_limma = inner_join(brain_exp, 
                               brain_limma, 
                               by = 'GeneID')



liver_exp = read_tsv('Liver_Normalized_expression.txt')
liver_exp_gene = read_tsv("Liver_Normalized_expression_gene_list.txt")

liver_exp = bind_cols(liver_exp_gene, 
                      liver_exp)


liver_limma = read_tsv("Liver_limma_gene_list.txt")

liver_count_limma = inner_join(liver_exp, 
                               liver_limma, 
                               by = 'GeneID')

# liver_all expression network ---------------------------------------------------------------

liver_Genes = liver_count_limma %>% select(GeneID)

liver_input_mat = liver_count_limma %>% 
  select(-GeneID) %>% 
  t()

# liver_input_mat = t(liver_count_limma) %>% 
#   as.data.frame()
# 
# liver_input_mat = liver_input_mat[-1,]

powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

liver_sft = pickSoftThreshold(
  liver_input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)


par(mfrow = c(1,2));
cex1 = 0.9


plot(liver_sft$fitIndices[, 1],
     -sign(liver_sft$fitIndices[, 3]) * liver_sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(liver_sft$fitIndices[, 1],
     -sign(liver_sft$fitIndices[, 3]) * liver_sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(liver_sft$fitIndices[, 1],
     liver_sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(liver_sft$fitIndices[, 1],
     liver_sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")


## pick 6 as it's where the plot asymptotes

picked_power = 6
temp_cor <- cor       
cor <- WGCNA::cor 

liver_netwk <- blockwiseModules(liver_input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)


mergedColors = labels2colors(liver_netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  liver_netwk$dendrograms[[1]],
  mergedColors[liver_netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )


liver_netwk$colors[liver_netwk$blockGenes[[1]]]
table(liver_netwk$colors)


liver_module_df <- data.frame(
  gene_id = names(liver_netwk$colors),
  colors = labels2colors(liver_netwk$colors))

liver_module_df = bind_cols(liver_Genes, 
                            liver_module_df) %>% 
  select(-gene_id) %>% 
  mutate(cluster = as.character(case_when(
    # colors == 'turquoise' ~ 'Cluster1', 
    colors == 'grey' ~ 'Cluster1'
  ))) %>% 
  select(-colors)



# write_delim(liver_module_df,
#             file = "liver_gene_modules.txt",
#             delim = "\t")

mod_eigen = moduleEigengenes(liver_input_mat, mergedColors)$eigengenes %>%
  rownames_to_column()
mod_eigen$ecotemp = row.names(mod_eigen)

mod_eigen = orderMEs(mod_eigen)

module_order = names(mod_eigen) %>% gsub("ME","", .)


liver_mME = mod_eigen %>%
  pivot_longer(-ecotemp) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )


liver_mME = liver_mME %>% 
  separate(col = ecotemp, 
                       into = c('ecotype', 
                                'temp', 
                                'family', 
                                'sample'), 
                       sep = '_', 
           remove = F) %>% 
  separate(col = ecotype, 
           into = c('sample_num', 
                    'ecotype'), 
           sep = '-', 
           remove = F) %>% 
  unite(col = ecotemp2, 
        c('ecotype',
          'temp'),
        sep = '_',
        remove = F) %>% 
  mutate(cluster = as.character(case_when(
          # name == 'turquoise' ~ 'Cluster1', 
          name == 'grey' ~ 'Cluster1'
        )))
  

liver_mME %>% ggplot(., aes(x=cluster, 
                            y=ecotype, 
                            fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  facet_grid(~temp)+
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", 
       y = "Modules", 
       fill="corr")


liver_details = liver_mME %>% 
  select(1:3) %>% 
  rename(name = ecotemp) %>% 
  distinct(name)


liver_exp_pattern = liver_count_limma %>%
  pivot_longer(-GeneID) %>% 
  inner_join(., 
             liver_module_df, 
             by = 'GeneID') %>% 
  separate(col = name, 
           into = c('ecotype', 
                    'temp', 
                    'family', 
                    'sample'), 
           sep = '_', 
           remove = F) %>% 
  separate(col = ecotype, 
           into = c('sample_num', 
                    'ecotype'), 
           sep = '-', 
           remove = F) %>% 
  unite(col = ecotemp, 
        c('ecotype',
          'temp'),
        sep = '_',
        remove = F)
  

liver_exp_pattern %>% 
  ggplot(., aes(x=ecotemp, 
                  y=value, 
                  group=GeneID)) +
  geom_line(aes(color = GeneID),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(cluster)) +
  labs(x = "treatment",
       y = "normalized expression")


liver_exp_pattern %>% 
  # filter(GeneID == 'ENSGACG00000010111') %>% 
  ggplot(., 
         aes(x = ecotype, 
             y = value))+
  geom_violin(aes(fill = temp), 
              color = 'black')+
  facet_wrap(~GeneID)




# brain all expression network --------------------------------------------

brain_Genes = brain_count_limma %>% select(GeneID)

brain_input_mat = brain_count_limma %>% 
  select(-GeneID) %>% 
  t()

# brain_input_mat = t(brain_count_limma) %>% 
#   as.data.frame()
# 
# brain_input_mat = brain_input_mat[-1,]

powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

brain_sft = pickSoftThreshold(
  brain_input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)


par(mfrow = c(1,2));
cex1 = 0.9


plot(brain_sft$fitIndices[, 1],
     -sign(brain_sft$fitIndices[, 3]) * brain_sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(brain_sft$fitIndices[, 1],
     -sign(brain_sft$fitIndices[, 3]) * brain_sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(brain_sft$fitIndices[, 1],
     brain_sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(brain_sft$fitIndices[, 1],
     brain_sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")


## pick 3 as it's where the plot asymptotes

picked_power = 3
temp_cor <- cor       
cor <- WGCNA::cor 

brain_netwk <- blockwiseModules(brain_input_mat,                # <= input here
                                
                                # == Adjacency Function ==
                                power = picked_power,                # <= power here
                                networkType = "signed",
                                
                                # == Tree and Block Options ==
                                deepSplit = 2,
                                pamRespectsDendro = F,
                                # detectCutHeight = 0.75,
                                minModuleSize = 30,
                                maxBlockSize = 4000,
                                
                                # == Module Adjustments ==
                                reassignThreshold = 0,
                                mergeCutHeight = 0.25,
                                
                                # == TOM == Archive the run results in TOM file (saves time)
                                saveTOMs = T,
                                saveTOMFileBase = "ER",
                                
                                # == Output Options
                                numericLabels = T,
                                verbose = 3)


mergedColors = labels2colors(brain_netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  brain_netwk$dendrograms[[1]],
  mergedColors[brain_netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )


brain_netwk$colors[brain_netwk$blockGenes[[1]]]
table(brain_netwk$colors)


brain_module_df <- data.frame(
  gene_id = names(brain_netwk$colors),
  colors = labels2colors(brain_netwk$colors))

brain_module_df = bind_cols(brain_Genes, 
                            brain_module_df) %>% 
  select(-gene_id) %>% 
  mutate(cluster = as.character(case_when(
    # colors == 'turquoise' ~ 'Cluster1', 
    colors == 'black' ~ 'Cluster 1', 
    colors == 'blue' ~ 'Cluster 2', 
    colors == 'brown' ~ 'Cluster 3', 
    colors == 'cyan' ~ 'Cluster 4', 
    colors == 'green' ~ 'Cluster 5', 
    colors == 'greenyellow' ~ 'Cluster 6', 
    colors == 'grey' ~ 'Cluster 7', 
    colors == 'magenta' ~ 'Cluster 8', 
    colors == 'midnightblue' ~ 'Cluster 9', 
    colors == 'pink' ~ 'Cluster 10', 
    colors == 'purple' ~ 'Cluster 11', 
    colors == 'red' ~ 'Cluster 12', 
    colors == 'salmon' ~ 'Cluster 13', 
    colors == 'tan' ~ 'Cluster 14', 
    colors == 'turquoise' ~ 'Cluster 15', 
    colors == 'yellow' ~ 'Cluster 16'
  ))) %>% 
  select(-colors)



# write_delim(brain_module_df,
#             file = "brain_gene_modules.txt",
#             delim = "\t")

mod_eigen = moduleEigengenes(brain_input_mat, mergedColors)$eigengenes 
# %>%
#   rownames_to_column()
mod_eigen$ecotemp = row.names(mod_eigen)

mod_eigen = orderMEs(mod_eigen)

module_order = names(mod_eigen) %>% gsub("ME","", .)


brain_mME = mod_eigen %>%
  pivot_longer(-ecotemp) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )


brain_mME = brain_mME %>% 
  separate(col = ecotemp, 
           into = c('ecotype', 
                    'temp', 
                    'family', 
                    'sample'), 
           sep = '_', 
           remove = F) %>% 
  separate(col = ecotype, 
           into = c('sample_num', 
                    'ecotype'), 
           sep = '-', 
           remove = F) %>% 
  unite(col = ecotemp2, 
        c('ecotype',
          'temp'),
        sep = '_',
        remove = F) %>% 
  mutate(cluster = as.character(case_when(
    # colors == 'turquoise' ~ 'Cluster1', 
    name == 'black' ~ 'Cluster 1', 
    name == 'blue' ~ 'Cluster 2', 
    name == 'brown' ~ 'Cluster 3', 
    name == 'cyan' ~ 'Cluster 4', 
    name == 'green' ~ 'Cluster 5', 
    name == 'greenyellow' ~ 'Cluster 6', 
    name == 'grey' ~ 'Cluster 7', 
    name == 'magenta' ~ 'Cluster 8', 
    name == 'midnightblue' ~ 'Cluster 9', 
    name == 'pink' ~ 'Cluster 10', 
    name == 'purple' ~ 'Cluster 11', 
    name == 'red' ~ 'Cluster 12', 
    name == 'salmon' ~ 'Cluster 13', 
    name == 'tan' ~ 'Cluster 14', 
    name == 'turquoise' ~ 'Cluster 15', 
    name == 'yellow' ~ 'Cluster 16'
  )))


brain_mME %>% ggplot(., aes(x=cluster, 
                            y=ecotemp2, 
                            fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  # facet_grid(~temp)+
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", 
       y = "Modules", 
       fill="corr")

brain_details = brain_mME %>% 
  select(ecotemp, 
         ecotemp2, 
         ecotype, 
         temp, 
         cluster) %>% 
  rename(name = ecotemp) %>% 
  distinct(name)


brain_exp_pattern = brain_count_limma %>%
  pivot_longer(-GeneID) %>% 
  inner_join(., 
             brain_module_df, 
             by = 'GeneID') %>% 
  separate(col = name, 
           into = c('ecotype', 
                    'temp', 
                    'family', 
                    'sample'), 
           sep = '_', 
           remove = F) %>% 
  separate(col = ecotype, 
           into = c('sample_num', 
                    'ecotype'), 
           sep = '-', 
           remove = F) %>% 
  unite(col = ecotemp, 
        c('ecotype',
          'temp'),
        sep = '_',
        remove = F)


brain_exp_pattern %>% 
  filter(cluster == 'Cluster 16') %>% 
  ggplot(., aes(x=ecotemp, 
                y=value, 
                group=GeneID)) +
  geom_line(aes(color = cluster),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(cluster)) +
  labs(x = "treatment",
       y = "normalized expression")






