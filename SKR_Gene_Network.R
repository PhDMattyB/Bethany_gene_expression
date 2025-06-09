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

brain_exp = read_tsv('brain_gene_read_counts_table_all_final.tsv')
brain_limma = read_tsv('brain_limma_gene_list.txt')

brain_count_limma = inner_join(brain_exp, 
           brain_limma, 
           by = 'GeneID')


liver_exp = read_tsv('liver_gene_read_counts_table_all_final.tsv')
liver_limma = read_tsv("Liver_limma_gene_list.txt")

liver_count_limma = inner_join(liver_exp, 
                               liver_limma, 
                               by = 'GeneID')


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




# liver_all ---------------------------------------------------------------

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
    colors == 'turquoise' ~ 'Cluster1', 
    colors == 'grey' ~ 'Cluster2'
  ))) %>% 
  select(-colors)



write_delim(liver_module_df,
            file = "liver_gene_modules.txt",
            delim = "\t")

mod_eigen = moduleEigengenes(liver_input_mat, mergedColors)$eigengenes

mod_eigen = orderMEs(mod_eigen)

module_order = names(mod_eigen) %>% gsub("ME","", .)

MEs0$treatment = row.names(MEs0)

mod_eigen$ecotemp = row.names(mod_eigen)

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
                       sep = '_') %>% 
  separate(col = ecotype, 
           into = c('sample_num', 
                    'ecotype'), 
           sep = '-') %>% 
  unite(col = ecotemp, 
        c('ecotype',
          'temp'),
        sep = '_',
        remove = F) %>% 
  mutate(cluster = as.character(case_when(
          name == 'turquoise' ~ 'Cluster1', 
          name == 'grey' ~ 'Cluster2'
        )))
  

liver_mME %>% ggplot(., aes(x=ecotemp, 
                            y=cluster, 
                            fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", 
       y = "Modules", 
       fill="corr")
