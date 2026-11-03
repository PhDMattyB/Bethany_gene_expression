setwd('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/topGO_data/')

library(tidyverse)
library(topGO)

topDiffGenes <- function(allScore) {
  return(allScore == '1') # adjust based on your p-value cutoff
}


stick_GO_df = read_tsv('GCF_964276395.1-RS_2025_12_gene_ontology.gaf')


stick_GO_df %>%
  group_by(DB) %>% 
  dplyr::summarize(n = n())


View(stick_GO_df)



amb_plast_uniq_genes = read_tsv("Ambient_plastic_unique_genes.txt", 
                                col_names = c('Gene', 
                                              'Expression')) 

amb_plast_uniq_genelist = amb_plast_uniq_genes$Expression
names(amb_plast_uniq_genelist) <- amb_plast_uniq_genes$Gene


Trans_amb_hyb_12_genes = read_tsv("Trans_amb_hyb_12_gene_names_fixed.txt", 
                                  col_names = 'Gene')
Trans_amb_hyb_18_genes = read_tsv("Trans_amb_hyb_18_gene_names_fixed.txt", 
                                  col_names = 'Gene')

Trans_geo_hyb_12_genes = read_tsv("Trans_geo_hyb_12_gene_names_fixed.txt", 
                                  col_names = 'Gene')
Trans_geo_hyb_18_genes = read_tsv("Trans_geo_hyb_18_gene_names_fixed.txt", 
                                  col_names = 'Gene')

hyb_plast = read_tsv('Hybrid_plastic_genes.txt', 
                     col_names = 'Gene')

hyb_plast_genelist = hyb_plast$Gene

is.vector(hyb_plast_genelist)



# TopGO - AMB GEO DIVERGENCE ----------------------------------------------


# Divergence 12 degrees ---------------------------------------------------


div_12_genes = read_tsv("Brain_eco_div_12_genes.txt", 
                                col_names = c('Gene', 
                                              'Expression')) 

div_12_genelist = div_12_genes$Gene

geneID2GO = readMappings(file = "Stickleback_Gene_GO_Map.txt") 
all_genes_universe = names(geneID2GO)

all_genelist = rep(0, length(all_genes_universe))
names(all_genelist) = all_genes_universe
all_genelist[names(all_genelist) %in% div_12_genelist] <- 1
div_12_genelist <- as.factor(all_genelist)

div_12_GOdata_BP <- new("topGOdata",
                           ontology = 'BP',              # Choose "BP", "MF", or "CC"
                           allGenes = div_12_genelist,          # Named vector of genes
                           geneSelectionFun = topDiffGenes,
                           nodeSize = 5,                # Minimum number of annotated genes per GO term (e.g., >= 10)
                           annot = annFUN.gene2GO,       # Use custom gene-to-GO mapping
                           gene2GO = geneID2GO)

div_12_GOdata_MF <- new("topGOdata",
                           ontology = 'MF',              # Choose "BP", "MF", or "CC"
                           allGenes = div_12_genelist,          # Named vector of genes
                           geneSelectionFun = topDiffGenes,
                           nodeSize = 5,                # Minimum number of annotated genes per GO term (e.g., >= 10)
                           annot = annFUN.gene2GO,       # Use custom gene-to-GO mapping
                           gene2GO = geneID2GO)

div_12_GOdata_CC <- new("topGOdata",
                           ontology = 'CC',              # Choose "BP", "MF", or "CC"
                           allGenes = div_12_genelist,          # Named vector of genes
                           geneSelectionFun = topDiffGenes,
                           nodeSize = 5,                # Minimum number of annotated genes per GO term (e.g., >= 10)
                           annot = annFUN.gene2GO,       # Use custom gene-to-GO mapping
                           gene2GO = geneID2GO)

div_12_resultFisher_BP <- runTest(
  div_12_GOdata_BP,
  algorithm = "weight01",
  statistic = "fisher")

div_12_resultFisher_MF <- runTest(
  div_12_GOdata_MF,
  algorithm = "weight01",
  statistic = "fisher")

div_12_resultFisher_CC <- runTest(
  div_12_GOdata_CC,
  algorithm = "weight01",
  statistic = "fisher")

div_12_GO_BP <- GenTable(
  div_12_GOdata_BP,
  weight01Fisher = div_12_resultFisher_BP,
  orderBy = "weight01Fisher",
  ranksOf = "weight01Fisher",
  topNodes = length(score(div_12_resultFisher_BP)))

# Molecular Function
div_12_GO_MF <- GenTable(
  div_12_GOdata_MF,
  weight01Fisher = div_12_resultFisher_MF,
  orderBy = "weight01Fisher",
  ranksOf = "weight01Fisher",
  topNodes = length(score(div_12_resultFisher_MF)))

# Cellular Component
div_12_GO_CC <- GenTable(
  div_12_GOdata_CC,
  weight01Fisher = div_12_resultFisher_CC,
  orderBy = "weight01Fisher",
  ranksOf = "weight01Fisher",
  topNodes = length(score(div_12_resultFisher_CC)))


div_12_GO_BP = div_12_GO_BP %>% 
  as_tibble() %>% 
  # mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
  #         FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(weight01Fisher <= 0.05)

div_12_GO_MF = div_12_GO_MF %>% 
  as_tibble() %>% 
  # mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
  #         FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(weight01Fisher <= 0.05) 

div_12_GO_CC = div_12_GO_CC %>% 
  as_tibble() %>% 
  # mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
  #         FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(weight01Fisher <= 0.05)

div_12_BP_GO_genes = genesInTerm(
  div_12_GOdata_BP,
  div_12_GO_BP$GO.ID
)

div_12_BP_GO_genes_sig <- lapply(
  div_12_BP_GO_genes,
  function(genes) dplyr::intersect(genes, div_12_genes$Gene))


div_12_MF_GO_genes = genesInTerm(
  div_12_GOdata_MF,
  div_12_GO_MF$GO.ID
)

div_12_MF_GO_genes_sig <- lapply(
  div_12_MF_GO_genes,
  function(genes) dplyr::intersect(genes, div_12_genes$Gene))


div_12_CC_GO_genes = genesInTerm(
  div_12_GOdata_CC,
  div_12_GO_CC$GO.ID
)

div_12_CC_GO_genes_sig <- lapply(
  div_12_CC_GO_genes,
  function(genes) dplyr::intersect(genes, div_12_genes$Gene))



# Divergence 18 degrees ---------------------------------------------------

div_18_genes = read_tsv("Brain_eco_div_18_genes.txt", 
                        col_names = c('Gene', 
                                      'Expression')) 

div_18_genelist = div_18_genes$Gene

geneID2GO = readMappings(file = "Stickleback_Gene_GO_Map.txt") 
all_genes_universe = names(geneID2GO)

all_genelist = rep(0, length(all_genes_universe))
names(all_genelist) = all_genes_universe
all_genelist[names(all_genelist) %in% div_18_genelist] <- 1
div_18_genelist <- as.factor(all_genelist)

div_18_GOdata_BP <- new("topGOdata",
                        ontology = 'BP',              # Choose "BP", "MF", or "CC"
                        allGenes = div_18_genelist,          # Named vector of genes
                        geneSelectionFun = topDiffGenes,
                        nodeSize = 5,                # Minimum number of annotated genes per GO term (e.g., >= 10)
                        annot = annFUN.gene2GO,       # Use custom gene-to-GO mapping
                        gene2GO = geneID2GO)

div_18_GOdata_MF <- new("topGOdata",
                        ontology = 'MF',              # Choose "BP", "MF", or "CC"
                        allGenes = div_18_genelist,          # Named vector of genes
                        geneSelectionFun = topDiffGenes,
                        nodeSize = 5,                # Minimum number of annotated genes per GO term (e.g., >= 10)
                        annot = annFUN.gene2GO,       # Use custom gene-to-GO mapping
                        gene2GO = geneID2GO)

div_18_GOdata_CC <- new("topGOdata",
                        ontology = 'CC',              # Choose "BP", "MF", or "CC"
                        allGenes = div_18_genelist,          # Named vector of genes
                        geneSelectionFun = topDiffGenes,
                        nodeSize = 5,                # Minimum number of annotated genes per GO term (e.g., >= 10)
                        annot = annFUN.gene2GO,       # Use custom gene-to-GO mapping
                        gene2GO = geneID2GO)

div_18_resultFisher_BP <- runTest(
  div_18_GOdata_BP,
  algorithm = "weight01",
  statistic = "fisher")

div_18_resultFisher_MF <- runTest(
  div_18_GOdata_MF,
  algorithm = "weight01",
  statistic = "fisher")

div_18_resultFisher_CC <- runTest(
  div_18_GOdata_CC,
  algorithm = "weight01",
  statistic = "fisher")

div_18_GO_BP <- GenTable(
  div_18_GOdata_BP,
  weight01Fisher = div_18_resultFisher_BP,
  orderBy = "weight01Fisher",
  ranksOf = "weight01Fisher",
  topNodes = length(score(div_18_resultFisher_BP)))

# Molecular Function
div_18_GO_MF <- GenTable(
  div_18_GOdata_MF,
  weight01Fisher = div_18_resultFisher_MF,
  orderBy = "weight01Fisher",
  ranksOf = "weight01Fisher",
  topNodes = length(score(div_18_resultFisher_MF)))

# Cellular Component
div_18_GO_CC <- GenTable(
  div_18_GOdata_CC,
  weight01Fisher = div_18_resultFisher_CC,
  orderBy = "weight01Fisher",
  ranksOf = "weight01Fisher",
  topNodes = length(score(div_18_resultFisher_CC)))


div_18_GO_BP = div_18_GO_BP %>% 
  as_tibble() %>% 
  # mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
  #         FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(weight01Fisher <= 0.05)

div_18_GO_MF = div_18_GO_MF %>% 
  as_tibble() %>% 
  # mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
  #         FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(weight01Fisher <= 0.05) 

div_18_GO_CC = div_18_GO_CC %>% 
  as_tibble() %>% 
  # mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
  #         FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(weight01Fisher <= 0.05)

div_18_BP_GO_genes = genesInTerm(
  div_18_GOdata_BP,
  div_18_GO_BP$GO.ID
)

div_18_BP_GO_genes_sig <- lapply(
  div_18_BP_GO_genes,
  function(genes) dplyr::intersect(genes, div_18_genes$Gene))


div_18_MF_GO_genes = genesInTerm(
  div_18_GOdata_MF,
  div_18_GO_MF$GO.ID
)

div_18_MF_GO_genes_sig <- lapply(
  div_18_MF_GO_genes,
  function(genes) dplyr::intersect(genes, div_18_genes$Gene))


div_18_CC_GO_genes = genesInTerm(
  div_18_GOdata_CC,
  div_18_GO_CC$GO.ID
)

div_18_CC_GO_genes_sig <- lapply(
  div_18_CC_GO_genes,
  function(genes) dplyr::intersect(genes, div_18_genes$Gene))

# topGO plastic genes -----------------------------------------------------

# Ambient plasticity ------------------------------------------------------
amb_plast_uniq_genes = read_tsv("Ambient_plastic_unique_genes.txt", 
                                col_names = c('Gene', 
                                              'Expression')) 


amb_plast_genelist = amb_plast_uniq_genes$Gene

geneID2GO = readMappings(file = "Stickleback_Gene_GO_Map.txt") 
all_genes_universe = names(geneID2GO)

all_genelist = rep(0, length(all_genes_universe))
names(all_genelist) = all_genes_universe
all_genelist[names(all_genelist) %in% amb_plast_genelist] <- 1
amb_plast_genelist <- as.factor(all_genelist)

amb_plast_GOdata_BP <- new("topGOdata",
                           ontology = 'BP',              # Choose "BP", "MF", or "CC"
                           allGenes = amb_plast_genelist,          # Named vector of genes
                           geneSelectionFun = topDiffGenes,
                           nodeSize = 5,                # Minimum number of annotated genes per GO term (e.g., >= 10)
                           annot = annFUN.gene2GO,       # Use custom gene-to-GO mapping
                           gene2GO = geneID2GO)

amb_plast_GOdata_MF <- new("topGOdata",
                           ontology = 'MF',              # Choose "BP", "MF", or "CC"
                           allGenes = amb_plast_genelist,          # Named vector of genes
                           geneSelectionFun = topDiffGenes,
                           nodeSize = 5,                # Minimum number of annotated genes per GO term (e.g., >= 10)
                           annot = annFUN.gene2GO,       # Use custom gene-to-GO mapping
                           gene2GO = geneID2GO)

amb_plast_GOdata_CC <- new("topGOdata",
                           ontology = 'CC',              # Choose "BP", "MF", or "CC"
                           allGenes = amb_plast_genelist,          # Named vector of genes
                           geneSelectionFun = topDiffGenes,
                           nodeSize = 5,                # Minimum number of annotated genes per GO term (e.g., >= 10)
                           annot = annFUN.gene2GO,       # Use custom gene-to-GO mapping
                           gene2GO = geneID2GO)

amb_plast_resultFisher_BP <- runTest(
  amb_plast_GOdata_BP,
  algorithm = "weight01",
  statistic = "fisher")

amb_plast_resultFisher_MF <- runTest(
  amb_plast_GOdata_MF,
  algorithm = "weight01",
  statistic = "fisher")

amb_plast_resultFisher_CC <- runTest(
  amb_plast_GOdata_CC,
  algorithm = "weight01",
  statistic = "fisher")

amb_plast_GO_BP <- GenTable(
  amb_plast_GOdata_BP,
  weight01Fisher = amb_plast_resultFisher_BP,
  orderBy = "weight01Fisher",
  ranksOf = "weight01Fisher",
  topNodes = length(score(amb_plast_resultFisher_BP)))

# Molecular Function
amb_plast_GO_MF <- GenTable(
  amb_plast_GOdata_MF,
  weight01Fisher = amb_plast_resultFisher_MF,
  orderBy = "weight01Fisher",
  ranksOf = "weight01Fisher",
  topNodes = length(score(amb_plast_resultFisher_MF)))

# Cellular Component
amb_plast_GO_CC <- GenTable(
  amb_plast_GOdata_CC,
  weight01Fisher = amb_plast_resultFisher_CC,
  orderBy = "weight01Fisher",
  ranksOf = "weight01Fisher",
  topNodes = length(score(amb_plast_resultFisher_CC)))


amb_plast_GO_BP = amb_plast_GO_BP %>% 
  as_tibble() %>% 
  # mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
  #         FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(weight01Fisher <= 0.05)

amb_plast_GO_MF = amb_plast_GO_MF %>% 
  as_tibble() %>% 
  # mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
  #         FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(weight01Fisher <= 0.05) 

amb_plast_GO_CC = amb_plast_GO_CC %>% 
  as_tibble() %>% 
  # mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
  #         FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(weight01Fisher <= 0.05)

amb_plast_BP_GO_genes = genesInTerm(
  amb_plast_GOdata_BP,
  amb_plast_GO_BP$GO.ID
)

amb_plast_BP_GO_genes_sig <- lapply(
  amb_plast_BP_GO_genes,
  function(genes) intersect(genes, amb_plast_uniq_genes$Gene))


amb_plast_MF_GO_genes = genesInTerm(
  amb_plast_GOdata_MF,
  amb_plast_GO_MF$GO.ID
)

amb_plast_MF_GO_genes_sig <- lapply(
  amb_plast_MF_GO_genes,
  function(genes) intersect(genes, amb_plast_uniq_genes$Gene))


amb_plast_CC_GO_genes = genesInTerm(
  amb_plast_GOdata_CC,
  amb_plast_GO_CC$GO.ID
)

amb_plast_CC_GO_genes_sig <- lapply(
  amb_plast_CC_GO_genes,
  function(genes) intersect(genes, amb_plast_uniq_genes$Gene))




# Geothermal plasticity ---------------------------------------------------
geo_plast_uniq_genes = read_tsv("Geothermal_plastic_genes_unique.txt", 
                                col_names = c('Gene', 
                                              'Expression')) 

geo_plast_genelist = geo_plast_uniq_genes$Gene

geneID2GO = readMappings(file = "Stickleback_Gene_GO_Map.txt") 
all_genes_universe = names(geneID2GO)

all_genelist = rep(0, length(all_genes_universe))
names(all_genelist) = all_genes_universe
all_genelist[names(all_genelist) %in% geo_plast_genelist] <- 1
geo_plast_genelist <- as.factor(all_genelist)

geo_plast_GOdata_BP <- new("topGOdata",
                        ontology = 'BP',              # Choose "BP", "MF", or "CC"
                        allGenes = geo_plast_genelist,          # Named vector of genes
                        geneSelectionFun = topDiffGenes,
                        nodeSize = 5,                # Minimum number of annotated genes per GO term (e.g., >= 10)
                        annot = annFUN.gene2GO,       # Use custom gene-to-GO mapping
                        gene2GO = geneID2GO)

geo_plast_GOdata_MF <- new("topGOdata",
                        ontology = 'MF',              # Choose "BP", "MF", or "CC"
                        allGenes = geo_plast_genelist,          # Named vector of genes
                        geneSelectionFun = topDiffGenes,
                        nodeSize = 5,                # Minimum number of annotated genes per GO term (e.g., >= 10)
                        annot = annFUN.gene2GO,       # Use custom gene-to-GO mapping
                        gene2GO = geneID2GO)

geo_plast_GOdata_CC <- new("topGOdata",
                           ontology = 'CC',              # Choose "BP", "MF", or "CC"
                           allGenes = geo_plast_genelist,          # Named vector of genes
                           geneSelectionFun = topDiffGenes,
                           nodeSize = 5,                # Minimum number of annotated genes per GO term (e.g., >= 10)
                           annot = annFUN.gene2GO,       # Use custom gene-to-GO mapping
                           gene2GO = geneID2GO)

geo_plast_resultFisher_BP <- runTest(
  geo_plast_GOdata_BP,
  algorithm = "weight01",
  statistic = "fisher")

geo_plast_resultFisher_MF <- runTest(
  geo_plast_GOdata_MF,
  algorithm = "weight01",
  statistic = "fisher")

geo_plast_resultFisher_CC <- runTest(
  geo_plast_GOdata_CC,
  algorithm = "weight01",
  statistic = "fisher")

geo_plast_GO_BP <- GenTable(
  geo_plast_GOdata_BP,
  weight01Fisher = geo_plast_resultFisher_BP,
  orderBy = "weight01Fisher",
  ranksOf = "weight01Fisher",
  topNodes = length(score(geo_plast_resultFisher_BP)))

# Molecular Function
geo_plast_GO_MF <- GenTable(
  geo_plast_GOdata_MF,
  weight01Fisher = geo_plast_resultFisher_MF,
  orderBy = "weight01Fisher",
  ranksOf = "weight01Fisher",
  topNodes = length(score(geo_plast_resultFisher_MF)))

# Cellular Component
geo_plast_GO_CC <- GenTable(
  geo_plast_GOdata_CC,
  weight01Fisher = geo_plast_resultFisher_CC,
  orderBy = "weight01Fisher",
  ranksOf = "weight01Fisher",
  topNodes = length(score(geo_plast_resultFisher_CC)))


geo_plast_GO_BP = geo_plast_GO_BP %>% 
  as_tibble() %>% 
  # mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
  #         FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(weight01Fisher <= 0.05)

geo_plast_GO_MF = geo_plast_GO_MF %>% 
  as_tibble() %>% 
  # mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
  #         FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(weight01Fisher <= 0.05) 

geo_plast_GO_CC = geo_plast_GO_CC %>% 
  as_tibble() %>% 
  # mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
  #         FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(weight01Fisher <= 0.05)

geo_plast_BP_GO_genes = genesInTerm(
  geo_plast_GOdata_BP,
  geo_plast_GO_BP$GO.ID
)

geo_plast_BP_GO_genes_sig <- lapply(
  geo_plast_BP_GO_genes,
  function(genes) intersect(genes, geo_plast_uniq_genes$Gene))


geo_plast_MF_GO_genes = genesInTerm(
  geo_plast_GOdata_MF,
  geo_plast_GO_MF$GO.ID
)

geo_plast_MF_GO_genes_sig <- lapply(
  geo_plast_MF_GO_genes,
  function(genes) intersect(genes, geo_plast_uniq_genes$Gene))


geo_plast_CC_GO_genes = genesInTerm(
  geo_plast_GOdata_CC,
  geo_plast_GO_CC$GO.ID
)

geo_plast_CC_GO_genes_sig <- lapply(
  geo_plast_CC_GO_genes,
  function(genes) intersect(genes, geo_plast_uniq_genes$Gene))


# hybrid plasticity -------------------------------------------------------

hyb_plast_uniq_genes = read_tsv("Hybrid_plastic_genes_unique.txt", 
                                col_names = c('Gene', 
                                              'Expression')) 

hyb_plast_genelist = hyb_plast_uniq_genes$Gene

geneID2GO = readMappings(file = "Stickleback_Gene_GO_Map.txt") 
all_genes_universe = names(geneID2GO)

all_genelist = rep(0, length(all_genes_universe))
names(all_genelist) = all_genes_universe
all_genelist[names(all_genelist) %in% hyb_plast_genelist] <- 1
hyb_plast_genelist <- as.factor(all_genelist)

hyb_plast_GOdata_BP <- new("topGOdata",
                           ontology = 'BP',              # Choose "BP", "MF", or "CC"
                           allGenes = hyb_plast_genelist,          # Named vector of genes
                           geneSelectionFun = topDiffGenes,
                           nodeSize = 5,                # Minimum number of annotated genes per GO term (e.g., >= 10)
                           annot = annFUN.gene2GO,       # Use custom gene-to-GO mapping
                           gene2GO = geneID2GO)

hyb_plast_GOdata_MF <- new("topGOdata",
                           ontology = 'MF',              # Choose "BP", "MF", or "CC"
                           allGenes = hyb_plast_genelist,          # Named vector of genes
                           geneSelectionFun = topDiffGenes,
                           nodeSize = 5,                # Minimum number of annotated genes per GO term (e.g., >= 10)
                           annot = annFUN.gene2GO,       # Use custom gene-to-GO mapping
                           gene2GO = geneID2GO)

hyb_plast_GOdata_CC <- new("topGOdata",
                           ontology = 'CC',              # Choose "BP", "MF", or "CC"
                           allGenes = hyb_plast_genelist,          # Named vector of genes
                           geneSelectionFun = topDiffGenes,
                           nodeSize = 5,                # Minimum number of annotated genes per GO term (e.g., >= 10)
                           annot = annFUN.gene2GO,       # Use custom gene-to-GO mapping
                           gene2GO = geneID2GO)

hyb_plast_resultFisher_BP <- runTest(
  hyb_plast_GOdata_BP,
  algorithm = "weight01",
  statistic = "fisher")

hyb_plast_resultFisher_MF <- runTest(
  hyb_plast_GOdata_MF,
  algorithm = "weight01",
  statistic = "fisher")

hyb_plast_resultFisher_CC <- runTest(
  hyb_plast_GOdata_CC,
  algorithm = "weight01",
  statistic = "fisher")

# Classical Fisher's exact test
# Biological Process
hyb_plast_GO_BP <- GenTable(
  hyb_plast_GOdata_BP,
  weight01Fisher = hyb_plast_resultFisher_BP,
  orderBy = "weight01Fisher",
  ranksOf = "weight01Fisher",
  topNodes = length(score(hyb_plast_resultFisher_BP)))

# Molecular Function
hyb_plast_GO_MF <- GenTable(
  hyb_plast_GOdata_MF,
  weight01Fisher = hyb_plast_resultFisher_MF,
  orderBy = "weight01Fisher",
  ranksOf = "weight01Fisher",
  topNodes = length(score(hyb_plast_resultFisher_MF)))

# Cellular Component
hyb_plast_GO_CC <- GenTable(
  hyb_plast_GOdata_CC,
  weight01Fisher = hyb_plast_resultFisher_CC,
  orderBy = "weight01Fisher",
  ranksOf = "weight01Fisher",
  topNodes = length(score(hyb_plast_resultFisher_CC)))


hyb_plast_GO_BP = hyb_plast_GO_BP %>% 
  as_tibble() %>% 
  # mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
  #         FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(weight01Fisher <= 0.05)

hyb_plast_GO_MF = hyb_plast_GO_MF %>% 
  as_tibble() %>% 
  # mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
  #         FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(weight01Fisher <= 0.05) 

hyb_plast_GO_CC = hyb_plast_GO_CC %>% 
  as_tibble() %>% 
  # mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
  #         FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(weight01Fisher <= 0.05)

hyb_plast_BP_GO_genes = genesInTerm(
  hyb_plast_GOdata_BP,
  hyb_plast_GO_BP$GO.ID
)

hyb_plast_BP_GO_genes_sig <- lapply(
  hyb_plast_BP_GO_genes,
  function(genes) intersect(genes, hyb_plast_uniq_genes$Gene))


hyb_plast_MF_GO_genes = genesInTerm(
  hyb_plast_GOdata_MF,
  hyb_plast_GO_MF$GO.ID
)

hyb_plast_MF_GO_genes_sig <- lapply(
  hyb_plast_MF_GO_genes,
  function(genes) intersect(genes, hyb_plast_uniq_genes$Gene))


hyb_plast_CC_GO_genes = genesInTerm(
  hyb_plast_GOdata_CC,
  hyb_plast_GO_CC$GO.ID
)

hyb_plast_CC_GO_genes_sig <- lapply(
  hyb_plast_CC_GO_genes,
  function(genes) intersect(genes, hyb_plast_uniq_genes$Gene))



# Ecotype overlap plasticity ----------------------------------------------

overlap_plast_uniq_genes = read_tsv("Overlapping_plasticity_genes_allecotypes.txt", 
                                col_names = c('Gene', 
                                              'Expression')) 

overlap_plast_genelist = overlap_plast_uniq_genes$Gene

geneID2GO = readMappings(file = "Stickleback_Gene_GO_Map.txt") 
all_genes_universe = names(geneID2GO)

all_genelist = rep(0, length(all_genes_universe))
names(all_genelist) = all_genes_universe
all_genelist[names(all_genelist) %in% overlap_plast_genelist] <- 1
overlap_plast_genelist <- as.factor(all_genelist)

overlap_plast_GOdata_BP <- new("topGOdata",
                           ontology = 'BP',              # Choose "BP", "MF", or "CC"
                           allGenes = overlap_plast_genelist,          # Named vector of genes
                           geneSelectionFun = topDiffGenes,
                           nodeSize = 5,                # Minimum number of annotated genes per GO term (e.g., >= 10)
                           annot = annFUN.gene2GO,       # Use custom gene-to-GO mapping
                           gene2GO = geneID2GO)

overlap_plast_GOdata_MF <- new("topGOdata",
                           ontology = 'MF',              # Choose "BP", "MF", or "CC"
                           allGenes = overlap_plast_genelist,          # Named vector of genes
                           geneSelectionFun = topDiffGenes,
                           nodeSize = 5,                # Minimum number of annotated genes per GO term (e.g., >= 10)
                           annot = annFUN.gene2GO,       # Use custom gene-to-GO mapping
                           gene2GO = geneID2GO)

overlap_plast_GOdata_CC <- new("topGOdata",
                           ontology = 'CC',              # Choose "BP", "MF", or "CC"
                           allGenes = overlap_plast_genelist,          # Named vector of genes
                           geneSelectionFun = topDiffGenes,
                           nodeSize = 5,                # Minimum number of annotated genes per GO term (e.g., >= 10)
                           annot = annFUN.gene2GO,       # Use custom gene-to-GO mapping
                           gene2GO = geneID2GO)

overlap_plast_resultFisher_BP <- runTest(
  overlap_plast_GOdata_BP,
  algorithm = "weight01",
  statistic = "fisher")

overlap_plast_resultFisher_MF <- runTest(
  overlap_plast_GOdata_MF,
  algorithm = "weight01",
  statistic = "fisher")

overlap_plast_resultFisher_CC <- runTest(
  overlap_plast_GOdata_CC,
  algorithm = "weight01",
  statistic = "fisher")


# Biological Process
overlap_plast_GO_BP <- GenTable(
  overlap_plast_GOdata_BP,
  weight01Fisher = overlap_plast_resultFisher_BP,
  orderBy = "weight01Fisher",
  ranksOf = "weight01Fisher",
  topNodes = length(score(overlap_plast_resultFisher_BP)))

# Molecular Function
overlap_plast_GO_MF <- GenTable(
  overlap_plast_GOdata_MF,
  weight01Fisher = overlap_plast_resultFisher_MF,
  orderBy = "weight01Fisher",
  ranksOf = "weight01Fisher",
  topNodes = length(score(overlap_plast_resultFisher_MF)))

# Cellular Component
overlap_plast_GO_CC <- GenTable(
  overlap_plast_GOdata_CC,
  weight01Fisher = overlap_plast_resultFisher_CC,
  orderBy = "weight01Fisher",
  ranksOf = "weight01Fisher",
  topNodes = length(score(overlap_plast_resultFisher_CC)))


overlap_plast_GO_BP = overlap_plast_GO_BP %>% 
  as_tibble() %>% 
  # mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
  #         FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(weight01Fisher <= 0.05)

overlap_plast_GO_MF = overlap_plast_GO_MF %>% 
  as_tibble() %>% 
  # mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
  #         FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(weight01Fisher <= 0.05) 

overlap_plast_GO_CC = overlap_plast_GO_CC %>% 
  as_tibble() %>% 
  # mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
  #         FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(weight01Fisher <= 0.05)

overlap_plast_BP_GO_genes = genesInTerm(
  overlap_plast_GOdata_BP,
  overlap_plast_GO_BP$GO.ID
)

overlap_plast_BP_GO_genes_sig <- lapply(
  overlap_plast_BP_GO_genes,
  function(genes) intersect(genes, overlap_plast_uniq_genes$Gene))


overlap_plast_MF_GO_genes = genesInTerm(
  overlap_plast_GOdata_MF,
  overlap_plast_GO_MF$GO.ID
)

overlap_plast_MF_GO_genes_sig <- lapply(
  overlap_plast_MF_GO_genes,
  function(genes) intersect(genes, overlap_plast_uniq_genes$Gene))


overlap_plast_CC_GO_genes = genesInTerm(
  overlap_plast_GOdata_CC,
  overlap_plast_GO_CC$GO.ID
)

overlap_plast_CC_GO_genes_sig <- lapply(
  overlap_plast_CC_GO_genes,
  function(genes) intersect(genes, overlap_plast_uniq_genes$Gene))


