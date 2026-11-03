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

# Classical Fisher's exact test
amb_plast_resultFisher_BP = runTest(amb_plast_GOdata_BP, 
                                    algorithm = "classic", 
                                    statistic = "fisher")
amb_plast_resultFisher_MF = runTest(amb_plast_GOdata_MF, 
                                    algorithm = "classic", 
                                    statistic = "fisher")
amb_plast_resultFisher_CC = runTest(amb_plast_GOdata_CC, 
                                    algorithm = "classic", 
                                    statistic = "fisher")

amb_plast_GO_BP <- GenTable(amb_plast_GOdata_BP, 
                            classicFisher = amb_plast_resultFisher_BP, 
                            # weight01Fisher = resultWeight01,
                            orderBy = "classicFisher", 
                            ranksOf = "classicFisher", 
                            topNodes = length(score(amb_plast_resultFisher_BP)))

amb_plast_GO_MF <- GenTable(amb_plast_GOdata_MF, 
                            classicFisher = amb_plast_resultFisher_MF, 
                            # weight01Fisher = resultWeight01,
                            orderBy = "classicFisher", 
                            ranksOf = "classicFisher", 
                            topNodes = length(score(amb_plast_resultFisher_MF)))

amb_plast_GO_CC <- GenTable(amb_plast_GOdata_CC, 
                            classicFisher = amb_plast_resultFisher_CC, 
                            # weight01Fisher = resultWeight01,
                            orderBy = "classicFisher", 
                            ranksOf = "classicFisher", 
                            topNodes = length(score(amb_plast_resultFisher_CC)))

amb_plast_GO_BP = amb_plast_GO_BP %>% 
  as_tibble() %>% 
  mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
          FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(FDR <= 0.05)

amb_plast_GO_MF = amb_plast_GO_MF %>% 
  as_tibble() %>% 
  mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
          FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(FDR <= 0.05) 

amb_plast_GO_CC = amb_plast_GO_CC %>% 
  as_tibble() %>% 
  mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
          FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(FDR <= 0.05)




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

# Classical Fisher's exact test
geo_plast_resultFisher_BP = runTest(geo_plast_GOdata_BP, 
                                 algorithm = "classic", 
                                 statistic = "fisher")
geo_plast_resultFisher_MF = runTest(geo_plast_GOdata_MF, 
                                 algorithm = "classic", 
                                 statistic = "fisher")
geo_plast_resultFisher_CC = runTest(geo_plast_GOdata_CC, 
                                 algorithm = "classic", 
                                 statistic = "fisher")

geo_plast_GO_BP <- GenTable(geo_plast_GOdata_BP, 
                         classicFisher = geo_plast_resultFisher_BP, 
                         # weight01Fisher = resultWeight01,
                         orderBy = "classicFisher", 
                         ranksOf = "classicFisher", 
                         topNodes = length(score(geo_plast_resultFisher_BP)))

geo_plast_GO_MF <- GenTable(geo_plast_GOdata_MF, 
                         classicFisher = geo_plast_resultFisher_MF, 
                         # weight01Fisher = resultWeight01,
                         orderBy = "classicFisher", 
                         ranksOf = "classicFisher", 
                         topNodes = length(score(geo_plast_resultFisher_MF)))

geo_plast_GO_CC <- GenTable(geo_plast_GOdata_CC, 
                         classicFisher = geo_plast_resultFisher_CC, 
                         # weight01Fisher = resultWeight01,
                         orderBy = "classicFisher", 
                         ranksOf = "classicFisher", 
                         topNodes = length(score(geo_plast_resultFisher_CC)))

geo_plast_GO_BP = geo_plast_GO_BP %>% 
  as_tibble() %>% 
  mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
    FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(FDR <= 0.05)

geo_plast_GO_MF = geo_plast_GO_MF %>% 
  as_tibble() %>% 
  mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
          FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(FDR <= 0.05) 

geo_plast_GO_CC = geo_plast_GO_CC %>% 
  as_tibble() %>% 
  mutate( pvalue = as.numeric(sub("^<\\s*", "", classicFisher)),
          FDR = p.adjust(pvalue, method = "BH")) %>% 
  filter(FDR <= 0.05)


GO_genes = genesInTerm(
  geo_plast_GOdata_MF,
  geo_plast_GO_MF$GO.ID
)

GO_sig_genes <- lapply(
  GO_genes,
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


GO_genes = genesInTerm(
  hyb_plast_GOdata_MF,
  hyb_plast_GO_MF$GO.ID
)

GO_sig_genes <- lapply(
  GO_genes,
  function(genes) intersect(genes, hyb_plast_uniq_genes$Gene))



# Ecotype overlap plasticity ----------------------------------------------






# stick_GO_df %>% 
#   dplyr::select(Symbol, 
#                 GO_ID) %>% 
#   write_tsv('Stickleback_Gene_GO_Map.txt')


geneID2GO <- readMappings(file = "Stickleback_Gene_GO_Map.txt") 

all_genes_universe <- names(geneID2GO)
geneList <- rep(0, length(all_genes_universe))
names(geneList) <- all_genes_universe
geneList[names(geneList) %in% hyb_plast_genelist] <- 1
geneList <- as.factor(geneList)




GOdata <- new("topGOdata",
              ontology = "BP",              # Choose "BP", "MF", or "CC"
              allGenes = geneList,          # Named vector of genes
              geneSelectionFun = topDiffGenes,
              nodeSize = 5,                # Minimum number of annotated genes per GO term (e.g., >= 10)
              annot = annFUN.gene2GO,       # Use custom gene-to-GO mapping
              gene2GO = geneID2GO)



# Classical Fisher's exact test
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

# Or a more conservative tree-aware algorithm like weight01
resultWeight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")



allRes <- GenTable(GOdata, 
                   classicFisher = resultFisher, 
                   weight01Fisher = resultWeight01,
                   orderBy = "weight01Fisher", 
                   ranksOf = "classicFisher", 
                   topNodes = 100)

allRes %>% 
  as_tibble() %>% 
  filter(weight01Fisher <= 0.05) %>% View()


head(allRes)
write.table(allRes, file = "stickleback_topGO_results.txt", sep = "\t", row.names = FALSE)