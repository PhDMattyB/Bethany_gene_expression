setwd('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/topGO_data/')


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


topDiffGenes <- function(allScore) {
  return(allScore < 0.05) # adjust based on your p-value cutoff
}

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
                   topNodes = 50)



head(allRes)
write.table(allRes, file = "stickleback_topGO_results.txt", sep = "\t", row.names = FALSE)