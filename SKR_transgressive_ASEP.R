##############################
## Transgressive allele specific expression 
##
## Matt Brachmann (PhDMattyB)
##
## 15.07.2024
##
##############################


setwd('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/')

library(tidyverse)


# brain_12_snps = read.table('brain_12oC_sigASEPSNPs.csv', 
#            header = T) %>% 
#   as_tibble() %>% 
#   dplyr::select(SNP, 
#                 p.value) %>% 
#   rename(Chrom_pos_allele_REF_ALT = SNP) %>% 
#   inner_join(., 
#              ecotype_snps, 
#              by = 'Chrom_pos_allele_REF_ALT')

# gene annotations --------------------------------------------------------

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


gene_metadata = gene_annotation %>% 
  # filter(feature == 'gene') %>% 
  dplyr::select(position,
                chromosome,
                feature, 
                start, 
                end)

ensemlbe_annotation_data = gene_annotation %>% 
  # filter(feature == 'gene') %>% 
  pull(gene_id) %>% 
  as_tibble() %>% 
  separate(value, 
           into = c('ensemble_id', 
                    'gene_name',
                    'relationship'), 
           sep = ';')

ensemble_annotation_genes = ensemlbe_annotation_data %>% 
  separate(ensemble_id, 
           into = c('trash', 
                    'ensemble_id'), 
           sep = '=') %>%
  separate(gene_name, 
           into = c('trash', 
                    'gene_name'), 
           sep = '=') %>% 
  dplyr::select(-trash) %>% 
  separate(ensemble_id, 
           into = c('ensemble_name', 
                    'trash'), 
           sep = '.CDS') %>% 
  dplyr::select(-trash) 

annotation_data = bind_cols(gene_metadata, 
                            ensemble_annotation_genes) %>% 
  rename(CHR = chromosome, 
         BP = position) %>% 
  rename(gene_ensembl = ensemble_name)


# ecotype snps ------------------------------------------------------------

ecotype_snps = read_csv('ecotype_unique_snps_all.csv') %>% 
  separate(col = Chrom_pos_allele_REF_ALT, 
           into = c('CHR', 
                    'BP', 
                    'REF1', 
                    'ALT1'), 
           sep = '_') %>% 
  inner_join(., 
             annotation_data, 
             by = 'gene_ensembl', 
             relationship = 'many-to-many')

ecotype_snps_fixed = ecotype_snps %>% 
  dplyr::select(CHR.x, 
         CHR.y, 
         BP.x, 
         BP.y, 
         start, 
         end,
         REF, 
         ALT, 
         relationship, 
         feature, 
         other_affected_genes, 
         ecotype) %>% 
  na.omit() %>% 
  separate(relationship, 
           into = c('trash', 
                    'GeneID'), 
           sep = '=') %>% 
  dplyr::select(-trash)

# Brain ASEP data ---------------------------------------------------------


brain_12_snps = read.table('brain_12oC_sigASEPSNPs.csv', 
           header = T) %>% 
  as_tibble() %>% 
  separate(gene, 
           into = c('gene_ensembl', 
                    'the_rest'),
           sep = ',') %>% 
  dplyr::select(-the_rest) %>% 
  inner_join(.,
             ecotype_snps,
             by = 'gene_ensembl',
             relationship = 'many-to-many')%>%
  # inner_join(., 
  #            ecotype_snps, 
  #            by = 'gene_ensembl')%>%
  # distinct(gene_ensembl) %>% 
  separate(col = SNP, 
           into = c('CHR', 
                    'BP', 
                    'REF', 
                    'ALT'), 
           sep = '_')

  
brain_12_snps %>% 
  group_by(ecotype) %>% 
  summarize(n = n())


# brain_18_snps = read.table('brain_18oC_sigASEPSNPs.csv', 
#                            header = T) %>% 
#   as_tibble() %>% 
#   dplyr::select(SNP, 
#                 p.value) %>% 
#   rename(Chrom_pos_allele_REF_ALT = SNP) %>% 
#   inner_join(., 
#              ecotype_snps, 
#              by = 'Chrom_pos_allele_REF_ALT')

brain_18_snps = read.table('brain_18oC_sigASEPSNPs.csv', 
                           header = T) %>% 
  as_tibble() %>% 
  separate(gene, 
           into = c('gene_ensembl', 
                    'the_rest'),
           sep = ',') %>% 
  dplyr::select(-the_rest) %>% 
  inner_join(., 
             ecotype_snps, 
             by = 'gene_ensembl', 
             relationship = 'many-to-many')%>% 
  separate(col = SNP, 
           into = c('CHR', 
                    'BP', 
                    'REF', 
                    'ALT'), 
           sep = '_')


brain_18_snps %>% 
  group_by(ecotype) %>% 
  summarize(n = n())


# Transgressive expression in hybrids -------------------------------------

brain_exp = read_tsv('Brain_Normalized_expression_counts.txt')
brain_genes = read_tsv('Brain_Normalized_expression_gene_list.txt')

brain_exp = bind_cols(brain_genes, 
                      brain_exp)


brain_limma = read_tsv('brain_limma_gene_list.txt')

brain_count_limma = inner_join(brain_exp, 
                               brain_limma, 
                               by = 'GeneID')

Gene_ID_12 = read_csv("Brain_Transgressive_expression_12degrees.csv")

### amb vs hyb expression at 12 degrees
Trans_amb_hyb_12 = read_csv('Brain_amb_hyb_12_div.csv')%>% 
  filter(adj.P.Val <= 0.05) %>% 
  # mutate(status = 'Outlier') %>% 
  inner_join(., 
             Gene_ID_12) %>% 
  dplyr::select(GeneID, 
                logFC, 
                adj.P.Val) %>% 
  rename(gene_ensembl = GeneID) %>% 
  inner_join(., 
           annotation_data) %>% 
  rename(GeneID = gene_name)

Trans_amb_hyb_12$BP = as.character(Trans_amb_hyb_12$BP)

## combine with the 12 degree SNPS with allele specific expression



# ecotype_snps$BP = as.character(ecotype_snps$BP)
annotation_data$BP = as.character(annotation_data$BP)

# inner_join(ecotype_snps, 
#            annotation_data, 
#            by = c('CHR', 
#                   'BP'))

Trans_amb_hyb_12_snps = inner_join(Trans_amb_hyb_12, 
           ecotype_snps_fixed, 
           by = 'GeneID', 
           relationship = 'many-to-many') %>% 
  dplyr::select(gene_ensembl, 
         logFC, 
         adj.P.Val, 
         CHR, 
         GeneID,
         feature.x, 
         start.x, 
         end.x, 
         BP.x, 
         REF, 
         ALT, 
         other_affected_genes, 
         ecotype)

Trans_amb_hyb_12_snps %>% 
  # group_by(GeneID, 
  #          ecotype) %>% 
  group_by(GeneID) %>% 
  summarize(n = n()) 


Trans_geo_hyb_12 = read_csv('Brain_geo_hyb_12_div.csv')%>% 
  filter(adj.P.Val <= 0.05) %>% 
  inner_join(., 
             Gene_ID_12) %>% 
  dplyr::select(GeneID, 
                logFC, 
                adj.P.Val) %>% 
  rename(gene_ensembl = GeneID)%>% 
  inner_join(., 
             annotation_data) %>% 
  rename(GeneID = gene_name)

trans_geo_hyb_12_snps = inner_join(Trans_geo_hyb_12, 
           ecotype_snps_fixed, 
           by = 'GeneID', 
           relationship = 'many-to-many')%>% 
  dplyr::select(gene_ensembl, 
                logFC, 
                adj.P.Val, 
                CHR, 
                GeneID,
                feature.x, 
                start.x, 
                end.x, 
                BP.x, 
                REF, 
                ALT, 
                other_affected_genes, 
                ecotype)

trans_geo_hyb_12_snps %>% 
  # group_by(GeneID, 
  #          ecotype) %>% 
  group_by(GeneID) %>% 
  summarize(n = n()) 

## lets try 18 degrees and see what happens

Gene_ID_18 = read.csv('Brain_Transgressive_expression_18degrees.csv')

Trans_amb_hyb_18 = read_csv('Brain_amb_hyb_18_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  inner_join(., 
             Gene_ID_18) %>% 
  dplyr::select(GeneID, 
                logFC, 
                adj.P.Val) %>% 
  rename(gene_ensembl = GeneID)%>% 
  inner_join(., 
             annotation_data) %>% 
  rename(GeneID = gene_name)
  
Trans_amb_hyb_18_snps = inner_join(Trans_amb_hyb_18, 
           ecotype_snps_fixed, 
           by = 'GeneID', 
           relationship = 'many-to-many')%>% 
  dplyr::select(gene_ensembl, 
                logFC, 
                adj.P.Val, 
                CHR, 
                GeneID,
                feature.x, 
                start.x, 
                end.x, 
                BP.x, 
                REF, 
                ALT, 
                other_affected_genes, 
                ecotype)

Trans_amb_hyb_18_snps %>% 
  # group_by(GeneID, 
  #          ecotype) %>% 
  group_by(GeneID) %>% 
  summarize(n = n()) 


Trans_geo_hyb_18 = read_csv('Brain_geo_hyb_18_div.csv') %>% 
  filter(adj.P.Val <= 0.05) %>% 
  inner_join(., 
             Gene_ID_18) %>% 
  dplyr::select(GeneID, 
                logFC, 
                adj.P.Val) %>% 
  rename(gene_ensembl = GeneID)%>% 
  inner_join(., 
             annotation_data) %>% 
  rename(GeneID = gene_name)

Trans_geo_hyb_18_snps = inner_join(Trans_geo_hyb_18, 
           ecotype_snps_fixed, 
           by = 'GeneID', 
           relationship = 'many-to-many')%>% 
  dplyr::select(gene_ensembl, 
                logFC, 
                adj.P.Val, 
                CHR, 
                GeneID,
                feature.x, 
                start.x, 
                end.x, 
                BP.x, 
                REF, 
                ALT, 
                other_affected_genes, 
                ecotype)

Trans_geo_hyb_18_snps %>% 
  # group_by(GeneID, 
  #          ecotype) %>% 
  group_by(GeneID) %>% 
  summarize(n = n()) 




