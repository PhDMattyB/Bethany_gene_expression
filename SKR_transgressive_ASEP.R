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
         effect,
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

# Trans_amb_hyb_12$BP = as.character(Trans_amb_hyb_12$BP)
# 
# ## combine with the 12 degree SNPS with allele specific expression
# 
# 
# 
# # ecotype_snps$BP = as.character(ecotype_snps$BP)
# annotation_data$BP = as.character(annotation_data$BP)

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
         effect,
         ecotype)

# Trans_amb_hyb_12_snps %>% 
#   write_csv('Trans_amb_hyb_12_TRANSGRESSIVE_EXP_snps.csv')


Trans_amb_hyb_12_snps %>% 
  group_by(GeneID,
           ecotype) %>%
  # group_by(GeneID) %>% 
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
                effect,
                ecotype)

# trans_geo_hyb_12_snps %>%
#   write_csv('Trans_geo_hyb_12_TRANSGRESSIVE_EXP_snps.csv')


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
                effect,
                ecotype)

# Trans_amb_hyb_18_snps %>%
#   write_csv('Trans_amb_hyb_18_TRANSGRESSIVE_EXP_snps.csv')


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
                effect,
                ecotype)

# Trans_geo_hyb_18_snps %>%
#   write_csv('Trans_geo_hyb_18_TRANSGRESSIVE_EXP_snps.csv')


Trans_geo_hyb_18_snps %>% 
  # group_by(GeneID, 
  #          ecotype) %>% 
  group_by(GeneID) %>% 
  summarize(n = n()) 



# Compare ASEP between groups ---------------------------------------------

Trans_amb_hyb_12_snps = read_csv('Trans_amb_hyb_12_TRANSGRESSIVE_EXP_snps.csv') %>% 
  mutate(temp = '12')
Trans_amb_hyb_18_snps = read_csv('Trans_amb_hyb_18_TRANSGRESSIVE_EXP_snps.csv') %>% 
  mutate(temp = '18')
Trans_geo_hyb_12_snps = read_csv('Trans_geo_hyb_12_TRANSGRESSIVE_EXP_snps.csv') %>% 
  mutate(temp = '12')
Trans_geo_hyb_18_snps = read_csv('Trans_geo_hyb_18_TRANSGRESSIVE_EXP_snps.csv') %>% 
  mutate(temp = '18')



# Ambient plastic ASEP ----------------------------------------------------
## ambient plasticit ASEP
amb_plast_snps = inner_join(Trans_amb_hyb_12_snps, 
                            Trans_amb_hyb_18_snps, 
                            by = c('gene_ensembl', 
                                   'CHR', 
                                   'GeneID', 
                                   'start.x', 
                                   'end.x', 
                                   'BP.x', 
                                   'REF', 
                                   'ALT', 
                                   'ecotype'),
                            relationship = 'many-to-many') 

amb_plast_snps %>% 
  group_by(GeneID, 
           ecotype) %>% 
  summarize(n = n()) 

## trangressive expressed snps Unique to 12 degress
anti_join(Trans_amb_hyb_12_snps, 
           Trans_amb_hyb_18_snps, 
           by = c('gene_ensembl', 
                  'CHR', 
                  'GeneID', 
                  'start.x', 
                  'end.x', 
                  'BP.x', 
                  'REF', 
                  'ALT', 
                  'ecotype')) 

## trasngresssive expressed snps unique to 18 degrees
anti_join(Trans_amb_hyb_18_snps, 
          Trans_amb_hyb_12_snps, 
          by = c('gene_ensembl', 
                 'CHR', 
                 'GeneID', 
                 'start.x', 
                 'end.x', 
                 'BP.x', 
                 'REF', 
                 'ALT', 
                 'ecotype')) 


# Geothermal plastic ASEP -------------------------------------------------

geo_plast_snps = inner_join(Trans_geo_hyb_12_snps, 
                            Trans_geo_hyb_18_snps, 
                            by = c('gene_ensembl', 
                                   'CHR', 
                                   'GeneID', 
                                   'start.x', 
                                   'end.x', 
                                   'BP.x', 
                                   'REF', 
                                   'ALT', 
                                   'ecotype'),
                            relationship = 'many-to-many') 

geo_plast_snps %>% 
  group_by(GeneID, 
           ecotype) %>% 
  summarize(n = n()) 

## trangressive expressed snps Unique to 12 degress
anti_join(Trans_geo_hyb_12_snps, 
          Trans_geo_hyb_18_snps, 
          by = c('gene_ensembl', 
                 'CHR', 
                 'GeneID', 
                 'start.x', 
                 'end.x', 
                 'BP.x', 
                 'REF', 
                 'ALT', 
                 'ecotype')) 

## trasngresssive expressed snps unique to 18 degrees
anti_join(Trans_geo_hyb_18_snps, 
          Trans_geo_hyb_12_snps, 
          by = c('gene_ensembl', 
                 'CHR', 
                 'GeneID', 
                 'start.x', 
                 'end.x', 
                 'BP.x', 
                 'REF', 
                 'ALT', 
                 'ecotype')) 


# divergence at 12 degrees ------------------------------------------------

div_12_snps = inner_join(Trans_amb_hyb_12_snps, 
                            Trans_geo_hyb_12_snps, 
                            by = c('gene_ensembl', 
                                   'CHR', 
                                   'GeneID', 
                                   'start.x', 
                                   'end.x', 
                                   'BP.x', 
                                   'REF', 
                                   'ALT', 
                                   'ecotype'),
                            relationship = 'many-to-many') 

div_12_snps %>% 
  group_by(GeneID, 
           ecotype) %>% 
  summarize(n = n()) 

## trangressive expressed snps Unique to ambient 
anti_join(Trans_amb_hyb_12_snps, 
          Trans_geo_hyb_12_snps, 
          by = c('gene_ensembl', 
                 'CHR', 
                 'GeneID', 
                 'start.x', 
                 'end.x', 
                 'BP.x', 
                 'REF', 
                 'ALT', 
                 'ecotype')) 
# %>% 
#   dplyr::select(GeneID, 
#                 BP.x, 
#                 REF, 
#                 ALT)

## trasngresssive expressed snps unique to geo
anti_join(Trans_geo_hyb_12_snps, 
          Trans_amb_hyb_12_snps, 
          by = c('gene_ensembl', 
                 'CHR', 
                 'GeneID', 
                 'start.x', 
                 'end.x', 
                 'BP.x', 
                 'REF', 
                 'ALT', 
                 'ecotype')) 
# %>% 
#   dplyr::select(GeneID, 
#                 BP.x, 
#                 REF, 
#                 ALT)



# divergence at 18 degrees ------------------------------------------------

div_18_snps = inner_join(Trans_amb_hyb_18_snps, 
                         Trans_geo_hyb_18_snps, 
                         by = c('gene_ensembl', 
                                'CHR', 
                                'GeneID', 
                                'start.x', 
                                'end.x', 
                                'BP.x', 
                                'REF', 
                                'ALT', 
                                'ecotype'),
                         relationship = 'many-to-many') 

div_18_snps %>% 
  group_by(GeneID, 
           ecotype) %>% 
  summarize(n = n()) 

## trangressive expressed snps Unique to ambient 
anti_join(Trans_amb_hyb_18_snps, 
          Trans_geo_hyb_18_snps, 
          by = c('gene_ensembl', 
                 'CHR', 
                 'GeneID', 
                 'start.x', 
                 'end.x', 
                 'BP.x', 
                 'REF', 
                 'ALT', 
                 'ecotype')) 
# %>% 
#   dplyr::select(GeneID, 
#                 BP.x, 
#                 REF, 
#                 ALT)

## trasngresssive expressed snps unique to geo
anti_join(Trans_geo_hyb_18_snps, 
          Trans_amb_hyb_18_snps, 
          by = c('gene_ensembl', 
                 'CHR', 
                 'GeneID', 
                 'start.x', 
                 'end.x', 
                 'BP.x', 
                 'REF', 
                 'ALT', 
                 'ecotype')) 
# %>% 
#   dplyr::select(GeneID, 
#                 BP.x, 
#                 REF, 
#                 ALT)



# BIG OVERLAP -------------------------------------------------------------

inner_join(amb_plast_snps, 
           geo_plast_snps, 
           by = c('gene_ensembl', 
                  'CHR', 
                  'GeneID', 
                  'start.x', 
                  'end.x', 
                  'BP.x', 
                  'REF', 
                  'ALT', 
                  'ecotype'))  
inner_join(div_12_snps, 
           div_18_snps, 
           by = c('gene_ensembl', 
                  'CHR', 
                  'GeneID', 
                  'start.x', 
                  'end.x', 
                  'BP.x', 
                  'REF', 
                  'ALT', 
                  'ecotype'))


# plot expression changes per ecotype!!!! ---------------------------------

theme_set(theme_bw())

Trans_amb_hyb_12_snps_plot = Trans_amb_hyb_12_snps %>% 
  dplyr::select(GeneID, 
                logFC, 
                adj.P.Val, 
                CHR, 
                BP.x, 
                REF, 
                ALT, 
                ecotype, 
                effect,
                temp) %>% 
  mutate(ecotemp = 'Ambient 12 degrees')
Trans_amb_hyb_18_snps_plot = Trans_amb_hyb_18_snps%>% 
  dplyr::select(GeneID, 
                logFC, 
                adj.P.Val, 
                CHR, 
                BP.x, 
                REF, 
                ALT, 
                ecotype, 
                effect,
                temp) %>% 
  mutate(ecotemp = 'Ambient 18 degrees')


Trans_amb_hyb_12_snps_plot = Trans_amb_hyb_12_snps %>% 
  filter(GeneID %in% Trans_amb_hyb_18_snps$GeneID) %>% 
  arrange(GeneID) %>% 
  rowid_to_column() %>% 
  rename(allele_id = rowid) %>% 
  unite(col = Gene_allele_id, 
        c(GeneID, 
          allele_id), 
        sep = '-', 
        remove = F) %>% 
  unite(col = ensembl_allele_id, 
        c(gene_ensembl, 
          allele_id), 
        sep = '-', 
        remove = F)

Trans_amb_hyb_18_snps_plot = Trans_amb_hyb_18_snps %>% 
  filter(GeneID %in% Trans_amb_hyb_12_snps$GeneID) %>% 
  arrange(GeneID) %>% 
  rowid_to_column() %>% 
  rename(allele_id = rowid)%>% 
  unite(col = Gene_allele_id, 
        c(GeneID, 
          allele_id), 
        sep = '-', 
        remove = F)%>% 
  unite(col = ensembl_allele_id, 
        c(gene_ensembl, 
          allele_id), 
        sep = '-', 
        remove = F)

trans_amb_hyb =  bind_rows(Trans_amb_hyb_12_snps_plot, 
                               Trans_amb_hyb_18_snps_plot) %>% 
   arrange(GeneID, 
           CHR, 
           BP.x)  

trans_exp_plot = c('#457b9d',
                   '#c1121f')

trans_amb_hyp_allele_plast = trans_amb_hyb %>% 
  ggplot(aes(x = temp, 
             y = logFC, 
             group = ecotype,
             col = ecotype))+
  geom_line(col = 'black')+
  geom_point() +
  scale_color_manual(values = trans_exp_plot)+
  facet_grid(~Gene_allele_id)+
  theme(panel.grid = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.position = 'none', 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'))


ggsave('Amb_hyb_12_transgressive_ASEP.tiff', 
       plot = trans_amb_hyp_allele_plast, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20, 
       height = 10)

